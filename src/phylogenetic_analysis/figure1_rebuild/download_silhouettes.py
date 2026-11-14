#!/usr/bin/env python3
"""
Download representative bird-family silhouettes from PhyloPic (https://www.phylopic.org).

Usage
-----
    python3 download_silhouettes.py inno_db_25_clootl_overlap.csv  [-o silhouettes] [--min-tips 15]

Needs only the Python standard library (3.8+) and internet access.

What it does
------------
1. Reads the trait table (columns: Order, Family, Genus, Species) and picks the
   families to draw:
     * the largest family of every Order (so every Order is represented), plus
     * every family with at least --min-tips species.
   Edit FAMILIES_EXTRA / FAMILIES_SKIP below to add or drop families by hand.
2. For each family it asks the PhyloPic API for a silhouette, in this order:
     family node image -> any image in the family clade -> most frequent genus
     in your table -> most frequent species of that genus.
   Only images that may be reused in a paper are accepted: CC0, Public Domain
   Mark, CC-BY and CC-BY-SA. NonCommercial images are skipped.
3. Saves   <out>/<FAMILY>.svg   (vector, preferred)  and  <out>/<FAMILY>.png
   and writes  <out>/attribution.csv  with author, licence and source URL.
   CC-BY / CC-BY-SA silhouettes MUST be credited in the figure legend; copy
   the credits from attribution.csv.

Then zip / upload the output folder.
"""
import argparse
import csv
import json
import re
import sys
import time
import urllib.error
import urllib.parse
import urllib.request
from collections import Counter, defaultdict
from pathlib import Path

API = "https://api.phylopic.org"
UA = {"User-Agent": "silhouette-downloader/1.0 (figure preparation)"}

# Families that always get a silhouette, even if small (UPPER CASE as in the table).
FAMILIES_EXTRA = []
# Families you do not want.
FAMILIES_SKIP = []

# Licences allowed for reuse (substring match on the licence URL).
OK_LICENSES = ("/publicdomain/zero/", "/publicdomain/mark/", "/licenses/by/", "/licenses/by-sa/")


def get_json(url, tries=4):
    for k in range(tries):
        try:
            req = urllib.request.Request(url, headers={**UA, "Accept": "application/json"})
            with urllib.request.urlopen(req, timeout=30) as r:
                return json.load(r)
        except urllib.error.HTTPError as e:
            if e.code == 404:
                return None
            err = e
        except Exception as e:  # noqa: BLE001
            err = e
        time.sleep(1.5 * (k + 1))
    print(f"   ! request failed: {url} ({err})", file=sys.stderr)
    return None


def get_bytes(url, tries=4):
    for k in range(tries):
        try:
            with urllib.request.urlopen(urllib.request.Request(url, headers=UA), timeout=60) as r:
                return r.read()
        except Exception as e:  # noqa: BLE001
            err = e
            time.sleep(1.5 * (k + 1))
    print(f"   ! download failed: {url} ({err})", file=sys.stderr)
    return None


def abs_url(href):
    return href if href.startswith("http") else API + href


def license_ok(lic_url):
    return bool(lic_url) and any(s in lic_url for s in OK_LICENSES)


def image_record(img):
    """Pull the useful fields out of a PhyloPic image object."""
    links = img.get("_links", {})
    lic = links.get("license", {}).get("href", "")
    vec = links.get("vectorFile", {}).get("href")
    rasters = links.get("rasterFiles", [])
    png = None
    if rasters:  # largest raster
        def px(r):
            m = re.match(r"(\d+)x(\d+)", r.get("sizes", "0x0"))
            return int(m.group(1)) * int(m.group(2)) if m else 0
        png = max(rasters, key=px).get("href")
    return {
        "image_uuid": img.get("uuid", ""),
        "license": lic,
        "attribution": img.get("attribution") or "",
        "contributor": links.get("contributor", {}).get("title", ""),
        "svg": vec,
        "png": png,
        "page": f"https://www.phylopic.org/images/{img.get('uuid', '')}",
    }


def find_nodes(build, name):
    q = urllib.parse.urlencode({"build": build, "filter_name": name.lower(), "page": 0, "embed_items": "true"})
    d = get_json(f"{API}/nodes?{q}")
    if not d:
        return []
    items = d.get("_embedded", {}).get("items", [])
    # exact (case-insensitive) name match first
    def names(n):
        return [x.get("name", "").lower() for x in n.get("names", [])] or [""]
    exact = [n for n in items if name.lower() in [s.lower() for s in sum(
        ([part.get("text", "") for part in nm] if isinstance(nm, list) else [str(nm)] for nm in n.get("names", [])), [])]]
    return exact or items[:1]


def image_for_node(build, node):
    """1) the node's own primary image  2) any reusable image in its clade."""
    uuid = node.get("uuid")
    prim = node.get("_links", {}).get("primaryImage", {}).get("href")
    if prim:
        sep = "&" if "?" in prim else "?"
        img = get_json(abs_url(prim) + ("" if "build=" in prim else f"{sep}build={build}"))
        if img:
            rec = image_record(img)
            if license_ok(rec["license"]) and rec["svg"]:
                return rec
    if uuid:
        q = urllib.parse.urlencode({
            "build": build, "filter_clade": uuid, "filter_license_nc": "false",
            "page": 0, "embed_items": "true"})
        d = get_json(f"{API}/images?{q}")
        for img in (d or {}).get("_embedded", {}).get("items", []):
            rec = image_record(img)
            if license_ok(rec["license"]) and rec["svg"]:
                return rec
    return None


def pick_families(rows, min_tips):
    fam_n = Counter(r["Family"] for r in rows)
    by_order = defaultdict(Counter)
    for r in rows:
        by_order[r["Order"]][r["Family"]] += 1
    chosen = {c.most_common(1)[0][0] for c in by_order.values()}
    chosen |= {f for f, n in fam_n.items() if n >= min_tips}
    chosen |= {f.upper() for f in FAMILIES_EXTRA}
    chosen -= {f.upper() for f in FAMILIES_SKIP}
    return sorted(chosen, key=lambda f: -fam_n[f]), fam_n


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("table", help="CSV with Order, Family, Genus, Species columns")
    ap.add_argument("-o", "--out", default="silhouettes")
    ap.add_argument("--min-tips", type=int, default=15)
    a = ap.parse_args()

    with open(a.table, newline="", encoding="utf-8") as fh:
        rows = [{k: (v or "").strip() for k, v in r.items()} for r in csv.DictReader(fh)]
    for r in rows:
        r["Family"], r["Order"] = r["Family"].upper(), r["Order"].upper()
    fams, fam_n = pick_families(rows, a.min_tips)
    print(f"{len(fams)} families selected")

    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)

    root = get_json(f"{API}/")
    if not root:
        sys.exit("Cannot reach api.phylopic.org - check the network / allowed domains.")
    build = root["build"]
    print("PhyloPic build", build)

    report, missing = [], []
    for fam in fams:
        genus_n = Counter(r["Genus"] for r in rows if r["Family"] == fam and r["Genus"])
        sp_n = Counter(r["Species"].replace("_", " ") for r in rows if r["Family"] == fam)
        candidates = [("family", fam.title())]
        candidates += [("genus", g) for g, _ in genus_n.most_common(3)]
        candidates += [("species", s) for s, _ in sp_n.most_common(3)]
        rec, used = None, None
        for level, name in candidates:
            for node in find_nodes(build, name):
                rec = image_for_node(build, node)
                if rec:
                    used = (level, name)
                    break
            if rec:
                break
        if not rec:
            print(f"  [--] {fam:<22} no reusable silhouette found")
            missing.append(fam)
            continue
        ok = False
        if rec["svg"]:
            data = get_bytes(rec["svg"])
            if data:
                (out / f"{fam}.svg").write_bytes(data)
                ok = True
        if rec["png"]:
            data = get_bytes(rec["png"])
            if data:
                (out / f"{fam}.png").write_bytes(data)
                ok = True
        if ok:
            print(f"  [ok] {fam:<22} <- {used[0]}: {used[1]}  ({rec['license'].split('/')[-3:-1] and rec['license']})")
            report.append({"Family": fam, "n_species_in_table": fam_n[fam],
                           "matched_level": used[0], "matched_name": used[1],
                           "contributor": rec["contributor"], "attribution": rec["attribution"],
                           "license": rec["license"], "source": rec["page"]})
        else:
            missing.append(fam)
        time.sleep(0.3)  # be polite to the API

    with open(out / "attribution.csv", "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=["Family", "n_species_in_table", "matched_level", "matched_name",
                                           "contributor", "attribution", "license", "source"])
        w.writeheader()
        w.writerows(report)
    print(f"\nDone: {len(report)} silhouettes in {out}/   missing: {len(missing)}")
    if missing:
        print("Missing (draw/pick by hand and save as <FAMILY>.svg):", ", ".join(missing))


if __name__ == "__main__":
    main()
