#!/usr/bin/env python3
"""Download and reformat InterProScan 6 member databases.

Extracts everything under ./interproscan/ so that each database follows the
path IPS6 FIND_MISSING_DATA expects:
    interproscan/{db_dir}/{version}/{db_subdir}

Run from the Nextflow task work-directory; publishDir copies the result.
"""
import urllib.request, json, tarfile, os, hashlib, sys, re, shutil

FTP         = "https://ftp.ebi.ac.uk/pub/software/unix/iprscan/6"
IPRSCAN_VER = "6.0"
OUTDIR      = "interproscan"


def log(msg):
    print(msg, flush=True)


def md5sum(path):
    h = hashlib.md5()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def download(url, dest):
    log("  GET " + url)
    urllib.request.urlretrieve(url, dest)


def norm(s):
    """Normalise a key the same way FIND_MISSING_DATA does in the IPS6 subworkflow."""
    return re.sub(r"[\s\-]+", "", s).lower()


def verify_and_fix(arcname, version, db_subdir):
    """Ensure the extracted layout matches what IPS6 FIND_MISSING_DATA expects.

    IPS6 expects: OUTDIR/{arcname}/{version}/{db_subdir}
    If the tarball extracted without the version layer, this function creates it.
    """
    ver_dir  = os.path.join(OUTDIR, arcname, version)
    expected = os.path.join(ver_dir, db_subdir) if db_subdir else ver_dir
    if os.path.exists(expected):
        return
    arc_root = os.path.join(OUTDIR, arcname)
    if not os.path.isdir(arc_root):
        log("  WARNING: " + arc_root + " not found after extraction")
        return
    if not os.path.isdir(ver_dir):
        os.makedirs(ver_dir, exist_ok=True)
        for item in os.listdir(arc_root):
            if item != version:
                shutil.move(os.path.join(arc_root, item), os.path.join(ver_dir, item))
        log("  Restructured: " + arcname + "/" + version + "/")
    if db_subdir and not os.path.exists(expected):
        log("  WARNING: expected path " + expected + " still missing")


def download_and_extract(arcname, version, db_subdir=""):
    arc   = arcname + "-" + version + ".tar.gz"
    md5_f = arc + ".md5"
    url   = FTP + "/" + IPRSCAN_VER + "/" + arcname + "/" + arc
    download(url, arc)
    download(url + ".md5", md5_f)
    with open(md5_f) as fh:
        expected_md5 = fh.read().split()[0]
    if md5sum(arc) != expected_md5:
        raise ValueError("MD5 mismatch for " + arc)
    with tarfile.open(arc) as t:
        t.extractall(OUTDIR)
    os.remove(arc)
    os.remove(md5_f)
    verify_and_fix(arcname, version, db_subdir)


# Step 1: resolve latest compatible InterPro version
log("Fetching compatible InterPro versions...")
with urllib.request.urlopen(FTP + "/" + IPRSCAN_VER + "/versions.json") as r:
    versions = json.loads(r.read())
ipro_ver = str(versions["interpro"][-1])
log("Using InterPro version: " + ipro_ver)

# Step 2: interpro metadata archive
log("Downloading InterPro metadata archive...")
os.makedirs(OUTDIR, exist_ok=True)
if not os.path.isdir(OUTDIR + "/interpro/" + ipro_ver):
    download_and_extract("interpro", ipro_ver)
else:
    log("  Already present: interpro/" + ipro_ver)

# Step 3: read databases.json; normalise keys like FIND_MISSING_DATA does
dbs_path = OUTDIR + "/interpro/" + ipro_ver + "/databases.json"
with open(dbs_path) as fh:
    raw_dbs = json.load(fh)
norm_dbs = {norm(k): v for k, v in raw_dbs.items()}
log("databases.json normalised keys (first 8): " + str(sorted(norm_dbs)[:8]))

# Step 4: build download list
# arcname    = FTP directory name / tarball prefix
# norm_key   = normalised databases.json key (matches FIND_MISSING_DATA logic)
# db_subdir  = second path component from applications.config dir field
# CATH-Gene3D and CATH-FunFam share one tarball ("cath")
db_list = [
    ("AntiFam",     "antifam",     "antifam",         ""),
    ("CATH",        "cath",        "cathgene3d",       "gene3d"),
    ("CDD",         "cdd",         "cdd",             ""),
    ("HAMAP",       "hamap",       "hamap",           ""),
    ("NCBIFAM",     "ncbifam",     "ncbifam",         ""),
    ("PANTHER",     "panther",     "panther",         ""),
    ("Pfam",        "pfam",        "pfam",            ""),
    ("PIRSF",       "pirsf",       "pirsf",           ""),
    ("PIRSR",       "pirsr",       "pirsr",           ""),
    ("PRINTS",      "prints",      "prints",          ""),
    ("PROSITE",     "prosite",     "prositeprofiles",  ""),
    ("SFLD",        "sfld",        "sfld",            ""),
    ("SMART",       "smart",       "smart",           ""),
    ("SUPERFAMILY", "superfamily", "superfamily",     ""),
]

errors     = []
downloaded = set()
for label, arcname, norm_key, db_subdir in db_list:
    version = norm_dbs.get(norm_key)
    if version is None:
        log("Skipping " + label + ": key " + norm_key + " not in databases.json")
        continue
    if db_subdir:
        expected_path = os.path.join(OUTDIR, arcname, version, db_subdir)
    else:
        expected_path = os.path.join(OUTDIR, arcname, version)
    if os.path.exists(expected_path):
        log("Skipping " + label + " " + version + ": already at " + expected_path)
        continue
    if arcname in downloaded:
        verify_and_fix(arcname, version, db_subdir)
        continue
    log("Downloading " + label + " " + version + "...")
    try:
        download_and_extract(arcname, version, db_subdir)
        downloaded.add(arcname)
        log("  Done: " + label + " " + version)
    except Exception as e:
        msg = "WARNING: Failed to download " + label + " " + version + ": " + str(e)
        log(msg)
        errors.append(msg)

if errors:
    print("\nDatabases that could not be downloaded:", file=sys.stderr)
    for e in errors:
        print("  " + e, file=sys.stderr)

log("InterProScan database setup complete.")
