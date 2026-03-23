#!/usr/bin/env python3
"""
Main download script: bulk RNA-seq, ChIP-seq, CUT&RUN, and external datasets.

Run this on any machine (no Cell Ranger needed). Uses nf-core pipelines for processing.
Data: Alexanian et al. Nature 2024 | GEO GSE221699 | BioProject PRJNA915384

Disk space: ~150 GB (38 runs)

Includes:
  - Bulk RNA-seq (iPSC-CF): Unstim, IL1B, TGFB, TGFB+IL1B → Fig 5
  - ChIP-seq (iPSC-CF): H3K27Ac, p65/RELA, input → Fig 5
  - CUT&RUN: BRD4 in Cx3cr1+ Sham/TAC → Fig 4
  - Sorted fibroblasts bulk RNA-seq (GSE247261)
  - Optional: Kuppe et al. human cardiac scATAC (GSE183852) for Fig 4N

Run from project root: python scripts/download_data.py
"""

import argparse
import hashlib
import os
import subprocess
import sys
from pathlib import Path

# ENA run metadata for non-Cell-Ranger data (Bulk RNA-seq, ChIP-seq, CUT&RUN, Sorted FBs)
OTHER_RUNS = {
    # Bulk RNA-seq (16 runs)
    "SRR22882084": ("Bulk RNA-Seq Unstim-1", "paired"),
    "SRR22882085": ("Bulk RNA-Seq Unstim-2", "paired"),
    "SRR22882086": ("Bulk RNA-Seq Unstim-3", "paired"),
    "SRR22882087": ("Bulk RNA-Seq Unstim-4", "paired"),
    "SRR22882088": ("Bulk RNA-Seq IL1B-1", "paired"),
    "SRR22882089": ("Bulk RNA-Seq IL1B-2", "paired"),
    "SRR22882090": ("Bulk RNA-Seq IL1B-3", "paired"),
    "SRR22882091": ("Bulk RNA-Seq IL1B-4", "paired"),
    "SRR22882092": ("Bulk RNA-Seq TGFB-1", "paired"),
    "SRR22882093": ("Bulk RNA-Seq TGFB-2", "paired"),
    "SRR22882094": ("Bulk RNA-Seq TGFB-3", "paired"),
    "SRR22882095": ("Bulk RNA-Seq TGFB-4", "paired"),
    "SRR22882096": ("Bulk RNA-Seq TGFB+IL1B-1", "paired"),
    "SRR22882097": ("Bulk RNA-Seq TGFB+IL1B-2", "paired"),
    "SRR22882098": ("Bulk RNA-Seq TGFB+IL1B-3", "paired"),
    "SRR22882099": ("Bulk RNA-Seq TGFB+IL1B-4", "paired"),
    # ChIP-seq (6 runs)
    "SRR22882100": ("ChIP-Seq Uns_iPS-CF_input_S1_R1_001", "single"),
    "SRR22882101": ("ChIP-Seq Uns_iPS-CF_K27Ac_A_S4_R1_001", "single"),
    "SRR22882102": ("ChIP-Seq Uns_iPS-CF_K27Ac_B_S5_R1_001", "single"),
    "SRR22882103": ("ChIP-Seq Uns_iPS-CF_RELA_S10_R1_001", "single"),
    "SRR22882104": ("ChIP-Seq Inp_Uns_S1_R1_001_Nov2022", "single"),
    "SRR22882105": ("ChIP-Seq Uns_AM_RelA_S11_R1_001_Nov2022", "single"),
    # CUT&RUN (4 runs)
    "SRR22882106": ("Cut and Run Sham_CX3CR1pos_FLAG", "paired"),
    "SRR22882107": ("Cut and Run TAC1_CX3CR1pos_FLAG", "paired"),
    "SRR22882108": ("Cut and Run TAC2_CX3CR1pos_FLAG", "paired"),
    "SRR22882109": ("Cut and Run TAC2_CX3CR1neg_IgG", "paired"),
    # Sorted fibroblasts bulk RNA-seq (GSE247261)
    "SRR26716616": ("Sorted FBs TAC Il1b-KO_Rep3", "paired"),
    "SRR26716617": ("Sorted FBs TAC Il1b-KO_Rep3", "paired"),
    "SRR26716618": ("Sorted FBs TAC WT _Rep3", "paired"),
    "SRR26716619": ("Sorted FBs TAC WT _Rep3", "paired"),
    "SRR26716620": ("Sorted FBs TAC Il1b-KO_Rep2", "paired"),
    "SRR26716621": ("Sorted FBs TAC Il1b-KO_Rep2", "paired"),
    "SRR26716622": ("Sorted FBs TAC Il1b-KO_Rep1", "paired"),
    "SRR26716623": ("Sorted FBs TAC Il1b-KO_Rep1", "paired"),
    "SRR26716624": ("Sorted FBs TAC WT _Rep2", "paired"),
    "SRR26716625": ("Sorted FBs TAC WT _Rep2", "paired"),
    "SRR26716626": ("Sorted FBs TAC WT _Rep1", "paired"),
    "SRR26716627": ("Sorted FBs TAC WT _Rep1", "paired"),
}

# External datasets (optional)
EXTERNAL = {
    "GSE183852": {
        "name": "Kuppe et al. human cardiac scATAC-seq",
        "url": "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE183852&format=file",
        "note": "Fig 4N - IL1B locus in control vs MI. Large; check GEO for processed data.",
    },
    # Link et al. RELA ChIP (BMDM) - accession TBD, add when found
}


def download_with_curl(url: str, out_path: Path, expected_md5: str | None = None) -> bool:
    """Download file using curl. Uses 2h timeout and resume (-C -) for large files."""
    out_path.parent.mkdir(parents=True, exist_ok=True)
    cmd = [
        "curl", "-f", "-L", "-C", "-",  # -C - enables resume for partial downloads
        "--max-time", "7200",             # 2 hours per file
        "--connect-timeout", "120",       # 2 min to establish connection
        "-o", str(out_path), url,
    ]
    try:
        subprocess.run(cmd, check=True, capture_output=True)
    except subprocess.CalledProcessError as e:
        print(f"  Error: {e.stderr.decode() if e.stderr else e}", file=sys.stderr)
        if out_path.exists():
            out_path.unlink()  # Remove partial file so re-run will retry
        return False

    if expected_md5:
        with open(out_path, "rb") as f:
            digest = hashlib.md5(f.read()).hexdigest()
        if digest != expected_md5:
            print(f"  MD5 mismatch: got {digest}, expected {expected_md5}", file=sys.stderr)
            return False
    return True


def download_run_ena(run_id: str, out_dir: Path, md5_checks: bool = False) -> bool:
    """Fetch FASTQ URLs from ENA and download."""
    import urllib.request

    url = f"https://www.ebi.ac.uk/ena/portal/api/filereport?accession={run_id}&result=read_run&fields=fastq_ftp,fastq_md5"
    req = urllib.request.Request(url, headers={"User-Agent": "Chromatin-HF-Download/1.0"})
    with urllib.request.urlopen(req) as resp:
        text = resp.read().decode()
    lines = text.strip().split("\n")
    if len(lines) < 2:
        print(f"  No ENA data for {run_id}", file=sys.stderr)
        return False

    headers = lines[0].split("\t")
    vals = lines[1].split("\t")
    data = dict(zip(headers, vals))
    fastq_ftp = data.get("fastq_ftp", "")
    fastq_md5 = data.get("fastq_md5", "")

    if not fastq_ftp:
        print(f"  No FASTQ FTP for {run_id}", file=sys.stderr)
        return False

    urls = []
    for u in fastq_ftp.split(";"):
        u = u.strip()
        if not u:
            continue
        if u.startswith("http"):
            urls.append(u)
        else:
            urls.append(f"ftp://{u}")
    md5_list = (fastq_md5.split(";") if fastq_md5 else [])
    md5s = [md5_list[i] if i < len(md5_list) else None for i in range(len(urls))]

    run_dir = out_dir / run_id
    run_dir.mkdir(parents=True, exist_ok=True)
    ok = True
    for i, (u, m) in enumerate(zip(urls, md5s)):
        fname = Path(u).name
        out_path = run_dir / fname
        if out_path.exists() and (not m or hashlib.md5(out_path.read_bytes()).hexdigest() == m):
            print(f"  Skip (exists): {fname}")
            continue
        print(f"  Downloading {fname} ...")
        if not download_with_curl(u, out_path, m if md5_checks else None):
            ok = False
    return ok


def download_with_fasterq_dump(run_id: str, out_dir: Path) -> bool:
    """Use SRA fasterq-dump (requires sra-tools)."""
    run_dir = out_dir / run_id
    run_dir.mkdir(parents=True, exist_ok=True)
    cmd = [
        "fasterq-dump",
        run_id,
        "-O", str(run_dir),
        "-f",
        "-e", str(os.cpu_count() or 4),
    ]
    try:
        subprocess.run(cmd, check=True)
        return True
    except (subprocess.CalledProcessError, FileNotFoundError) as e:
        print(f"  fasterq-dump failed: {e}", file=sys.stderr)
        return False


def download_external(acc: str, out_dir: Path) -> bool:
    """Download external GEO dataset (placeholder - GEO package download may require auth)."""
    info = EXTERNAL.get(acc)
    if not info:
        print(f"  Unknown external: {acc}", file=sys.stderr)
        return False
    print(f"  {info['name']}")
    print(f"  Note: {info['note']}")
    print("  For GSE183852: use GEO webpage or gget/GEOfetch to download.")
    print("  Run: pip install gget && gget geo GSE183852 -o data/external/")
    return True


def main():
    parser = argparse.ArgumentParser(
        description="Download bulk RNA-seq, ChIP-seq, CUT&RUN, and optional external data"
    )
    parser.add_argument(
        "-o", "--output-dir",
        type=Path,
        default=Path("data"),
        help="Output directory (default: data)",
    )
    parser.add_argument(
        "--use-sra",
        action="store_true",
        help="Use fasterq-dump instead of ENA direct download (requires sra-tools)",
    )
    parser.add_argument(
        "--runs",
        nargs="+",
        default=None,
        help="Specific run IDs to download (default: all)",
    )
    parser.add_argument(
        "--md5",
        action="store_true",
        help="Verify MD5 checksums (slower)",
    )
    parser.add_argument(
        "--external",
        nargs="+",
        default=[],
        metavar="ACC",
        help="External GEO accessions (e.g. GSE183852). Prints instructions.",
    )
    args = parser.parse_args()

    runs = args.runs or list(OTHER_RUNS)
    for r in runs:
        if r not in OTHER_RUNS:
            print(f"Unknown run: {r}", file=sys.stderr)
            sys.exit(1)

    args.output_dir.mkdir(parents=True, exist_ok=True)
    print(f"Downloading {len(runs)} runs to {args.output_dir}")
    print("~150 GB total. Ensure sufficient disk space.\n")

    failed = []
    for i, run_id in enumerate(runs, 1):
        label = OTHER_RUNS[run_id][0]
        print(f"[{i}/{len(runs)}] {run_id} ({label})")
        if args.use_sra:
            ok = download_with_fasterq_dump(run_id, args.output_dir)
        else:
            ok = download_run_ena(run_id, args.output_dir, md5_checks=args.md5)
        if not ok:
            failed.append(run_id)

    if args.external:
        print("\n--- External datasets ---")
        ext_dir = args.output_dir / "external"
        ext_dir.mkdir(parents=True, exist_ok=True)
        for acc in args.external:
            download_external(acc, ext_dir)

    if failed:
        print(f"\nFailed: {failed}", file=sys.stderr)
        sys.exit(1)
    print("\nDone.")


if __name__ == "__main__":
    main()
