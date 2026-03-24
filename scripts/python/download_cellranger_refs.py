#!/usr/bin/env python3
"""
Download Cell Ranger reference genomes for mouse (mm10/GRCm39).

Run on x86 before scripts/python/run_cellranger.py. Saves to data/refs/ by default.

References needed for this project:
  - Mouse scRNA (GRCm39 2024-A): for cellranger count
  - Mouse scATAC (GRCm39 2024-A): for cellranger-atac count

Run from project root: python scripts/python/download_cellranger_refs.py
"""

import argparse
import subprocess
import sys
from pathlib import Path

# 10x Genomics reference URLs (mouse only - for Alexanian study)
REF_SCRNA = {
    "url": "https://cf.10xgenomics.com/supp/cell-exp/refdata-gex-GRCm39-2024-A.tar.gz",
    "dirname": "refdata-gex-GRCm39-2024-A",
}
REF_ATAC = {
    # Cell Ranger ATAC 2.2 uses cell-arc refs (from https://www.10xgenomics.com/support/software/cell-ranger-atac/downloads)
    "url": "https://cf.10xgenomics.com/supp/cell-arc/refdata-cellranger-arc-GRCm39-2024-A.tar.gz",
    "dirname": "refdata-cellranger-arc-GRCm39-2024-A",
    "manual_url": "https://www.10xgenomics.com/support/software/cell-ranger-atac/downloads",
}
REF_SIZES = {"scRNA": "9.7 GB", "scATAC": "~13 GB"}


def download_with_curl(url: str, out_path: Path) -> bool:
    """Download file using curl. Uses 2h timeout and resume (-C -) for large files."""
    out_path.parent.mkdir(parents=True, exist_ok=True)
    cmd = [
        "curl", "-f", "-L", "-C", "-",
        "-H", "User-Agent: Mozilla/5.0 (compatible; Chromatin-HF/1.0)",
        "--max-time", "7200",
        "--connect-timeout", "120",
        "-o", str(out_path), url,
    ]
    try:
        subprocess.run(cmd, check=True, capture_output=True)
        return True
    except subprocess.CalledProcessError as e:
        print(f"  Error: {e.stderr.decode() if e.stderr else e}", file=sys.stderr)
        if out_path.exists():
            out_path.unlink()
        return False


def extract_tar(tar_path: Path, out_dir: Path) -> bool:
    """Extract .tar.gz to out_dir (parent of archive)."""
    try:
        subprocess.run(
            ["tar", "-xzf", str(tar_path), "-C", str(out_dir)],
            check=True,
            capture_output=True,
        )
        return True
    except subprocess.CalledProcessError as e:
        print(f"  Error extracting: {e}", file=sys.stderr)
        return False


def download_ref(name: str, ref: dict, out_dir: Path, extract: bool = True) -> Path | None:
    """Download a reference tarball, optionally extract. Returns path to ref dir or None."""
    out_dir.mkdir(parents=True, exist_ok=True)
    tar_path = out_dir / Path(ref["url"]).name
    ref_dir = out_dir / ref["dirname"]

    if ref_dir.exists():
        print(f"  Skip (exists): {ref_dir}")
        return ref_dir

    if not tar_path.exists():
        print(f"  Downloading {name} ({REF_SIZES.get(name, '?')})...")
        if not download_with_curl(ref["url"], tar_path):
            return None
    else:
        print(f"  Using cached: {tar_path}")

    if extract:
        print(f"  Extracting...")
        if not extract_tar(tar_path, out_dir):
            return None
        # Optionally remove tarball to save space
        # tar_path.unlink()
    return ref_dir if ref_dir.exists() else None


def main():
    parser = argparse.ArgumentParser(
        description="Download Cell Ranger mouse references (scRNA + scATAC)"
    )
    parser.add_argument(
        "-o", "--output-dir",
        type=Path,
        default=Path("data/refs"),
        help="Output directory (default: data/refs)",
    )
    parser.add_argument(
        "--scrna-only",
        action="store_true",
        help="Download only scRNA reference",
    )
    parser.add_argument(
        "--atac-only",
        action="store_true",
        help="Download only scATAC reference",
    )
    parser.add_argument(
        "--no-extract",
        action="store_true",
        help="Download tarballs only, do not extract",
    )
    args = parser.parse_args()

    do_scrna = not args.atac_only
    do_atac = not args.scrna_only

    args.output_dir.mkdir(parents=True, exist_ok=True)
    print(f"Downloading Cell Ranger references to {args.output_dir}\n")

    paths = {}
    if do_scrna:
        p = download_ref("scRNA", REF_SCRNA, args.output_dir, extract=not args.no_extract)
        if p:
            paths["scrna"] = p
        else:
            print("  Failed to get scRNA reference", file=sys.stderr)
            sys.exit(1)

    if do_atac:
        p = download_ref("scATAC", REF_ATAC, args.output_dir, extract=not args.no_extract)
        if p:
            paths["atac"] = p
        else:
            print("\n  scATAC reference returned 403 (10x may require download from their site).")
            print("  Download manually from:")
            print(f"    {REF_ATAC['manual_url']}")
            print("  Look for 'Mouse reference (GRCm39) - 2024-A' under References.")
            print("  Extract to:", args.output_dir / REF_ATAC["dirname"])
            sys.exit(1)

    print("\nDone. Run Cell Ranger with:")
    if "scrna" in paths:
        print(f"  --ref-scrna {paths['scrna']}")
    if "atac" in paths:
        print(f"  --ref-atac {paths['atac']}")
    print("\nExample:")
    print(f"  python scripts/python/run_cellranger.py \\")
    print(f"    --ref-scrna {paths.get('scrna', '...')} \\")
    print(f"    --ref-atac {paths.get('atac', '...')}")


if __name__ == "__main__":
    main()
