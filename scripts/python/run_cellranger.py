#!/usr/bin/env python3
"""
Run Cell Ranger on downloaded scRNA-seq and scATAC-seq FASTQs.

Run on x86 where Cell Ranger is installed. Creates 10x-compatible symlinks
(ENA naming → Illumina-style) and runs cellranger count / cellranger-atac count.

Prerequisites:
  - scripts/python/download_cellranger_data.py has been run
  - Cell Ranger and Cell Ranger ATAC installed
  - mm10 reference genomes downloaded (see INSTALL.md or 10x Genomics support)

Run from project root:
  python scripts/python/run_cellranger.py --ref-scrna /path/to/refdata-gex-mm10-2020-A \\
                                  --ref-atac /path/to/refdata-cellranger-arc-GRCm39-2024-A
"""

import argparse
import subprocess
import sys
from pathlib import Path

# Ensure scripts/python/ is on path for import
_scripts_dir = Path(__file__).resolve().parent
if str(_scripts_dir) not in sys.path:
    sys.path.insert(0, str(_scripts_dir))

try:
    from download_cellranger_data import CELLRANGER_RUNS, SCATAC_RUNS
except ImportError:
    sys.exit("Could not import download_cellranger_data. Run from project root: python scripts/python/run_cellranger.py")


def create_symlinks(run_dir: Path, run_id: str) -> bool:
    """
    Create 10x-style symlinks for ENA FASTQs.
    ENA: SRRxxxxx_1.fastq.gz, SRRxxxxx_2.fastq.gz
    10x: SRRxxxxx_S1_L001_R1_001.fastq.gz, SRRxxxxx_S1_L001_R2_001.fastq.gz
    """
    fq1 = run_dir / f"{run_id}_1.fastq.gz"
    fq2 = run_dir / f"{run_id}_2.fastq.gz"
    if not fq1.exists() or not fq2.exists():
        return False

    ln1 = run_dir / f"{run_id}_S1_L001_R1_001.fastq.gz"
    ln2 = run_dir / f"{run_id}_S1_L001_R2_001.fastq.gz"
    for src, dst in [(fq1, ln1), (fq2, ln2)]:
        if not dst.exists():
            dst.symlink_to(src.name)
    return True


def create_symlinks_atac(run_dir: Path, run_id: str) -> bool:
    """
    Create 10x-style symlinks for scATAC FASTQs from SRA (--split-files --include-technical).
    SRA: SRRxxxxx_1.fastq, SRRxxxxx_2.fastq, SRRxxxxx_3.fastq
    10x Chromium: R1, I2 (16bp barcode), R2. Accepts .fastq or .fastq.gz.
    """
    ext = None
    for e in (".fastq.gz", ".fastq"):
        if (run_dir / f"{run_id}_1{e}").exists():
            ext = e
            break
    if not ext:
        return False
    fq1 = run_dir / f"{run_id}_1{ext}"
    fq2 = run_dir / f"{run_id}_2{ext}"
    fq3 = run_dir / f"{run_id}_3{ext}"
    if not fq1.exists() or not fq2.exists() or not fq3.exists():
        return False
    # 10x Chromium scATAC: R1, I2 (barcode), R2
    ln1 = run_dir / f"{run_id}_S1_L001_R1_001{ext}"
    ln2 = run_dir / f"{run_id}_S1_L001_I2_001{ext}"
    ln3 = run_dir / f"{run_id}_S1_L001_R2_001{ext}"
    for src, dst in [(fq1, ln1), (fq2, ln2), (fq3, ln3)]:
        if not dst.exists():
            dst.symlink_to(src.name)
    return True


def run_cellranger_count(
    run_id: str,
    fastq_dir: Path,
    out_dir: Path,
    ref: Path,
    include_introns: bool = True,
    localcores: int | None = None,
    localmem: int | None = None,
    dry_run: bool = False,
) -> bool:
    """Run cellranger count for scRNA/snRNA-seq."""
    out_path = out_dir / run_id
    if (out_path / "outs" / "filtered_feature_bc_matrix.h5").exists():
        print(f"  Skip (exists): {run_id}")
        return True

    cmd = [
        "cellranger", "count",
        "--id", run_id,
        "--transcriptome", str(ref),
        "--fastqs", str(fastq_dir),
        "--sample", run_id,
    ]
    # Cell Ranger v10 expects an explicit true/false value.
    # Older versions accepted a bare flag; passing the explicit value works reliably.
    cmd.extend(["--include-introns", "true" if include_introns else "false"])
    if localcores:
        cmd.extend(["--localcores", str(localcores)])
    if localmem:
        cmd.extend(["--localmem", str(localmem)])

    if dry_run:
        print(f"  [dry-run] {' '.join(cmd)}")
        return True

    try:
        subprocess.run(cmd, cwd=out_dir, check=True)
        return True
    except subprocess.CalledProcessError as e:
        print(f"  Error: cellranger count failed for {run_id}", file=sys.stderr)
        return False
    except FileNotFoundError:
        print("  Error: cellranger not found. Add to PATH or activate environment.", file=sys.stderr)
        raise


def run_cellranger_atac_count(
    run_id: str,
    fastq_dir: Path,
    out_dir: Path,
    ref: Path,
    localcores: int | None = None,
    localmem: int | None = None,
    dry_run: bool = False,
) -> bool:
    """Run cellranger-atac count for scATAC-seq."""
    out_path = out_dir / run_id
    if (out_path / "outs" / "fragments.tsv.gz").exists():
        print(f"  Skip (exists): {run_id}")
        return True

    cmd = [
        "cellranger-atac", "count",
        "--id", run_id,
        "--reference", str(ref),
        "--fastqs", str(fastq_dir),
        "--sample", run_id,
    ]
    if localcores:
        cmd.extend(["--localcores", str(localcores)])
    if localmem:
        cmd.extend(["--localmem", str(localmem)])

    if dry_run:
        print(f"  [dry-run] {' '.join(cmd)}")
        return True

    try:
        subprocess.run(cmd, cwd=out_dir, check=True)
        return True
    except subprocess.CalledProcessError as e:
        print(f"  Error: cellranger-atac count failed for {run_id}", file=sys.stderr)
        return False
    except FileNotFoundError:
        print("  Error: cellranger-atac not found. Add to PATH.", file=sys.stderr)
        raise


def main():
    parser = argparse.ArgumentParser(
        description="Run Cell Ranger on downloaded scRNA/scATAC FASTQs (x86)"
    )
    parser.add_argument(
        "-i", "--fastq-dir",
        type=Path,
        default=Path("data/cellranger"),
        help="Directory with per-run FASTQ subdirs (default: data/cellranger)",
    )
    parser.add_argument(
        "-o", "--output-dir",
        type=Path,
        default=Path("output/cellranger"),
        help="Output directory for Cell Ranger results (default: output/cellranger)",
    )
    parser.add_argument(
        "--ref-scrna",
        type=Path,
        required=False,
        help="Path to scRNA mm10 reference (e.g. refdata-gex-mm10-2020-A)",
    )
    parser.add_argument(
        "--ref-atac",
        type=Path,
        required=False,
        help="Path to scATAC mouse reference (e.g. refdata-cellranger-arc-GRCm39-2024-A)",
    )
    parser.add_argument(
        "--runs",
        nargs="+",
        default=None,
        help="Specific run IDs (default: all)",
    )
    parser.add_argument(
        "--no-include-introns",
        action="store_true",
        help="Omit --include-introns for scRNA (use only if pure scRNA, no snRNA)",
    )
    parser.add_argument(
        "--localcores",
        type=int,
        default=None,
        help="Limit cores for Cell Ranger (e.g. 16)",
    )
    parser.add_argument(
        "--localmem",
        type=int,
        default=None,
        help="Limit memory in GB (e.g. 64)",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Print commands without running",
    )
    args = parser.parse_args()

    runs = args.runs or list(CELLRANGER_RUNS)
    for r in runs:
        if r not in CELLRANGER_RUNS:
            print(f"Unknown run: {r}", file=sys.stderr)
            sys.exit(1)

    fastq_dir = args.fastq_dir.resolve()
    output_dir = args.output_dir.resolve()
    output_dir.mkdir(parents=True, exist_ok=True)

    scrna_runs = [r for r in runs if "scATAC" not in CELLRANGER_RUNS[r][0]]
    atac_runs = [r for r in runs if "scATAC" in CELLRANGER_RUNS[r][0]]

    if scrna_runs:
        if not args.ref_scrna:
            print("Error: --ref-scrna is required when running scRNA/snRNA runs", file=sys.stderr)
            sys.exit(1)
        if not args.ref_scrna.exists():
            print(f"Error: scRNA reference not found: {args.ref_scrna}", file=sys.stderr)
            sys.exit(1)
        ref_scrna = args.ref_scrna.resolve()
    else:
        ref_scrna = None

    if atac_runs:
        if not args.ref_atac:
            print("Error: --ref-atac is required when running scATAC runs", file=sys.stderr)
            sys.exit(1)
        if not args.ref_atac.exists():
            print(f"Error: scATAC reference not found: {args.ref_atac}", file=sys.stderr)
            sys.exit(1)
        ref_atac = args.ref_atac.resolve()
    else:
        ref_atac = None

    failed = []

    if scrna_runs:
        print(f"\n--- scRNA/snRNA-seq ({len(scrna_runs)} runs) ---")
        for i, run_id in enumerate(scrna_runs, 1):
            label = CELLRANGER_RUNS[run_id][0]
            run_dir = fastq_dir / run_id
            print(f"[{i}/{len(scrna_runs)}] {run_id} ({label})")
            if not create_symlinks(run_dir, run_id):
                print(f"  Skip: no FASTQs in {run_dir}")
                continue
            ok = run_cellranger_count(
                run_id, run_dir, output_dir, ref_scrna,
                include_introns=not args.no_include_introns,
                localcores=args.localcores, localmem=args.localmem, dry_run=args.dry_run,
            )
            if not ok:
                failed.append(run_id)

    if atac_runs:
        print(f"\n--- scATAC-seq ({len(atac_runs)} runs) ---")
        for i, run_id in enumerate(atac_runs, 1):
            label = CELLRANGER_RUNS[run_id][0]
            run_dir = fastq_dir / run_id
            print(f"[{i}/{len(atac_runs)}] {run_id} ({label})")
            if not create_symlinks_atac(run_dir, run_id):
                print(f"  Skip: no FASTQs in {run_dir} (scATAC needs 3 files from SRA --split-files)")
                continue
            ok = run_cellranger_atac_count(
                run_id, run_dir, output_dir, ref_atac,
                localcores=args.localcores, localmem=args.localmem, dry_run=args.dry_run,
            )
            if not ok:
                failed.append(run_id)

    if failed:
        print(f"\nFailed: {failed}", file=sys.stderr)
        sys.exit(1)
    print("\nDone.")


if __name__ == "__main__":
    main()
