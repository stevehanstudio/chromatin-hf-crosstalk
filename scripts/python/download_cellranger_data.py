#!/usr/bin/env python3
"""
Download scRNA-seq and scATAC-seq data for Cell Ranger processing.

Run this on an x86 machine where Cell Ranger is available.
For bulk RNA-seq, ChIP-seq, CUT&RUN → use scripts/python/download_data.py instead.
Data: Alexanian et al. Nature 2024 | GEO GSE221699 | BioProject PRJNA915384

Includes:
  - scRNA-seq (10x): Sham/TAC/JQ1, CD45+, TAC WT/Brd4KO, TAC IgG/anti-IL1B
  - snRNA-seq (10x): Whole heart TAC WT/Brd4KO
  - scATAC-seq (10x): Whole heart + CD45+ nuclei

Run from project root:
  python scripts/python/download_cellranger_data.py --minimal   # 15 runs; see disk note below
  python scripts/python/download_cellranger_data.py             # all 31 runs

Disk (uncompressed .fastq): scATAC runs are very large (often hundreds of GB each).
``--minimal`` still includes 5 scATAC libraries — plan on the order of 1–4 TB for FASTQs alone,
plus temp space for fasterq-dump. Older "~150 GB total" estimates were wrong for this layout.

If fasterq-dump exits with code 3 / ``rcNotFound`` during concat (large scATAC runs):
  free disk space (temp + output need headroom), remove ``fasterq.tmp.*`` dirs in CWD,
  pass ``--fasterq-temp-dir /path/to/big/scratch``, and retry failed SRR with ``--runs``.
"""

import argparse
import hashlib
import os
import subprocess
import sys
import tempfile
import time
from pathlib import Path

# ENA run metadata for Cell Ranger data (scRNA-seq, snRNA-seq, scATAC-seq)
# Excludes: Bulk RNA-seq, ChIP-seq, CUT&RUN, Sorted FBs (bulk)
# Source: https://www.ebi.ac.uk/ena/portal/api/filereport?accession=PRJNA915384&result=read_run
# scATAC run IDs - require SRA with --split-files (ENA only has 2 FASTQs; 10x needs R1,R2,I2)
SCATAC_RUNS = {
    "SRR22882159", "SRR22882160", "SRR22882161", "SRR22882162", "SRR22882163",
    "SRR22882164", "SRR22882165", "SRR22882166", "SRR22882167",
}

# Minimal subset (15 runs): one rep per condition for Figs 1-5 — includes 5 scATAC + 10 scRNA/snRNA
# (disk: order ~1–4 TB uncompressed FASTQs; scATAC dominates, not ~150 GB total)
# scATAC: whole heart TAC/Brd4KO, CD45+ Sham/TAC/Brd4KO
# scRNA: Sham/TAC/JQ1, CD45+, snRNA TAC/Brd4KO, IgG/anti-IL1B
MINIMAL_RUNS = [
    "SRR22882159", "SRR22882161",  # scATAC whole heart TAC_WT, TAC_BRD4KO
    "SRR22882163", "SRR22882164", "SRR22882166",  # scATAC CD45+ Sham, TAC_WT, TAC_BRD4KO
    "SRR22882168", "SRR22882171", "SRR22882174",  # scRNA Sham, TAC, TAC_JQ1
    "SRR22882177", "SRR22882178", "SRR22882180",  # scRNA CD45+ Sham, TAC, Brd4KO
    "SRR22882182", "SRR22882184",  # snRNA TAC_WT, TAC_BRD4KO
    "SRR22882186", "SRR22882188",  # scRNA TAC IgG, anti-IL1B
]

CELLRANGER_RUNS = {
    # scATAC-seq (10 runs)
    "SRR22882159": ("scATAC-Seq TAC_WT_rep1", "paired"),
    "SRR22882160": ("scATAC-Seq TAC_WT_rep2", "paired"),
    "SRR22882161": ("scATAC-Seq TAC_BRD4KO_rep1", "paired"),
    "SRR22882162": ("scATAC-Seq TAC_BRD4KO_rep2", "paired"),
    "SRR22882163": ("scATAC-Seq CD45pos_Sham", "paired"),
    "SRR22882164": ("scATAC-Seq CD45pos_TAC_WT_rep1", "paired"),
    "SRR22882165": ("scATAC-Seq CD45pos_TAC_WT_rep2", "paired"),
    "SRR22882166": ("scATAC-Seq CD45pos_TAC_BRD4KO_rep1", "paired"),
    "SRR22882167": ("scATAC-Seq CD45pos_TAC_BRD4KO_rep2", "paired"),
    # scRNA-seq / snRNA-seq (22 runs)
    "SRR22882168": ("scRNA-Seq Sham_rep1", "paired"),
    "SRR22882169": ("scRNA-Seq Sham_rep2", "paired"),
    "SRR22882170": ("scRNA-Seq Sham_rep3", "paired"),
    "SRR22882171": ("scRNA-Seq TAC_rep1", "paired"),
    "SRR22882172": ("scRNA-Seq TAC_rep2", "paired"),
    "SRR22882173": ("scRNA-Seq TAC_rep3", "paired"),
    "SRR22882174": ("scRNA-Seq TAC_JQ1_rep1", "paired"),
    "SRR22882175": ("scRNA-Seq TAC_JQ1_rep2", "paired"),
    "SRR22882176": ("scRNA-Seq TAC_JQ1_rep3", "paired"),
    "SRR22882177": ("scRNA-Seq CD45+ Sham", "paired"),
    "SRR22882178": ("scRNA-Seq CD45pos_TAC_WT_rep1", "paired"),
    "SRR22882179": ("scRNA-Seq CD45pos_TAC_WT_rep2", "paired"),
    "SRR22882180": ("scRNA-Seq CD45pos_TAC_BRD4KO_rep1", "paired"),
    "SRR22882181": ("scRNA-Seq CD45pos_TAC_BRD4KO_rep2", "paired"),
    "SRR22882182": ("scRNA-Seq TAC_WT_rep1", "paired"),
    "SRR22882183": ("scRNA-Seq TAC_WT_2", "paired"),
    "SRR22882184": ("scRNA-Seq TAC_BRD4KO_rep1", "paired"),
    "SRR22882185": ("scRNA-Seq TAC_BRD4KO_rep2", "paired"),
    "SRR22882186": ("scRNA-Seq TAC_IgG_Ab1", "paired"),
    "SRR22882187": ("scRNA-Seq TAC_IgG_Ab2", "paired"),
    "SRR22882188": ("scRNA-Seq TAC_antiIL1B_Ab1", "paired"),
    "SRR22882189": ("scRNA-Seq TAC_antiIL1B_Ab2", "paired"),
}

def download_with_curl(
    url: str,
    out_path: Path,
    expected_md5: str | None = None,
    *,
    max_time: int = 21600,
    retries: int = 3,
) -> bool:
    """Download file using curl. Uses 6h timeout, resume (-C -), and retries for large files."""
    out_path.parent.mkdir(parents=True, exist_ok=True)
    cmd = [
        "curl", "-f", "-L", "-C", "-",
        # ENA links for very large FASTQs can drop mid-transfer; these flags improve resiliency.
        "--retry-all-errors",
        "--retry-connrefused",
        "--http1.1",
        "--max-time", str(max_time),
        "--connect-timeout", "120",
        "--retry", "2", "--retry-delay", "30",
        "-o", str(out_path), url,
    ]
    for attempt in range(1, retries + 1):
        try:
            result = subprocess.run(cmd, capture_output=True)
            if result.returncode == 0:
                break
            raise subprocess.CalledProcessError(result.returncode, cmd, result.stdout, result.stderr)
        except subprocess.CalledProcessError as e:
            stderr = (e.stderr or b"").decode(errors="replace")
            print(f"  Error (attempt {attempt}/{retries}): {stderr.strip() or str(e)}", file=sys.stderr)
            if e.returncode in (18, 28) and out_path.exists() and out_path.stat().st_size > 1000:
                print(f"  Keeping partial for resume", file=sys.stderr)
            elif out_path.exists():
                out_path.unlink()
            if attempt < retries:
                wait = 60 * attempt
                print(f"  Retrying in {wait}s...", file=sys.stderr)
                time.sleep(wait)
            else:
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
        # Prefer HTTPS for ENA paths; FTP can fail on some networks with
        # "Server denied you to change to the given directory" (curl 9).
        if u.startswith("http://") or u.startswith("https://"):
            urls.append(u)
        elif u.startswith("ftp://"):
            urls.append("https://" + u[len("ftp://"):])
        else:
            urls.append(f"https://{u}")
    md5_list = (fastq_md5.split(";") if fastq_md5 else [])
    md5s = [md5_list[i] if i < len(md5_list) else None for i in range(len(urls))]

    run_dir = out_dir / run_id
    run_dir.mkdir(parents=True, exist_ok=True)
    ok = True
    MIN_VALID_SIZE = 1000  # bytes; smaller = likely partial/failed
    for i, (u, m) in enumerate(zip(urls, md5s)):
        fname = Path(u).name
        out_path = run_dir / fname
        if out_path.exists():
            if out_path.stat().st_size < MIN_VALID_SIZE:
                out_path.unlink()  # Remove tiny/empty partial so we re-download
            elif not m or hashlib.md5(out_path.read_bytes()).hexdigest() == m:
                print(f"  Skip (exists): {fname}")
                continue
            else:
                out_path.unlink()  # MD5 mismatch; remove so we re-download
        print(f"  Downloading {fname} ...")
        if not download_with_curl(u, out_path, m if md5_checks else None):
            ok = False
    return ok


def _clear_partial_fasterq_outputs(run_dir: Path, run_id: str) -> None:
    """Remove incomplete fasterq-dump .fastq so a retry is not skipped as 'done'."""
    for p in run_dir.glob(f"{run_id}_*.fastq"):
        try:
            p.unlink()
        except OSError:
            pass


def download_with_fasterq_dump(
    run_id: str,
    out_dir: Path,
    scatac: bool = False,
    *,
    temp_dir: Path | None = None,
    threads: int | None = None,
    retries: int = 2,
) -> bool:
    """Use SRA fasterq-dump (requires sra-tools).
    For scATAC: use --split-files --include-technical to get R1,R2,I2 (ENA has only 2 FASTQs).
    Skips if expected output files already exist (resume-friendly).

    ``temp_dir`` should be on a volume with ample free space (scATAC can use 100+ GB
    transient files during concat). If None, uses CHROMATIN_HF_FASTERQ_TEMP, TMPDIR,
    or system temp.

    Exit code 3 with ``KDirectoryFileSize`` / ``rcNotFound`` usually means the temp
    volume filled up or temp chunks were lost — use a larger ``--fasterq-temp-dir``.
    """
    run_dir = out_dir / run_id
    run_dir.mkdir(parents=True, exist_ok=True)
    min_size = 10_000  # bytes; avoid skipping tiny partials
    if scatac:
        expected = [run_dir / f"{run_id}_{i}.fastq" for i in (1, 2, 3)]
    else:
        expected = [run_dir / f"{run_id}_{i}.fastq" for i in (1, 2)]
    if all(p.exists() and p.stat().st_size >= min_size for p in expected):
        print(f"  Skip (exists): {', '.join(p.name for p in expected)}")
        return True

    if temp_dir is None:
        env_tmp = os.environ.get("CHROMATIN_HF_FASTERQ_TEMP") or os.environ.get("TMPDIR")
        if env_tmp:
            base = Path(env_tmp)
        else:
            base = Path(tempfile.gettempdir())
        temp_dir = base / "chromatin_hf_fasterq"
    else:
        temp_dir = Path(temp_dir)
    temp_dir.mkdir(parents=True, exist_ok=True)

    nthr = threads if threads is not None else min(16, max(1, os.cpu_count() or 4))

    cmd = [
        "fasterq-dump",
        run_id,
        "-O", str(run_dir),
        "-f",
        "-e", str(nthr),
        "--temp", str(temp_dir),
    ]
    if scatac:
        cmd.extend(["--split-files", "--include-technical"])

    env = os.environ.copy()
    env["TMPDIR"] = str(temp_dir)

    last_err: Exception | None = None
    for attempt in range(1, retries + 1):
        try:
            subprocess.run(cmd, check=True, env=env)
            if all(p.exists() and p.stat().st_size >= min_size for p in expected):
                return True
            print(
                f"  fasterq-dump finished but outputs missing or tiny (attempt {attempt}/{retries})",
                file=sys.stderr,
            )
        except (subprocess.CalledProcessError, FileNotFoundError) as e:
            last_err = e
            print(f"  fasterq-dump failed (attempt {attempt}/{retries}): {e}", file=sys.stderr)
        if attempt < retries:
            _clear_partial_fasterq_outputs(run_dir, run_id)
            wait = 30 * attempt
            print(f"  Retrying in {wait}s after clearing partial FASTQs...", file=sys.stderr)
            time.sleep(wait)

    if last_err:
        print(f"  fasterq-dump gave up after {retries} attempts.", file=sys.stderr)
    return False


def estimate_disk_hint(runs: list[str]) -> str:
    """Order-of-magnitude for uncompressed FASTQs (fasterq-dump output is not gzip)."""
    n_atac = sum(1 for r in runs if r in SCATAC_RUNS)
    n_other = len(runs) - n_atac
    low_gb = n_atac * 200 + n_other * 8
    high_gb = n_atac * 700 + n_other * 40
    tb_lo, tb_hi = low_gb / 1024, high_gb / 1024
    return (
        f"~{tb_lo:.1f}–{tb_hi:.1f} TB FASTQs ({n_atac} scATAC + {n_other} scRNA/snRNA; wide range per run). "
        "Budget extra temp space for fasterq-dump during scATAC merges."
    )


def main():
    parser = argparse.ArgumentParser(
        description="Download scRNA-seq and scATAC-seq data for Cell Ranger (x86)"
    )
    parser.add_argument(
        "-o", "--output-dir",
        type=Path,
        default=Path("data/cellranger"),
        help="Output directory (default: data/cellranger)",
    )
    parser.add_argument(
        "--use-sra",
        action="store_true",
        help="Use fasterq-dump instead of ENA (required for scATAC: ENA has only 2 FASTQs)",
    )
    parser.add_argument(
        "--minimal",
        action="store_true",
        help="Download minimal subset (15 runs, ~1–4 TB FASTQs — mostly scATAC) for Figs 1-5",
    )
    parser.add_argument(
        "--runs",
        nargs="+",
        default=None,
        help="Specific run IDs (overrides --minimal). Default: all 31 runs (multi-TB FASTQs)",
    )
    parser.add_argument(
        "--md5",
        action="store_true",
        help="Verify MD5 checksums (slower)",
    )
    parser.add_argument(
        "--fasterq-temp-dir",
        type=Path,
        default=None,
        help="Directory for fasterq-dump temp files (large scATAC runs need 100+ GB free; "
        "default: CHROMATIN_HF_FASTERQ_TEMP, TMPDIR, or system temp)",
    )
    parser.add_argument(
        "--fasterq-threads",
        type=int,
        default=None,
        help="Threads for fasterq-dump -e (default: min(16, CPUs); lower if OOM)",
    )
    parser.add_argument(
        "--fasterq-retries",
        type=int,
        default=2,
        help="Retries per SRR if fasterq-dump fails (default: 2)",
    )
    args = parser.parse_args()

    if args.runs:
        runs = args.runs
    elif args.minimal:
        runs = MINIMAL_RUNS
        print("Using --minimal: 15 runs for Figs 1-5 (includes 5 large scATAC libraries)\n")
    else:
        runs = list(CELLRANGER_RUNS)
    for r in runs:
        if r not in CELLRANGER_RUNS:
            print(f"Unknown run: {r}", file=sys.stderr)
            sys.exit(1)

    args.output_dir.mkdir(parents=True, exist_ok=True)
    print(f"Downloading {len(runs)} runs to {args.output_dir}")
    print(f"Estimated FASTQ disk: {estimate_disk_hint(runs)}\n")
    if args.fasterq_temp_dir:
        print(f"fasterq-dump temp: {args.fasterq_temp_dir.resolve()}\n")

    failed = []
    for i, run_id in enumerate(runs, 1):
        label = CELLRANGER_RUNS[run_id][0]
        is_scatac = run_id in SCATAC_RUNS
        print(f"[{i}/{len(runs)}] {run_id} ({label})")
        if args.use_sra or is_scatac:
            if is_scatac and not args.use_sra:
                print("  (scATAC requires SRA; ENA has only 2 FASTQs, cellranger-atac needs I2)")
            ok = download_with_fasterq_dump(
                run_id,
                args.output_dir,
                scatac=is_scatac,
                temp_dir=args.fasterq_temp_dir,
                threads=args.fasterq_threads,
                retries=max(1, args.fasterq_retries),
            )
        else:
            ok = download_run_ena(run_id, args.output_dir, md5_checks=args.md5)
        if not ok:
            failed.append(run_id)

    if failed:
        print(f"\nFailed: {failed}", file=sys.stderr)
        sys.exit(1)
    print("\nDone.")


if __name__ == "__main__":
    main()
