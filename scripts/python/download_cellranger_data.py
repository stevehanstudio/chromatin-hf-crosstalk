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

Very large ENA HTTPS downloads sometimes fail with ``curl: (56) OpenSSL ... unexpected eof``.
The downloader keeps partial files for resume (``curl -C -`` / ``wget -c``). If curl keeps failing,
try ``--ena-backend wget`` or ``--use-sra`` (fasterq-dump; uncompressed ``.fastq`` — run_cellranger
symlink helper supports that).
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
    expected_bytes: int | None = None,
    *,
    max_time: int = 0,
    retries: int = 20,
) -> bool:
    """
    Download file using curl. Uses resume (-C -) and many short retries (large ENA drops are
    mid-transfer TLS EOFs — exit 56, etc.). Default ``max_time=0`` means no per-invocation cap
    so multi-hour resumes are possible.

    If `expected_bytes` is provided, a size mismatch is treated as an incomplete/truncated
    download (even if `--md5` is not enabled).
    """
    out_path.parent.mkdir(parents=True, exist_ok=True)
    curl_max_time = str(max_time) if max_time > 0 else "0"
    cmd = [
        "curl", "-f", "-L", "-C", "-",
        # ENA links for very large FASTQs can drop mid-transfer; these flags improve resiliency.
        "--retry-all-errors",
        "--retry-connrefused",
        "--http1.1",
        "--max-time", curl_max_time,
        "--connect-timeout", "120",
        "--retry", "5",
        "--retry-delay", "10",
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
            # Never delete a substantive partial: curl -C - needs it to resume (old bug: only
            # exit 18/28 kept partials, so TLS EOF=56 restarted from byte 0 every time).
            if out_path.exists():
                sz = out_path.stat().st_size
                if sz <= 1000:
                    out_path.unlink()
                else:
                    print(f"  Keeping partial ({sz / 1024**3:.2f} GB) for resume", file=sys.stderr)
            if attempt < retries:
                wait = min(300, 30 * attempt)
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
    if expected_bytes is not None:
        got = out_path.stat().st_size
        if got != expected_bytes:
            print(f"  Size mismatch: got {got} bytes, expected {expected_bytes} bytes", file=sys.stderr)
            return False
    return True


def download_with_wget(
    url: str,
    out_path: Path,
    expected_md5: str | None = None,
    expected_bytes: int | None = None,
    *,
    retries: int = 20,
) -> bool:
    """
    Download with GNU wget ``-c`` (continue). Useful when curl hits repeated TLS EOF on ENA.
    """
    out_path.parent.mkdir(parents=True, exist_ok=True)
    cmd = [
        "wget",
        "-O", str(out_path),
        "-c",
        "--timeout=300",
        "--read-timeout=300",
        "--tries=1",
        url,
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
            if out_path.exists():
                sz = out_path.stat().st_size
                if sz <= 1000:
                    out_path.unlink()
                else:
                    print(f"  Keeping partial ({sz / 1024**3:.2f} GB) for resume", file=sys.stderr)
            if attempt < retries:
                wait = min(300, 30 * attempt)
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
    if expected_bytes is not None:
        got = out_path.stat().st_size
        if got != expected_bytes:
            print(f"  Size mismatch: got {got} bytes, expected {expected_bytes} bytes", file=sys.stderr)
            return False
    return True


def download_run_ena(
    run_id: str,
    out_dir: Path,
    md5_checks: bool = False,
    *,
    backend: str = "curl",
) -> bool:
    """Fetch FASTQ URLs from ENA and download."""
    import urllib.request

    url = (
        "https://www.ebi.ac.uk/ena/portal/api/filereport"
        f"?accession={run_id}&result=read_run&fields=fastq_ftp,fastq_md5,fastq_bytes"
    )
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
    fastq_bytes = data.get("fastq_bytes", "")

    if not fastq_ftp:
        print(f"  No FASTQ FTP for {run_id}", file=sys.stderr)
        return False

    urls = []
    for u in fastq_ftp.split(";"):
        u = u.strip()
        if not u:
            continue
        # Backend-aware URL normalization:
        # - curl backend: prefer HTTPS (avoids some ftp directory errors).
        # - wget backend: prefer native FTP because some ENA HTTPS endpoints return 403 to wget.
        if u.startswith("http://") or u.startswith("https://"):
            urls.append(u)
        elif u.startswith("ftp://"):
            if backend == "wget":
                urls.append(u)
            else:
                urls.append("https://" + u[len("ftp://"):])
        else:
            if backend == "wget":
                urls.append(f"ftp://{u}")
            else:
                urls.append(f"https://{u}")
    md5_list = (fastq_md5.split(";") if fastq_md5 else [])
    md5s = [md5_list[i] if i < len(md5_list) else None for i in range(len(urls))]
    bytes_list = (fastq_bytes.split(";") if fastq_bytes else [])
    sizes = [
        int(bytes_list[i]) if i < len(bytes_list) and bytes_list[i] else None
        for i in range(len(urls))
    ]

    run_dir = out_dir / run_id
    run_dir.mkdir(parents=True, exist_ok=True)
    ok = True
    MIN_VALID_SIZE = 1000  # bytes; smaller = likely partial/failed
    for i, (u, m) in enumerate(zip(urls, md5s)):
        fname = Path(u).name
        expected_size = sizes[i] if i < len(sizes) else None
        out_path = run_dir / fname
        if expected_size is not None:
            print(f"  Expected {fname}: {expected_size/1024**3:.2f} GB")
        if out_path.exists():
            if out_path.stat().st_size < MIN_VALID_SIZE:
                out_path.unlink()  # Remove tiny/empty partial so we re-download
            elif expected_size is not None and out_path.stat().st_size == expected_size:
                # When MD5 is not enabled, size is a strong signal for a complete download.
                print(f"  Skip (exists, size ok): {fname}")
                continue
            else:
                if md5_checks and m and hashlib.md5(out_path.read_bytes()).hexdigest() == m:
                    print(f"  Skip (exists, md5 ok): {fname}")
                    continue
                out_path.unlink()  # size or md5 mismatch; remove so we re-download
        print(f"  Downloading {fname} ...")
        dl = download_with_wget if backend == "wget" else download_with_curl
        if not dl(
            u,
            out_path,
            m if md5_checks else None,
            expected_bytes=expected_size,
        ):
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
        "--ena-backend",
        choices=("curl", "wget"),
        default="curl",
        help="HTTP client for ENA FASTQs (default: curl). Use wget if curl hits repeated TLS EOF.",
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
            ok = download_run_ena(
                run_id,
                args.output_dir,
                md5_checks=args.md5,
                backend=args.ena_backend,
            )
        if not ok:
            failed.append(run_id)

    if failed:
        print(f"\nFailed: {failed}", file=sys.stderr)
        sys.exit(1)
    print("\nDone.")


if __name__ == "__main__":
    main()
