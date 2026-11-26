"""
BAM-to-VCF variant calling web app
-----------------------------------
A small Flask front end over the same GATK4 germline variant-calling
logic used in the Snakemake pipeline (see ../workflow/Snakefile), for
someone who wants to call variants on a single BAM without setting up
the whole pipeline: upload a coordinate-sorted, indexed BAM aligned
against the bundled small reference, and get back a filtered VCF plus
a variant summary table in the browser.

This intentionally reuses the same reference and known-sites files as
the Snakemake pipeline (data/reference/, one level up) — the point is
one shared, tested variant-calling core with two front ends: a batch
pipeline for many samples, and a web app for a single ad-hoc one.

Run with:
    export FLASK_APP=app.app
    flask run --debug
or:
    python app/app.py
"""
import gzip
import shutil
import subprocess
import uuid
from pathlib import Path

from flask import Flask, flash, redirect, render_template, request, send_file, url_for

BASE_DIR = Path(__file__).resolve().parent.parent
REFERENCE = BASE_DIR / "data" / "reference" / "reference.fasta"
UPLOAD_DIR = BASE_DIR / "app" / "uploads"
RESULTS_DIR = BASE_DIR / "app" / "results"
ALLOWED_EXTENSIONS = {".bam"}
MAX_CONTENT_LENGTH = 200 * 1024 * 1024  # 200 MB — this is a demo app, not a production upload service

UPLOAD_DIR.mkdir(parents=True, exist_ok=True)
RESULTS_DIR.mkdir(parents=True, exist_ok=True)

app = Flask(__name__)
app.secret_key = "dev-only-not-for-production"  # replace via env var before any real deployment
app.config["MAX_CONTENT_LENGTH"] = MAX_CONTENT_LENGTH


def run(cmd, **kwargs):
    """Run a subprocess, raising with captured output on failure."""
    result = subprocess.run(cmd, capture_output=True, text=True, **kwargs)
    if result.returncode != 0:
        raise RuntimeError(
            f"Command failed ({' '.join(cmd)}):\n{result.stdout}\n{result.stderr}"
        )
    return result


def call_variants(bam_path: Path, job_id: str) -> Path:
    """
    Index the uploaded BAM if needed, run HaplotypeCaller against the
    bundled reference, and return the path to the resulting VCF.

    No BQSR here deliberately — this app is meant for a BAM someone
    already has (already aligned, possibly already recalibrated
    elsewhere), not a full raw-reads pipeline. See the Snakemake
    pipeline for the full BQSR-included workflow.
    """
    job_dir = RESULTS_DIR / job_id
    job_dir.mkdir(parents=True, exist_ok=True)

    bai_path = bam_path.with_suffix(bam_path.suffix + ".bai")
    if not bai_path.exists():
        run(["samtools", "index", str(bam_path)])

    vcf_path = job_dir / "variants.vcf.gz"
    run([
        "gatk", "HaplotypeCaller",
        "-I", str(bam_path),
        "-R", str(REFERENCE),
        "-O", str(vcf_path),
    ])

    filtered_path = job_dir / "variants.filtered.vcf.gz"
    run([
        "gatk", "VariantFiltration",
        "-R", str(REFERENCE),
        "-V", str(vcf_path),
        "--filter-expression", "QD < 2.0", "--filter-name", "QD_lt_2",
        "--filter-expression", "FS > 60.0", "--filter-name", "FS_gt_60",
        "-O", str(filtered_path),
    ])
    return filtered_path


def summarize_vcf(vcf_path: Path) -> dict:
    """Parse the VCF directly for a quick summary table (see the same
    logic in workflow/scripts/summarize_variants.py — kept here as a
    plain function rather than importing across the app/pipeline
    boundary, so this app has no dependency on the Snakemake side)."""
    total = passed = snps = indels = 0
    rows = []
    with gzip.open(vcf_path, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            chrom, pos, ref, alt, qual, flt = (
                fields[0], fields[1], fields[3], fields[4], fields[5], fields[6],
            )
            total += 1
            is_pass = flt == "PASS"
            passed += is_pass
            is_snp = len(ref) == 1 and all(len(a) == 1 for a in alt.split(","))
            snps += is_snp
            indels += not is_snp
            if is_pass:
                rows.append({"chrom": chrom, "pos": pos, "ref": ref, "alt": alt,
                              "qual": round(float(qual), 1), "type": "SNP" if is_snp else "INDEL"})
    return {
        "total": total, "passed": passed, "filtered_out": total - passed,
        "snps": snps, "indels": indels, "rows": rows[:50],  # cap the table for display
        "truncated": len(rows) > 50,
    }


@app.route("/", methods=["GET"])
def index():
    return render_template("index.html")


@app.route("/call", methods=["POST"])
def call():
    file = request.files.get("bam_file")
    if not file or file.filename == "":
        flash("Please choose a BAM file.")
        return redirect(url_for("index"))
    if Path(file.filename).suffix.lower() not in ALLOWED_EXTENSIONS:
        flash("Only .bam files are accepted. Sort and index it first (samtools sort / index).")
        return redirect(url_for("index"))

    job_id = uuid.uuid4().hex[:12]
    bam_path = UPLOAD_DIR / f"{job_id}.bam"
    file.save(bam_path)

    try:
        vcf_path = call_variants(bam_path, job_id)
        summary = summarize_vcf(vcf_path)
    except RuntimeError as exc:
        flash(f"Variant calling failed: {exc}")
        return redirect(url_for("index"))
    finally:
        bam_path.unlink(missing_ok=True)
        bai_path = bam_path.with_suffix(bam_path.suffix + ".bai")
        bai_path.unlink(missing_ok=True)

    return render_template("results.html", job_id=job_id, summary=summary)


@app.route("/download/<job_id>")
def download(job_id):
    vcf_path = RESULTS_DIR / job_id / "variants.filtered.vcf.gz"
    if not vcf_path.exists():
        flash("That result has expired or never existed.")
        return redirect(url_for("index"))
    return send_file(vcf_path, as_attachment=True, download_name=f"{job_id}.filtered.vcf.gz",
                      mimetype="application/gzip")


if __name__ == "__main__":
    app.run(debug=True)
