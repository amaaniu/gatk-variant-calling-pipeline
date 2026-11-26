#!/usr/bin/env python3
"""
Build a one-row-per-sample summary of the variant-calling run:
mapping rate (from flagstat, pre-BQSR is identical post-BQSR since BQSR
doesn't change alignments) and variant counts by type (from bcftools
stats, computed directly here to keep this script's only dependency
on file parsing, not on a specific bcftools stats section format).
"""
import argparse
import gzip
import re
from pathlib import Path

import pandas as pd


def parse_flagstat(path):
    text = Path(path).read_text()
    total = re.search(r"^(\d+) \+ \d+ in total", text, re.MULTILINE)
    mapped = re.search(r"^(\d+) \+ \d+ mapped \(", text, re.MULTILINE)
    total = int(total.group(1)) if total else None
    mapped = int(mapped.group(1)) if mapped else None
    rate = round(mapped / total * 100, 2) if total else None
    return total, mapped, rate


def count_variants(vcf_path):
    """Count PASS SNPs and indels directly from the VCF, ref vs alt length."""
    n_snp_pass = n_indel_pass = n_snp_filtered = n_indel_filtered = 0
    with gzip.open(vcf_path, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            ref, alt, flt = fields[3], fields[4], fields[6]
            is_snp = len(ref) == 1 and all(len(a) == 1 for a in alt.split(","))
            passed = flt == "PASS"
            if is_snp:
                n_snp_pass += passed
                n_snp_filtered += not passed
            else:
                n_indel_pass += passed
                n_indel_filtered += not passed
    return n_snp_pass, n_snp_filtered, n_indel_pass, n_indel_filtered


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--vcfs", nargs="+", required=True)
    parser.add_argument("--flagstats", nargs="+", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    flagstat_by_sample = {
        Path(f).name.replace(".flagstat.txt", ""): f for f in args.flagstats
    }

    rows = []
    for vcf in args.vcfs:
        sample = Path(vcf).name.replace(".filtered.vcf.gz", "")
        total, mapped, rate = parse_flagstat(flagstat_by_sample[sample])
        snp_p, snp_f, indel_p, indel_f = count_variants(vcf)
        rows.append({
            "sample": sample,
            "total_reads": total,
            "mapped_reads": mapped,
            "mapping_rate_pct": rate,
            "snps_pass": snp_p,
            "snps_filtered_out": snp_f,
            "indels_pass": indel_p,
            "indels_filtered_out": indel_f,
        })

    df = pd.DataFrame(rows)
    df.to_csv(args.output, sep="\t", index=False)
    print(df.to_string(index=False))


if __name__ == "__main__":
    main()
