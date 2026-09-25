#!/usr/bin/env python3
"""
kallisto_parser.py - Filter Kallisto abundance output by TPM threshold.
"""

import argparse
from pathlib import Path
import sys
import pandas as pd


def cli() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Filter Kallisto abundance results by TPM threshold (mean, median, or specific value).",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    required_args = parser.add_argument_group("required arguments")
    required_args.add_argument(
        "-ktfile", "--kallisto-file", dest="ktfile",
        metavar="[abundance.tsv]",
        required=True,
        help="Kallisto abundance TSV output file.",
    )
    required_args.add_argument(
        "-basename", dest="basename",
        metavar="[basename]",
        required=True,
        help="Basename prefix for output file.",
    )

    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument(
        "-tpmavg", "--tpm-average", dest="tpm_avg", action="store_true",
        help="Filter using average TPM as threshold.",
    )
    group.add_argument(
        "-tpmmd", "--tpm-median", dest="tpm_median", action="store_true",
        help="Filter using median TPM as threshold.",
    )
    group.add_argument(
        "-tpmval", "--tpm-threshold", dest="tpm_value", type=float,
        metavar="VALUE",
        help="Filter using a specific TPM numeric threshold.",
    )
    return parser


def parse_kallisto(
    kallisto_file: str | Path,
    basename: str,
    method: str = "median",
    threshold_value: float | None = None,
    output_dir: str | Path | None = None,
) -> Path:
    """Filter Kallisto abundance dataframe by TPM and save to TSV."""
    kt_path = Path(kallisto_file)
    out_dir = Path(output_dir) if output_dir else Path.cwd()
    out_dir.mkdir(parents=True, exist_ok=True)

    if not kt_path.exists():
        raise FileNotFoundError(f"Kallisto output file not found: {kt_path}")

    df = pd.read_csv(kt_path, sep="\t")

    if method == "mean":
        cutoff = float(df["tpm"].mean())
    elif method == "median":
        cutoff = float(df["tpm"].median())
    elif method == "value":
        if threshold_value is None:
            raise ValueError("threshold_value must be provided when method is 'value'")
        cutoff = float(threshold_value)
    else:
        raise ValueError(f"Unknown threshold method: {method}")

    filtered = df[df["tpm"] >= cutoff].copy()
    columns_to_drop = [c for c in ["length", "eff_length", "est_counts"] if c in filtered.columns]
    filtered.drop(columns=columns_to_drop, inplace=True)

    out_file = out_dir / f"{basename}_Transcript_Quantification.tsv"
    filtered.to_csv(out_file, sep="\t", index=False)
    return out_file


def main(argv: list[str] | None = None) -> int:
    parser = cli()
    args = parser.parse_args(argv)

    if args.tpm_avg:
        method = "mean"
        val = None
    elif args.tpm_median:
        method = "median"
        val = None
    else:
        method = "value"
        val = args.tpm_value

    parse_kallisto(
        kallisto_file=args.ktfile,
        basename=args.basename,
        method=method,
        threshold_value=val,
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())