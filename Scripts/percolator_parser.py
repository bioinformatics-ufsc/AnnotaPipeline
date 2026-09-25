#!/usr/bin/env python3
"""
percolator_parser.py - Filter and parse Percolator peptide identification outputs.

Filters PSMs by q-value threshold, discards decoy hits, and generates
a simplified ProteinID -> Peptide -> Spectrum mapping file.
"""

import argparse
from pathlib import Path
import sys
import warnings


def cli() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Filter Percolator output by q-value and extract protein-peptide-spectrum mappings.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    required = parser.add_argument_group("required arguments")
    required.add_argument(
        "-p", "--percolator-output", dest="proteom",
        metavar="[percolator.tsv]",
        required=True,
        help="Percolator output file.",
    )
    required.add_argument(
        "-b", "--basename", dest="basename",
        metavar="[basename]",
        required=True,
        help="Basename prefix for output files.",
    )
    required.add_argument(
        "-qv", "-q-value", dest="qvalue",
        type=float,
        default=0.05,
        help="Maximum q-value threshold to accept a peptide (0.0 to 1.0).",
    )
    return parser


def parse_percolator(
    percolator_file: str | Path,
    basename: str,
    qvalue: float = 0.05,
    output_dir: str | Path | None = None,
) -> tuple[Path, Path]:
    """Filter Percolator output by q-value cutoff and produce parsed TSV."""
    p_path = Path(percolator_file)
    out_dir = Path(output_dir) if output_dir else Path.cwd()
    out_dir.mkdir(parents=True, exist_ok=True)

    if not (0.0 <= qvalue <= 1.0):
        raise ValueError(f"q-value must be between 0.0 and 1.0, got: {qvalue}")

    filtered_file = out_dir / f"{basename}_qvalue_filtered.tsv"
    parsed_file = out_dir / f"{basename}_parsed.tsv"

    # Step 1: Filter lines by q-value (column index 2)
    filtered_rows = []
    with open(p_path, "r", encoding="utf-8", errors="replace") as f_in:
        header = f_in.readline()
        for line in f_in:
            parts = line.rstrip("\r\n").split("\t")
            if len(parts) > 2:
                try:
                    qv = float(parts[2])
                    if qv <= qvalue:
                        filtered_rows.append(parts)
                except ValueError:
                    continue

    with open(filtered_file, "w", encoding="utf-8") as f_qval:
        f_qval.write(header)
        for parts in filtered_rows:
            f_qval.write("\t".join(parts) + "\n")

    if not filtered_rows:
        warnings.warn(f"No valid PSMs found below q-value threshold {qvalue} in {p_path}", stacklevel=2)
        # Create empty parsed file with header
        with open(parsed_file, "w", encoding="utf-8") as f_out:
            f_out.write("ProteinID\tPeptide\tSpectrum\n")
        return filtered_file, parsed_file

    # Step 2: Extract ProteinID, Peptide, Spectrum
    # Format of PSM id in Percolator: path/basename_spectrum_charge_num
    protein_map: dict[str, list[tuple[str, str]]] = {}

    for row in filtered_rows:
        if len(row) < 5:
            continue
        psm_id = row[0].split("/")[-1]
        tokens = psm_id.rsplit("_", 3)
        spectrum = tokens[1] if len(tokens) >= 4 else tokens[0]

        peptide = row[4]
        protein_ids = row[5:] if len(row) > 5 else []

        for pid in protein_ids:
            clean_pid = pid.strip()
            if "DECOY" in clean_pid:
                continue
            protein_map.setdefault(clean_pid, []).append((peptide, spectrum))

    # Step 3: Write parsed output
    with open(parsed_file, "w", encoding="utf-8") as f_out:
        f_out.write("ProteinID\tPeptide\tSpectrum\n")
        for pid in sorted(protein_map.keys()):
            for pep, spec in protein_map[pid]:
                f_out.write(f"{pid}\t{pep}\t{spec}\n")

    return filtered_file, parsed_file


def main(argv: list[str] | None = None) -> int:
    parser = cli()
    args = parser.parse_args(argv)
    try:
        parse_percolator(
            percolator_file=args.proteom,
            basename=args.basename,
            qvalue=args.qvalue,
        )
    except Exception as exc:
        print(f"Error parsing Percolator file: {exc}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())