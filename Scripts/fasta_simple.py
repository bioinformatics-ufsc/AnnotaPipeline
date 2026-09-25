#!/usr/bin/env python3
"""
fasta_simple.py - Simple FASTA header annotator for AnnotaPipeline.

Formats FASTA sequence headers with organism and annotation information:
>{sequence_id} | Organism: {organism} | Description: {description}
"""

import argparse
from pathlib import Path
import sys
import warnings
from Bio import SeqIO


def cli() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Join sequences and annotations into FASTA headers.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "-annot", "--annotations", dest="annot",
        metavar="[annotations.txt]",
        required=True,
        help="File containing ID and annotation separated by tab (\\t).",
    )
    parser.add_argument(
        "-b", "--basename", dest="basename",
        metavar="[basename]",
        required=True,
        help="Output file basename prefix.",
    )
    parser.add_argument(
        "-faf", "--fasta", dest="seqs",
        metavar="[fasta.fasta]",
        required=True,
        help="FASTA file with protein sequences.",
    )
    parser.add_argument(
        "-org", "--organism", dest="org",
        metavar="[organism]",
        required=True,
        help="Organism name (e.g. 'Homo sapiens').",
    )
    return parser


def parse_annotations(annot_file: Path) -> dict[str, str]:
    """Parse tab-separated ID and annotation file into a dictionary."""
    annotations = {}
    with open(annot_file, "r", encoding="utf-8", errors="replace") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t", 1)
            if len(parts) == 2:
                annotations[parts[0].strip()] = parts[1].strip()
            elif len(parts) == 1:
                annotations[parts[0].strip()] = ""
    return annotations


def write_simple_fasta(
    annot_file: str | Path,
    fasta_file: str | Path,
    basename: str,
    organism: str,
    output_dir: str | Path | None = None,
) -> Path:
    """Format fasta headers with organism and description."""
    annot_path = Path(annot_file)
    fasta_path = Path(fasta_file)
    out_dir = Path(output_dir) if output_dir else Path.cwd()
    out_dir.mkdir(parents=True, exist_ok=True)

    org_clean = str(organism).strip("\"'")
    annotations = parse_annotations(annot_path)

    out_fasta = out_dir / f"AnnotaPipeline_{basename}_proteins.fasta"
    ids_without_annot = []

    with open(fasta_path, "r", encoding="utf-8", errors="replace") as fasta_handle, \
         open(out_fasta, "w", encoding="utf-8") as out_handle:
        for record in SeqIO.parse(fasta_handle, "fasta"):
            seq_id = record.id.strip()
            desc = annotations.get(seq_id)
            if desc is not None:
                new_header = f"{seq_id} | Organism: {org_clean} | Description: {desc}"
                out_handle.write(f">{new_header}\n{str(record.seq).upper()}\n")
            else:
                ids_without_annot.append(seq_id)

    if ids_without_annot:
        warn_file = out_dir / f"{basename}_ids_with_no_annotations.txt"
        with open(warn_file, "w", encoding="utf-8") as wf:
            for item in ids_without_annot:
                wf.write(f"{item}\n")
        warnings.warn(
            f"WARNING: {len(ids_without_annot)} sequence(s) lacked annotations. IDs saved to {warn_file}",
            stacklevel=2,
        )

    return out_fasta


def main(argv: list[str] | None = None) -> int:
    parser = cli()
    args = parser.parse_args(argv)
    write_simple_fasta(
        annot_file=args.annot,
        fasta_file=args.seqs,
        basename=args.basename,
        organism=args.org,
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
