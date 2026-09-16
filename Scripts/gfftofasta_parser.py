#!/usr/bin/env python3
"""
gfftofasta_parser.py - Add GFF genomic coordinates and annotations to FASTA headers.

Headers follow the format:
>{sequence_id} | Organism: {organism} | Location: {scaffold} | Start: {start} | End: {end} | Strand: {strand} | Description: {description}
"""

import argparse
from pathlib import Path
import re
import sys
import warnings
from Bio import SeqIO


def cli() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Write GFF coordinates and functional annotations into FASTA sequence headers.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "-gff", dest="gff",
        metavar="[augustus_prediction.gff]",
        required=True,
        help="Augustus prediction GFF file.",
    )
    parser.add_argument(
        "-annot", "--annotations", dest="annot",
        metavar="[annotations.txt]",
        required=True,
        help="Tab-separated annotations file.",
    )
    parser.add_argument(
        "-b", "--basename", dest="basename",
        metavar="[basename]",
        required=True,
        help="Basename for output files.",
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
        help="Organism name.",
    )
    return parser


def parse_gff_transcripts(gff_file: Path) -> dict[str, dict[str, str]]:
    """Extract transcript coordinates from GFF file."""
    transcripts = {}
    with open(gff_file, "r", encoding="utf-8", errors="replace") as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            parts = line.rstrip("\r\n").split("\t")
            if len(parts) >= 9 and parts[2].lower() in ("transcript", "mrna"):
                attr = re.sub(r";.*", "", parts[8].strip())
                tx_id = re.sub(r"^.*=", "", attr).strip()
                transcripts[tx_id] = {
                    "scaff": parts[0],
                    "start": parts[3],
                    "end": parts[4],
                    "strand": parts[6],
                }
    return transcripts


def parse_annotations(annot_file: Path) -> dict[str, str]:
    """Parse ID -> annotation mapping from file."""
    annotations = {}
    with open(annot_file, "r", encoding="utf-8", errors="replace") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t", 1)
            if len(parts) >= 2:
                annotations[parts[0].strip()] = parts[1].strip()
    return annotations


def annotate_fasta_with_gff(
    gff_file: str | Path,
    annot_file: str | Path,
    fasta_file: str | Path,
    basename: str,
    organism: str,
    output_dir: str | Path | None = None,
) -> Path:
    """Combine GFF locations and annotations into FASTA sequence headers."""
    gff_path = Path(gff_file)
    annot_path = Path(annot_file)
    fasta_path = Path(fasta_file)
    out_dir = Path(output_dir) if output_dir else Path.cwd()
    out_dir.mkdir(parents=True, exist_ok=True)

    org_clean = str(organism).strip("\"'")
    transcripts = parse_gff_transcripts(gff_path)
    annotations = parse_annotations(annot_path)

    out_file = out_dir / f"AnnotaPipeline_{basename}_proteins.fasta"
    ids_without_annot = []

    with open(fasta_path, "r", encoding="utf-8", errors="replace") as fasta_handle, \
         open(out_file, "w", encoding="utf-8") as out_handle:
        for record in SeqIO.parse(fasta_handle, "fasta"):
            seq_id = record.id.strip()
            tx_info = transcripts.get(seq_id, {
                "scaff": "unknown",
                "start": "0",
                "end": "0",
                "strand": ".",
            })
            desc = annotations.get(seq_id)

            if desc is not None:
                header = (
                    f"{seq_id} | Organism: {org_clean} | "
                    f"Location: {tx_info['scaff']} | "
                    f"Start: {tx_info['start']} | "
                    f"End: {tx_info['end']} | "
                    f"Strand: {tx_info['strand']} | "
                    f"Description: {desc}"
                )
                out_handle.write(f">{header}\n{str(record.seq).upper()}\n")
            else:
                ids_without_annot.append(seq_id)

    if ids_without_annot:
        warn_file = out_dir / f"{basename}_ids_with_no_annotations.txt"
        with open(warn_file, "w", encoding="utf-8") as wf:
            for item in ids_without_annot:
                wf.write(f"{item}\n")
        warnings.warn(
            f"WARNING: {len(ids_without_annot)} sequences lacked annotations. Stored in {warn_file}",
            stacklevel=2,
        )

    return out_file


def main(argv: list[str] | None = None) -> int:
    parser = cli()
    args = parser.parse_args(argv)
    annotate_fasta_with_gff(
        gff_file=args.gff,
        annot_file=args.annot,
        fasta_file=args.seqs,
        basename=args.basename,
        organism=args.org,
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())