#!/usr/bin/env python3
"""
fastatogff.py - Integrate functional annotations into GFF files.

Parses AUGUSTUS GFF3 files and replaces transcript attributes with annotations
from All_Annotated_Products.txt.
"""

import argparse
from pathlib import Path
import re
import sys


def cli() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Integrate protein annotations into an AUGUSTUS GFF file.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "-all", "-a", "--annotations", "--allah", dest="anot",
        metavar="[All_Annotated_Products.txt]",
        required=True,
        help="AnnotaPipeline annotation file (e.g. All_Annotated_Products.txt).",
    )
    parser.add_argument(
        "-gff", "--gff", dest="gff",
        metavar="[prediction.gff]",
        required=True,
        help="GFF file containing AUGUSTUS predictions.",
    )
    parser.add_argument(
        "-b", "--basename", dest="basename",
        metavar="[basename]",
        required=True,
        help="Basename for output files.",
    )
    return parser


def load_annotations(annot_file: Path) -> dict[str, str]:
    """Load annotations into a dictionary mapping sequence ID to description."""
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


def extract_transcript_id(attribute_str: str) -> str:
    """Extract transcript/gene ID from the GFF attribute column."""
    # Remove ;.* to keep only first attribute if formatted as ID=...;Parent=...
    attr = re.sub(r";.*", "", attribute_str.strip())
    # Remove key= prefix (e.g. ID=)
    attr = re.sub(r"^.*=", "", attr)
    return attr.strip()


def annotate_gff(
    annot_file: str | Path,
    gff_file: str | Path,
    basename: str,
    output_dir: str | Path | None = None,
) -> Path:
    """Insert annotations into GFF transcript lines for matching IDs."""
    annot_path = Path(annot_file)
    gff_path = Path(gff_file)
    out_dir = Path(output_dir) if output_dir else Path.cwd()
    out_dir.mkdir(parents=True, exist_ok=True)

    annotations = load_annotations(annot_path)
    output_file = out_dir / f"{basename}_Annotated_GFF.gff"

    with open(gff_path, "r", encoding="utf-8", errors="replace") as gff_in, \
         open(output_file, "w", encoding="utf-8") as gff_out:

        header_lines = []
        current_block = []
        has_gene_started = False
        gene_transcript_id = None

        for line in gff_in:
            if line.startswith("# start gene "):
                if current_block:
                    _write_block_if_annotated(current_block, gene_transcript_id, annotations, gff_out)
                    current_block = []
                    gene_transcript_id = None
                has_gene_started = True
                current_block.append(line)
            elif not has_gene_started:
                header_lines.append(line)
                gff_out.write(line)
            else:
                current_block.append(line)
                if "\ttranscript\t" in line:
                    parts = line.rstrip("\r\n").split("\t")
                    if len(parts) >= 9:
                        gene_transcript_id = extract_transcript_id(parts[8])

        if current_block:
            _write_block_if_annotated(current_block, gene_transcript_id, annotations, gff_out)

    return output_file


def _write_block_if_annotated(
    block: list[str],
    transcript_id: str | None,
    annotations: dict[str, str],
    out_handle,
) -> None:
    """Write the gene block if transcript_id is annotated, replacing transcript attribute."""
    if transcript_id and transcript_id in annotations:
        annot_val = annotations[transcript_id]
        for line in block:
            if "\ttranscript\t" in line:
                parts = line.rstrip("\r\n").split("\t")
                if len(parts) >= 9:
                    # Replace 9th column (attributes) with annotation
                    parts[8] = annot_val
                    out_handle.write("\t".join(parts) + "\n")
                else:
                    out_handle.write(line)
            else:
                out_handle.write(line)


def main(argv: list[str] | None = None) -> int:
    parser = cli()
    args = parser.parse_args(argv)
    annotate_gff(
        annot_file=args.anot,
        gff_file=args.gff,
        basename=args.basename,
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())