#!/usr/bin/env python3
"""
info_parser.py - Combine similarity annotations with InterPro and GO functional terms.

Integrates BLAST-derived annotations, hypothetical classifications, and no-hit sequences
with InterProScan terms into All_Annotated_Products.txt.
"""

import argparse
from pathlib import Path
import re
import sys

UNWANTED_DATABASES = {"coils", "gene3d", "mobidblite"}


def cli() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Combine BLAST annotations with InterProScan terms into All_Annotated_Products.txt.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    req = parser.add_argument_group("required arguments")
    req.add_argument(
        "-ipr1", "--ipr1", dest="ipr1",
        metavar="[annotated_ipr.gff3]",
        required=True,
        help="InterProScan output for annotated proteins.",
    )
    req.add_argument(
        "-ipr2", "--ipr2", dest="ipr2",
        metavar="[hypothetical_ipr.gff3]",
        required=True,
        help="InterProScan output for hypothetical proteins.",
    )
    req.add_argument(
        "-a", "--annotated", dest="annot",
        metavar="[annotated_products.txt]",
        required=True,
        help="Tab-separated annotated products from similarity analysis.",
    )
    req.add_argument(
        "-nh", "--no_hit", dest="nohit",
        metavar="[no_hit_products.txt]",
        required=True,
        help="Sequence IDs with no BLAST hits.",
    )
    req.add_argument(
        "-hy", "--hypothetical", dest="hypo",
        metavar="[hypothetical_products.txt]",
        required=True,
        help="Sequence IDs classified as hypothetical from BLAST.",
    )
    parser.add_argument(
        "-o", "--output", dest="output",
        default="All_Annotated_Products.txt",
        help="Output annotation file name.",
    )
    return parser


def natural_sort_key(s: str) -> list[int | str]:
    """Natural sort key for alphanumeric IDs (sort -V equivalent)."""
    return [int(text) if text.isdigit() else text.lower() for text in re.split(r"(\d+)", s)]


def parse_ipr_gff_terms(gff_file: Path) -> dict[str, tuple[set[str], set[str]]]:
    """Parse GFF3 file and return mapping of seq_id -> (set(ipr_terms), set(go_terms))."""
    terms_map: dict[str, tuple[set[str], set[str]]] = {}
    if not gff_file.exists() or gff_file.stat().st_size == 0:
        return terms_map

    with open(gff_file, "r", encoding="utf-8", errors="replace") as f:
        for line in f:
            line = line.strip()
            if line.startswith("##FASTA"):
                break
            if not line or line.startswith("#"):
                continue

            parts = line.split("\t")
            if len(parts) < 9:
                continue

            seq_id = parts[0].strip()
            db = parts[1].strip()
            if any(unwanted in db.lower() for unwanted in UNWANTED_DATABASES):
                continue

            ipr_set, go_set = terms_map.setdefault(seq_id, (set(), set()))

            attributes = parts[8].split(";")
            for attr in attributes:
                if attr.startswith("Ontology_term="):
                    vals = attr.replace("Ontology_term=", "").replace('"', "").split(",")
                    for v in vals:
                        v_clean = v.strip()
                        if v_clean and v_clean != "None":
                            go_set.add(v_clean)
                elif attr.startswith("Dbxref="):
                    vals = attr.replace("Dbxref=", "").replace('"', "").split(",")
                    for v in vals:
                        v_clean = v.strip()
                        if v_clean and v_clean != "None":
                            ipr_set.add(v_clean)

    return terms_map


def load_id_list(file_path: Path) -> list[str]:
    """Load line-separated IDs from a file."""
    if not file_path.exists() or file_path.stat().st_size == 0:
        return []
    with open(file_path, "r", encoding="utf-8", errors="replace") as f:
        return [line.strip().split("\t")[0] for line in f if line.strip() and not line.startswith("#")]


def load_annotated_products(file_path: Path) -> dict[str, str]:
    """Load annotated products mapping (ID -> description)."""
    annotated = {}
    if not file_path.exists() or file_path.stat().st_size == 0:
        return annotated
    with open(file_path, "r", encoding="utf-8", errors="replace") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t", 1)
            if len(parts) >= 2:
                annotated[parts[0].strip()] = parts[1].strip()
    return annotated


def _format_functional_suffix(ipr_terms: set[str], go_terms: set[str]) -> str:
    """Format (IPR..., GO...) suffix string."""
    sorted_ipr = sorted(ipr_terms, key=natural_sort_key)
    sorted_go = sorted(go_terms, key=natural_sort_key)

    if sorted_ipr and sorted_go:
        return f" ({','.join(sorted_ipr)},{','.join(sorted_go)})"
    if sorted_ipr:
        return f" ({','.join(sorted_ipr)})"
    if sorted_go:
        return f" ({','.join(sorted_go)})"
    return ""


def build_annotated_products(
    ipr1: str | Path,
    ipr2: str | Path,
    annot_file: str | Path,
    hypo_file: str | Path,
    nohit_file: str | Path,
    output_file: str | Path = "All_Annotated_Products.txt",
    output_dir: str | Path | None = None,
) -> Path:
    """Combine similarity annotations and functional terms, sorted by ID."""
    out_dir = Path(output_dir) if output_dir else Path.cwd()
    out_dir.mkdir(parents=True, exist_ok=True)
    out_path = out_dir / Path(output_file).name

    # 1. Parse InterProScan terms
    terms_map = parse_ipr_gff_terms(Path(ipr1))
    for seq_id, (ipr_set, go_set) in parse_ipr_gff_terms(Path(ipr2)).items():
        existing_ipr, existing_go = terms_map.setdefault(seq_id, (set(), set()))
        existing_ipr.update(ipr_set)
        existing_go.update(go_set)

    # 2. Load similarity categories
    annotated_dict = load_annotated_products(Path(annot_file))
    hypo_ids = set(load_id_list(Path(hypo_file)))
    nohit_ids = set(load_id_list(Path(nohit_file)))

    results: dict[str, str] = {}

    # Process annotated products
    for seq_id, base_desc in annotated_dict.items():
        ipr_terms, go_terms = terms_map.get(seq_id, (set(), set()))
        suffix = _format_functional_suffix(ipr_terms, go_terms)
        results[seq_id] = f"{base_desc}{suffix}"

    # Process hypothetical products
    for seq_id in hypo_ids:
        if seq_id not in results:
            ipr_terms, go_terms = terms_map.get(seq_id, (set(), set()))
            suffix = _format_functional_suffix(ipr_terms, go_terms)
            results[seq_id] = f"hypothetical protein{suffix}"

    # Process no-hit products
    for seq_id in nohit_ids:
        if seq_id not in results:
            ipr_terms, go_terms = terms_map.get(seq_id, (set(), set()))
            suffix = _format_functional_suffix(ipr_terms, go_terms)
            results[seq_id] = f"hypothetical protein{suffix}"

    # Any remaining proteins found only in InterProScan
    for seq_id, (ipr_terms, go_terms) in terms_map.items():
        if seq_id not in results:
            suffix = _format_functional_suffix(ipr_terms, go_terms)
            results[seq_id] = f"hypothetical protein{suffix}"

    # Write sorted output (natural version sort)
    sorted_ids = sorted(results.keys(), key=natural_sort_key)
    with open(out_path, "w", encoding="utf-8") as out:
        for seq_id in sorted_ids:
            out.write(f"{seq_id}\t{results[seq_id]}\n")

    return out_path


def main(argv: list[str] | None = None) -> int:
    parser = cli()
    args = parser.parse_args(argv)
    build_annotated_products(
        ipr1=args.ipr1,
        ipr2=args.ipr2,
        annot_file=args.annot,
        hypo_file=args.hypo,
        nohit_file=args.nohit,
        output_file=args.output,
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())