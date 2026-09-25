#!/usr/bin/env python3
"""
funcannotation_parser.py - Parse and integrate InterProScan, HMMER, and RPS-BLAST outputs.
"""

import argparse
import logging
from pathlib import Path
import re
import sys

UNWANTED_DATABASES = {"coils", "gene3d", "mobidblite"}


def cli() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Parse and integrate InterProScan, HMMER, and RPS-BLAST results.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    req = parser.add_argument_group("required arguments")
    req.add_argument(
        "-ipr_annot", dest="ipr_annot",
        metavar="[InterProScan_annotated.gff3]",
        required=True,
        help="InterProScan GFF3 for annotated proteins.",
    )
    req.add_argument(
        "-ipr_hyp", dest="ipr_hyp",
        metavar="[InterProScan_hypothetical.gff3]",
        required=True,
        help="InterProScan GFF3 for hypothetical proteins.",
    )
    req.add_argument(
        "-hmm", "--hmmscan", dest="hmm",
        metavar="[hmmscan_output.txt]",
        required=True,
        help="HMMSCAN tabular output file.",
    )
    req.add_argument(
        "-rpsblast", "--rpsblast", dest="rpsblast",
        metavar="[rpsblast_output.txt]",
        required=True,
        help="RPS-BLAST tabular output file.",
    )
    req.add_argument(
        "-basename", "--basename", dest="basename",
        metavar="[basename]",
        required=True,
        help="Basename prefix for output files.",
    )
    return parser


def parse_interpro_gff(gff_file: Path) -> list[tuple[str, str, str, str, str, str]]:
    """
    Parse InterProScan GFF3 output before ##FASTA.
    Returns list of tuples: (subject_id, db, name, description, ipr, ontology).
    """
    results = []
    if not gff_file.exists() or gff_file.stat().st_size == 0:
        return results

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

            seq_id = parts[0]
            db = parts[1]
            if any(unwanted in db.lower() for unwanted in UNWANTED_DATABASES):
                continue

            attributes = parts[8].split(";")
            ontology = "None"
            signature_desc = "None"
            name = "None"
            dbxref = "None"

            for attr in attributes:
                if attr.startswith("Ontology_term="):
                    ontology = attr.replace("Ontology_term=", "").strip('"')
                elif attr.startswith("signature_desc="):
                    signature_desc = attr.replace("signature_desc=", "").strip('"')
                elif attr.startswith("Name="):
                    name = attr.replace("Name=", "").strip('"')
                elif attr.startswith("Dbxref="):
                    dbxref = attr.replace("Dbxref=", "").strip('"')

            results.append((seq_id, db, name, signature_desc, dbxref, ontology))

    return results


def parse_hmmscan_tblout(hmm_file: Path) -> list[tuple[str, str, str, str, str, str]]:
    """Parse HMMER hmmscan --tblout tabular format."""
    results = []
    if not hmm_file.exists() or hmm_file.stat().st_size == 0:
        return results

    with open(hmm_file, "r", encoding="utf-8", errors="replace") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            # Regex split whitespace up to column 18
            parts = re.split(r"\s+", line, maxsplit=18)
            if len(parts) >= 19:
                target_name = parts[0]
                query_name = parts[2]
                description = parts[18]
                results.append((query_name, "Pfam", target_name, description, "None", "None"))
    return results


def parse_rpsblast_out(rps_file: Path) -> list[tuple[str, str, str, str, str, str]]:
    """Parse RPS-BLAST outfmt 6 format."""
    results = []
    if not rps_file.exists() or rps_file.stat().st_size == 0:
        return results

    with open(rps_file, "r", encoding="utf-8", errors="replace") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t")
            if len(parts) >= 9:
                query_id = parts[0]
                subject_acc = parts[1]
                title = parts[-1]
                results.append((query_id, "CDD", subject_acc, title, "None", "None"))
    return results


def parse_functional_annotations(
    ipr_annot: str | Path,
    ipr_hyp: str | Path,
    hmm_file: str | Path,
    rpsblast_file: str | Path,
    basename: str,
    output_dir: str | Path | None = None,
    logger: logging.Logger | None = None,
) -> dict[str, Path]:
    """Execute parsing of InterProScan, HMMER, and RPS-BLAST outputs."""
    out_dir = Path(output_dir) if output_dir else Path.cwd()
    out_dir.mkdir(parents=True, exist_ok=True)
    log = logger or logging.getLogger("FunctionalAnnotation")

    ipr_out_file = out_dir / f"InterProScan_Out_{basename}.tsv"
    hmm_out_file = out_dir / f"Hmmscan_Out_{basename}.tsv"
    rps_out_file = out_dir / f"RPSblast_Out_{basename}.tsv"
    grouped_file = out_dir / f"{basename}_Grouped_Hypothetical_Information.tsv"

    # Parse hypothetical proteins
    hyp_ipr_rows = parse_interpro_gff(Path(ipr_hyp))
    hmm_rows = parse_hmmscan_tblout(Path(hmm_file))
    rps_rows = parse_rpsblast_out(Path(rpsblast_file))

    # Parse annotated proteins
    annot_ipr_rows = parse_interpro_gff(Path(ipr_annot))

    # Write InterProScan combined output
    with open(ipr_out_file, "w", encoding="utf-8") as f:
        for r in hyp_ipr_rows + annot_ipr_rows:
            f.write("\t".join(r) + "\n")
    log.info("InterProScan output written to %s", ipr_out_file)

    # Write HMMER formatted output
    with open(hmm_out_file, "w", encoding="utf-8") as f:
        for r in hmm_rows:
            f.write("\t".join(r) + "\n")
    log.info("Hmmscan output written to %s", hmm_out_file)

    # Write RPS-BLAST output
    with open(rps_out_file, "w", encoding="utf-8") as f:
        for r in rps_rows:
            f.write("\t".join(r) + "\n")
    log.info("RPSblast output written to %s", rps_out_file)

    # Combine and sort hypothetical annotations
    hypothetical_combined = hyp_ipr_rows + hmm_rows + rps_rows
    hypothetical_combined.sort(key=lambda x: (x[0], x[1], x[2]))

    with open(grouped_file, "w", encoding="utf-8") as f:
        f.write("ID\tDB\tDB_ACCESS\tDESCRIPTION\tIPR\tGO\n")
        for r in hypothetical_combined:
            f.write("\t".join(r) + "\n")
    log.info("Grouped hypothetical information written to %s", grouped_file)

    return {
        "interpro": ipr_out_file,
        "hmmscan": hmm_out_file,
        "rpsblast": rps_out_file,
        "grouped_hypothetical": grouped_file,
    }


def main(argv: list[str] | None = None) -> int:
    parser = cli()
    args = parser.parse_args(argv)
    logging.basicConfig(level=logging.INFO, format="%(asctime)s | %(levelname)s | %(message)s")
    parse_functional_annotations(
        ipr_annot=args.ipr_annot,
        ipr_hyp=args.ipr_hyp,
        hmm_file=args.hmm,
        rpsblast_file=args.rpsblast,
        basename=args.basename,
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())