#!/usr/bin/env python3
"""
summary_parser.py - Aggregate similarity, functional, transcriptomic, and proteomic annotations.

Produces the comprehensive AnnotaPipeline_{basename}_Summary.tsv table.
"""

import argparse
from pathlib import Path
import re
import sys

UNWANTED_DATABASES = {"coils", "gene3d", "mobidblite"}


def cli() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Aggregate multi-omics annotations into a unified summary TSV.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    req = parser.add_argument_group("required arguments")
    req.add_argument(
        "-b", "--basename", dest="basename",
        metavar="[basename]",
        required=True,
        help="Pipeline basename.",
    )
    req.add_argument(
        "-annot", dest="annot",
        metavar="[All_Annotated_Products.txt]",
        required=True,
        help="All annotated products file.",
    )
    req.add_argument(
        "-ipr_annot", dest="ipr_annot",
        metavar="[ipr_annotated.gff3]",
        required=True,
        help="InterProScan output for annotated proteins.",
    )
    req.add_argument(
        "-ipr_hyp", dest="ipr_hyp",
        metavar="[ipr_hypothetical.gff3]",
        required=True,
        help="InterProScan output for hypothetical proteins.",
    )

    opt = parser.add_argument_group("optional arguments")
    opt.add_argument(
        "-tr", dest="tr",
        metavar="[transcript_quant.tsv]",
        default=None,
        help="Kallisto transcript quantification output TSV.",
    )
    opt.add_argument(
        "-proteomics", dest="proteomics",
        metavar="[proteomics_quant.tsv]",
        default=None,
        help="Proteomics quantification output TSV.",
    )
    return parser


def parse_interpro_gff_summary(
    gff_file: Path,
    superfamily_dict: dict[str, set[str]],
    ipr_dict: dict[str, set[str]],
    go_dict: dict[str, set[str]],
) -> None:
    """Parse InterProScan GFF3 to extract Superfamily, IPR, and GO terms per protein."""
    if not gff_file.exists() or gff_file.stat().st_size == 0:
        return

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

            superfamily_dict.setdefault(seq_id, set())
            ipr_dict.setdefault(seq_id, set())
            go_dict.setdefault(seq_id, set())

            attrs = parts[8].split(";")
            for attr in attrs:
                if attr.startswith("Ontology_term="):
                    terms = attr.replace("Ontology_term=", "").replace('"', "").split(",")
                    for t in terms:
                        t_clean = t.strip()
                        if t_clean and t_clean != "None":
                            go_dict[seq_id].add(t_clean)
                elif attr.startswith("Dbxref="):
                    terms = attr.replace("Dbxref=", "").replace('"', "").split(",")
                    for t in terms:
                        t_clean = t.strip()
                        if t_clean and t_clean != "None":
                            ipr_dict[seq_id].add(t_clean)
                elif "superfamily" in db.lower() and attr.startswith("Name="):
                    name_val = attr.replace("Name=", "").replace('"', "").strip()
                    if name_val and name_val != "None":
                        superfamily_dict[seq_id].add(name_val)


def parse_annotations(annot_file: Path) -> dict[str, str]:
    """Parse All_Annotated_Products.txt, extracting clean protein description."""
    annotations = {}
    if not annot_file.exists():
        return annotations

    with open(annot_file, "r", encoding="utf-8", errors="replace") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t", 1)
            if len(parts) >= 2:
                seq_id = parts[0].strip()
                # Remove appended (IPR..., InterPro..., GO...) to get base description
                base_desc = re.split(r"\s*\((?:IPR|InterPro|GO)", parts[1])[0].strip()
                annotations[seq_id] = base_desc
    return annotations


def parse_transcripts(tr_file: Path) -> dict[str, float]:
    """Parse transcript quantification TSV (target_id -> TPM)."""
    transcripts = {}
    if not tr_file or not tr_file.exists():
        return transcripts

    with open(tr_file, "r", encoding="utf-8", errors="replace") as f:
        header = f.readline()
        for line in f:
            line = line.strip()
            if not line:
                continue
            parts = line.split("\t")
            if len(parts) >= 2:
                try:
                    transcripts[parts[0].strip()] = round(float(parts[1]), 2)
                except ValueError:
                    continue
    return transcripts


def parse_proteomics(prot_file: Path) -> tuple[dict[str, str], dict[str, str]]:
    """Parse Total_Proteomics_Quantification.tsv to return (unique_peptides, total_peptides)."""
    unique_peptides = {}
    total_peptides = {}
    if not prot_file or not prot_file.exists():
        return unique_peptides, total_peptides

    with open(prot_file, "r", encoding="utf-8", errors="replace") as f:
        header = f.readline().split("\t")
        col_map = {name.strip(): idx for idx, name in enumerate(header)}

        pid_idx = col_map.get("ProteinID", 0)
        uniq_idx = col_map.get("Unique Peptide", 1 if len(header) > 1 else None)
        tot_idx = col_map.get("Total Peptide", 2 if len(header) > 2 else None)

        for line in f:
            line = line.strip()
            if not line:
                continue
            parts = line.split("\t")
            if len(parts) <= pid_idx:
                continue
            pid = parts[pid_idx].strip()
            if tot_idx is not None and len(parts) > tot_idx:
                total_peptides[pid] = parts[tot_idx].strip()
            if uniq_idx is not None and len(parts) > uniq_idx:
                val = parts[uniq_idx].strip()
                if val != "0":
                    unique_peptides[pid] = val

    return unique_peptides, total_peptides


def generate_summary(
    basename: str,
    annot_file: str | Path,
    ipr_annot: str | Path,
    ipr_hyp: str | Path,
    transcript_file: str | Path | None = None,
    proteomics_file: str | Path | None = None,
    output_dir: str | Path | None = None,
) -> Path:
    """Generate integrated summary TSV combining all analysis steps."""
    out_dir = Path(output_dir) if output_dir else Path.cwd()
    out_dir.mkdir(parents=True, exist_ok=True)

    superfamily_dict: dict[str, set[str]] = {}
    ipr_dict: dict[str, set[str]] = {}
    go_dict: dict[str, set[str]] = {}

    parse_interpro_gff_summary(Path(ipr_annot), superfamily_dict, ipr_dict, go_dict)
    parse_interpro_gff_summary(Path(ipr_hyp), superfamily_dict, ipr_dict, go_dict)

    annotations = parse_annotations(Path(annot_file))
    transcripts = parse_transcripts(Path(transcript_file)) if transcript_file else {}
    unique_pep, total_pep = parse_proteomics(Path(proteomics_file)) if proteomics_file else ({}, {})

    # Collect all unique protein IDs in deterministic order
    all_proteins = list(dict.fromkeys(
        list(annotations.keys()) +
        list(ipr_dict.keys()) +
        list(go_dict.keys()) +
        list(superfamily_dict.keys()) +
        list(transcripts.keys()) +
        list(total_pep.keys())
    ))

    summary_file = out_dir / f"AnnotaPipeline_{basename}_Summary.tsv"
    with open(summary_file, "w", encoding="utf-8") as out:
        out.write(
            "ProteinID\tAnnotation\tIPR\tGO\tSuperfamily\t"
            "Transcript\tTPM\tExpression\tTotal_Peptide\tUnique_Peptide\n"
        )

        for pid in all_proteins:
            annot_desc = annotations.get(pid, "None")

            ipr_val = ",".join(sorted(ipr_dict.get(pid, []))) or "None"
            ipr_val = ipr_val.replace("InterPro:", "").strip()

            go_val = ",".join(sorted(go_dict.get(pid, []))) or "None"
            sf_val = ",".join(sorted(superfamily_dict.get(pid, []))) or "None"

            has_tx = pid in transcripts
            tx_str = "True" if has_tx else "False"
            tpm_val = str(transcripts[pid]) if has_tx else "None"

            has_prot = pid in total_pep
            expr_str = "True" if has_prot else "False"
            tot_pep_val = total_pep.get(pid, "None")
            uniq_pep_val = unique_pep.get(pid, "None")

            out.write(
                f"{pid}\t{annot_desc}\t{ipr_val}\t{go_val}\t{sf_val}\t"
                f"{tx_str}\t{tpm_val}\t{expr_str}\t{tot_pep_val}\t{uniq_pep_val}\n"
            )

    return summary_file


def main(argv: list[str] | None = None) -> int:
    parser = cli()
    args = parser.parse_args(argv)
    generate_summary(
        basename=args.basename,
        annot_file=args.annot,
        ipr_annot=args.ipr_annot,
        ipr_hyp=args.ipr_hyp,
        transcript_file=args.tr,
        proteomics_file=args.proteomics,
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())