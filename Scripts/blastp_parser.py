#!/usr/bin/env python3
"""
blastp_parser.py - Similarity analysis using SwissProt and a secondary database.

Executes BLASTp against SwissProt, filters annotations, extracts unannotated sequences,
and performs a secondary search against TrEMBL, NR, EuPathDB, or a custom database.
"""

import argparse
from collections import defaultdict
import logging
from pathlib import Path
import re
import shutil
import subprocess
import sys
from Bio import SeqIO


def cli() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Run and parse BLASTp searches against SwissProt and a secondary database.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    req = parser.add_argument_group("required arguments")
    req.add_argument(
        "-s", "--seq", dest="seq",
        metavar="[protein_file]",
        required=True,
        help="Input protein FASTA file.",
    )
    req.add_argument(
        "-sp", "--swissprot", dest="spdb",
        metavar="[SwissProt_db]",
        required=True,
        help="Path to SwissProt BLAST database.",
    )
    req.add_argument(
        "-basename", dest="basename",
        metavar="[basename]",
        required=True,
        help="Basename prefix for output files.",
    )

    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument(
        "-trbl", "--Trembl", dest="trembl",
        metavar="[Trembl_db]",
        help="Path to TrEMBL database.",
    )
    group.add_argument(
        "-spdb", "--specificdb", dest="specificdb",
        metavar="[Specific_db]",
        help="Path to specific database (e.g. EuPathDB format).",
    )
    group.add_argument(
        "-nr", "--nrdb", dest="nr",
        metavar="[NR_db]",
        help="Path to NR database.",
    )

    opt = parser.add_argument_group("optional arguments")
    opt.add_argument(
        "-id", "--identity", dest="id",
        type=float, default=40.0,
        help="Minimum percentage identity threshold.",
    )
    opt.add_argument(
        "-pos", "--positivity", dest="pos",
        type=float, default=60.0,
        help="Minimum percentage positivity threshold.",
    )
    opt.add_argument(
        "-cov", "--coverage", dest="cov",
        type=float, default=30.0,
        help="Minimum query coverage threshold.",
    )
    opt.add_argument(
        "-kw", "--keywords", dest="keywords",
        default="hypothetical,unspecified,fragment,partial,unknown,fragemnt",
        help="Comma-separated keywords indicating uninformative annotations.",
    )
    opt.add_argument(
        "-t", "--threads", dest="threads",
        type=int, default=8,
        help="Number of threads for BLASTp.",
    )
    opt.add_argument(
        "-max_target_seqs", dest="hsps",
        type=int, default=10,
        help="Maximum target sequences per query.",
    )
    opt.add_argument(
        "-evalue", dest="evalue",
        type=float, default=1e-5,
        help="E-value cutoff for BLASTp.",
    )
    opt.add_argument(
        "-blastp", dest="blastp",
        default="blastp",
        help="Path to blastp executable.",
    )
    opt.add_argument(
        "-customsep", dest="customsep",
        default="|",
        help="Field separator for custom database header.",
    )
    opt.add_argument(
        "-customcolumn", dest="customcolumn",
        type=int, default=5,
        help="Zero-indexed column containing annotation in custom database.",
    )
    return parser


def run_blast(
    blast_bin: str,
    query_file: Path,
    out_file: Path,
    db: str,
    threads: int,
    evalue: float,
    max_targets: int,
    logger: logging.Logger,
) -> None:
    """Execute blastp subprocess with standardized format 6 output."""
    fmt = "6 qseqid sseqid sacc bitscore evalue ppos pident qcovs stitle"
    cmd = [
        str(blast_bin),
        "-query", str(query_file),
        "-out", str(out_file),
        "-db", str(db),
        "-evalue", str(evalue),
        "-outfmt", fmt,
        "-max_target_seqs", str(max_targets),
        "-num_threads", str(threads),
    ]
    logger.info("Executing BLAST command: %s", " ".join(cmd))
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        logger.error("BLASTp execution failed: %s", result.stderr)
        raise RuntimeError(f"BLASTp failed with return code {result.returncode}: {result.stderr}")


def extract_description(
    title: str,
    dbtype: str,
    customsep: str = "|",
    customcolumn: int = 5,
) -> str:
    """Extract informative description string based on database format."""
    if dbtype == "nr":
        # Match '... [organism]'
        m = re.search(r"\s(.*?)\s\[", title)
        if m:
            return m.group(1).strip()
        return title.split(" ", 1)[-1].strip()

    if dbtype == "trembl" or dbtype == "swissprot":
        # Format: 'sp|P12345|ID Description OS=Species ...'
        desc_part = title.split(" ", 1)[-1] if " " in title else title
        return desc_part.split("OS=")[0].strip()

    if dbtype == "specificdb":
        # EuPathDB or custom format
        tokens = title.split(customsep)
        if len(tokens) > customcolumn:
            val = tokens[customcolumn].replace("transcript_product=", "").strip()
            return val
        return title.strip()

    return title.strip()


def parse_blast_records(outfmt6_file: Path) -> dict[str, list[dict]]:
    """Parse format 6 tabular output into hits grouped by query ID."""
    grouped_hits: dict[str, list[dict]] = defaultdict(list)
    if not outfmt6_file.exists() or outfmt6_file.stat().st_size == 0:
        return grouped_hits

    with open(outfmt6_file, "r", encoding="utf-8", errors="replace") as f:
        for line in f:
            parts = line.rstrip("\r\n").split("\t")
            if len(parts) < 9:
                continue
            try:
                grouped_hits[parts[0]].append({
                    "qseqid": parts[0],
                    "sseqid": parts[1],
                    "sacc": parts[2],
                    "bitscore": float(parts[3]),
                    "evalue": float(parts[4]),
                    "ppos": float(parts[5]),
                    "pident": float(parts[6]),
                    "qcovs": float(parts[7]),
                    "stitle": parts[8],
                })
            except (ValueError, IndexError):
                continue
    return grouped_hits


def process_swissprot_results(
    basename: str,
    query_fasta: Path,
    swiss_out: Path,
    min_identity: float,
    min_positivity: float,
    min_coverage: float,
    keywords: list[str],
    out_dir: Path,
) -> set[str]:
    """Filter SwissProt BLAST hits and export unannotated sequences."""
    grouped = parse_blast_records(swiss_out)
    annotated_file = out_dir / f"{basename}_annotated_products.txt"
    swiss_annot_file = out_dir / f"{basename}_SwissProt_annotations.txt"
    unannotated_fasta = out_dir / f"{basename}_BLASTp_AA_SwissProted.fasta"

    annotated_ids: set[str] = set()

    with open(annotated_file, "w", encoding="utf-8") as f_annot, \
         open(swiss_annot_file, "w", encoding="utf-8") as f_swiss:

        for qid, hits in grouped.items():
            valid_hits = []
            all_descriptions = []

            for hit in hits:
                desc = extract_description(hit["stitle"], "swissprot")
                all_descriptions.append(desc)

                has_kw = any(kw.lower() in desc.lower() for kw in keywords)
                if (hit["qcovs"] >= min_coverage and not has_kw and
                        hit["ppos"] >= min_positivity and hit["pident"] >= min_identity):
                    valid_hits.append((hit["bitscore"], desc))

            if valid_hits:
                # Select best bitscore hit
                valid_hits.sort(key=lambda x: x[0], reverse=True)
                best_desc = valid_hits[0][1]
                f_annot.write(f"{qid}\t{best_desc}\n")
                f_swiss.write(f"{qid}\t{len(all_descriptions)} Annotation(s): [{';'.join(all_descriptions)}]\n")
                annotated_ids.add(qid)

    # Write non-annotated queries to separate FASTA for secondary search
    unannotated_count = 0
    with open(unannotated_fasta, "w", encoding="utf-8") as f_unann:
        for record in SeqIO.parse(str(query_fasta), "fasta"):
            if record.id.strip() not in annotated_ids:
                SeqIO.write(record, f_unann, "fasta-2line")
                unannotated_count += 1

    return annotated_ids


def process_secondary_db_results(
    basename: str,
    secondary_out: Path,
    dbtype: str,
    min_identity: float,
    min_positivity: float,
    min_coverage: float,
    keywords: list[str],
    customsep: str,
    customcolumn: int,
    out_dir: Path,
) -> set[str]:
    """Parse secondary BLAST hits and classify into annotated vs hypothetical."""
    grouped = parse_blast_records(secondary_out)
    annotated_file = out_dir / f"{basename}_annotated_products.txt"
    hypo_file = out_dir / f"{basename}_hypothetical_products.txt"
    spec_annot_file = out_dir / f"{basename}_SpecifiedDB_annotations.txt"

    secondary_annotated_ids: set[str] = set()

    with open(annotated_file, "a", encoding="utf-8") as f_annot, \
         open(hypo_file, "w", encoding="utf-8") as f_hypo, \
         open(spec_annot_file, "w", encoding="utf-8") as f_spec:

        for qid, hits in grouped.items():
            valid_hits = []
            all_descriptions = []

            for hit in hits:
                desc = extract_description(hit["stitle"], dbtype, customsep, customcolumn)
                all_descriptions.append(desc)

                has_kw = any(kw.lower() in desc.lower() for kw in keywords)
                if (hit["qcovs"] >= min_coverage and not has_kw and
                        hit["ppos"] >= min_positivity and hit["pident"] >= min_identity):
                    valid_hits.append((hit["bitscore"], desc))

            if valid_hits:
                valid_hits.sort(key=lambda x: x[0], reverse=True)
                best_desc = valid_hits[0][1]
                f_annot.write(f"{qid}\t{best_desc}\n")
                f_spec.write(f"{qid}\t{len(all_descriptions)} Annotation(s): [{';'.join(all_descriptions)}]\n")
                secondary_annotated_ids.add(qid)
            else:
                f_hypo.write(f"{qid}\n")

    return secondary_annotated_ids


def identify_no_hits(
    basename: str,
    unannotated_fasta: Path,
    secondary_out: Path,
    out_dir: Path,
) -> Path:
    """Identify sequences with no alignments against SwissProt or secondary database."""
    hit_ids = set()
    if secondary_out.exists() and secondary_out.stat().st_size > 0:
        with open(secondary_out, "r", encoding="utf-8", errors="replace") as f:
            for line in f:
                parts = line.split("\t")
                if parts:
                    hit_ids.add(parts[0].strip())

    no_hit_file = out_dir / f"{basename}_no_hit_products.txt"
    with open(no_hit_file, "w", encoding="utf-8") as f_nohit:
        if unannotated_fasta.exists():
            for record in SeqIO.parse(str(unannotated_fasta), "fasta"):
                seq_id = record.id.strip()
                if seq_id not in hit_ids:
                    f_nohit.write(f"{seq_id}\n")

    return no_hit_file


def run_blastp_similarity(
    query_fasta: str | Path,
    swissprot_db: str,
    secondary_db: str,
    secondary_type: str,
    basename: str,
    blastp_bin: str = "blastp",
    threads: int = 8,
    evalue: float = 1e-5,
    max_targets: int = 10,
    min_identity: float = 40.0,
    min_positivity: float = 60.0,
    min_coverage: float = 30.0,
    keywords: list[str] | str | None = None,
    customsep: str = "|",
    customcolumn: int = 5,
    output_dir: str | Path | None = None,
    logger: logging.Logger | None = None,
) -> dict[str, Path]:
    """Complete modular similarity pipeline execution."""
    out_dir = Path(output_dir) if output_dir else Path.cwd()
    out_dir.mkdir(parents=True, exist_ok=True)
    log = logger or logging.getLogger("SimilarityAnalysis")

    if isinstance(keywords, str):
        kw_list = [w.strip() for w in keywords.split(",") if w.strip()]
    elif isinstance(keywords, list):
        kw_list = keywords
    else:
        kw_list = ["hypothetical", "unspecified", "fragment", "partial", "unknown", "fragemnt"]

    query_path = Path(query_fasta)
    swiss_out = out_dir / f"{basename}_BLASTp_AAvsSwissProt.outfmt6"
    secondary_out = out_dir / f"{basename}_BLASTp_AAvsSpecifiedDB.outfmt6"
    unannotated_fasta = out_dir / f"{basename}_BLASTp_AA_SwissProted.fasta"

    # 1. Search against SwissProt
    log.info("Running BLASTp against SwissProt...")
    run_blast(blastp_bin, query_path, swiss_out, swissprot_db, threads, evalue, max_targets, log)

    # 2. Process SwissProt hits
    log.info("Processing SwissProt annotations...")
    process_swissprot_results(
        basename, query_path, swiss_out, min_identity, min_positivity, min_coverage, kw_list, out_dir
    )

    # 3. Search against Secondary DB
    log.info("Running BLASTp against secondary database (%s)...", secondary_type)
    run_blast(blastp_bin, unannotated_fasta, secondary_out, secondary_db, threads, evalue, max_targets, log)

    # 4. Process Secondary DB hits
    log.info("Processing secondary database hits...")
    process_secondary_db_results(
        basename, secondary_out, secondary_type, min_identity, min_positivity, min_coverage,
        kw_list, customsep, customcolumn, out_dir
    )

    # 5. Extract no-hit queries
    log.info("Identifying queries with no hits in either database...")
    identify_no_hits(basename, unannotated_fasta, secondary_out, out_dir)

    log.info("Similarity analysis completed successfully.")
    return {
        "swissprot_out": swiss_out,
        "secondary_out": secondary_out,
        "annotated": out_dir / f"{basename}_annotated_products.txt",
        "hypothetical": out_dir / f"{basename}_hypothetical_products.txt",
        "no_hit": out_dir / f"{basename}_no_hit_products.txt",
    }


def main(argv: list[str] | None = None) -> int:
    parser = cli()
    args = parser.parse_args(argv)

    logging.basicConfig(level=logging.INFO, format="%(asctime)s | %(levelname)s | %(message)s")
    logger = logging.getLogger("SimilarityAnalysis")

    if args.nr:
        dbtype = "nr"
        secondary_db = args.nr
    elif args.trembl:
        dbtype = "trembl"
        secondary_db = args.trembl
    elif args.specificdb:
        dbtype = "specificdb"
        secondary_db = args.specificdb
    else:
        logger.error("No secondary database specified.")
        return 1

    try:
        run_blastp_similarity(
            query_fasta=args.seq,
            swissprot_db=args.spdb,
            secondary_db=secondary_db,
            secondary_type=dbtype,
            basename=args.basename,
            blastp_bin=args.blastp,
            threads=args.threads,
            evalue=args.evalue,
            max_targets=args.hsps,
            min_identity=args.id,
            min_positivity=args.pos,
            min_coverage=args.cov,
            keywords=args.keywords,
            customsep=args.customsep,
            customcolumn=args.customcolumn,
            logger=logger,
        )
    except Exception as exc:
        logger.error("BLAST pipeline failed: %s", exc, exc_info=True)
        return 1

    return 0


if __name__ == "__main__":
    sys.exit(main())