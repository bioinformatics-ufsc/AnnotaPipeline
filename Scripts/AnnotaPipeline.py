#!/usr/bin/env python3
"""
AnnotaPipeline.py - Integrated Multi-Omics Protein Annotation Pipeline.

Coordinates gene prediction (AUGUSTUS), similarity analysis (BLASTp),
functional annotation (InterProScan, HMMER, RPS-BLAST), transcriptomics (Kallisto),
proteomics (Comet MS/MS, Percolator), and comprehensive summary generation.
"""

import argparse
import logging
from pathlib import Path
import re
import shutil
import subprocess
import sys
import yaml
from Bio import SeqIO
import pandas as pd

# Ensure repo root is available in sys.path for direct script execution
REPO_ROOT = Path(__file__).resolve().parent.parent
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

# Import modular pipeline components from Scripts
from Scripts.blastp_parser import run_blastp_similarity
from Scripts.funcannotation_parser import parse_functional_annotations
from Scripts.info_parser import build_annotated_products
from Scripts.fasta_simple import write_simple_fasta
from Scripts.fastatogff import annotate_gff
from Scripts.gfftofasta_parser import annotate_fasta_with_gff
from Scripts.kallisto_parser import parse_kallisto
from Scripts.percolator_parser import parse_percolator
from Scripts.summary_parser import generate_summary

STEPS = [
    "prediction",
    "similarity",
    "functional",
    "transcriptomics",
    "proteomics",
    "summary",
]


def cli() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="AnnotaPipeline: Integrated Multi-Omics Protein Annotation Pipeline.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    input_group = parser.add_mutually_exclusive_group(required=True)
    input_group.add_argument(
        "-s", "--seq", dest="seq",
        metavar="[genomic_fasta]",
        help="Input genomic sequence FASTA file (runs AUGUSTUS gene prediction).",
    )
    input_group.add_argument(
        "-p", "--prot", dest="protein",
        metavar="[protein_fasta]",
        help="Input protein sequence FASTA file (skips AUGUSTUS gene prediction).",
    )

    req = parser.add_argument_group("required arguments")
    req.add_argument(
        "-c", "--config", dest="config",
        metavar="[AnnotaPipeline.yaml]",
        required=True,
        help="Pipeline configuration YAML file.",
    )

    opt = parser.add_argument_group("optional arguments")
    opt.add_argument(
        "-gff", dest="gff",
        metavar="[prediction.gff]",
        help="Corresponding AUGUSTUS GFF3 file for the input protein FASTA.",
    )
    opt.add_argument(
        "--step", dest="single_step",
        choices=["all"] + STEPS,
        default="all",
        help="Execute only a single specified step.",
    )
    opt.add_argument(
        "--start-from", dest="start_from",
        choices=STEPS,
        default=None,
        help="Resume or start execution from this step onwards.",
    )
    opt.add_argument(
        "--resume", dest="resume",
        action="store_true",
        help="Skip steps whose output files already exist.",
    )
    opt.add_argument(
        "--force", dest="force",
        action="store_true",
        help="Force rerun of all steps, ignoring existing outputs.",
    )
    return parser


def setup_logger(log_file: Path | None = None) -> logging.Logger:
    logger = logging.getLogger("AnnotaPipeline")
    logger.setLevel(logging.DEBUG)
    logger.handlers.clear()

    formatter = logging.Formatter("%(asctime)s | %(name)s | %(levelname)s | %(message)s", datefmt="%d/%m/%Y %H:%M:%S")

    # Console Handler
    ch = logging.StreamHandler(sys.stdout)
    ch.setLevel(logging.INFO)
    ch.setFormatter(formatter)
    logger.addHandler(ch)

    # File Handler
    if log_file:
        fh = logging.FileHandler(log_file, mode="a", encoding="utf-8")
        fh.setLevel(logging.DEBUG)
        fh.setFormatter(formatter)
        logger.addHandler(fh)

    return logger


def check_tool_exists(tool_name: str, logger: logging.Logger) -> bool:
    """Check whether a binary is available on $PATH."""
    if shutil.which(tool_name) is None:
        logger.warning("Tool '%s' was not found on PATH.", tool_name)
        return False
    return True


def sequence_cleaner(fasta_file: Path, min_length: int = 0, output_file: Path | None = None) -> Path:
    """Filter sequences in FASTA by minimum length."""
    out_path = output_file or fasta_file.parent / f"Clear_{fasta_file.name}"
    with open(out_path, "w", encoding="utf-8") as out_handle:
        for record in SeqIO.parse(str(fasta_file), "fasta"):
            seq_str = str(record.seq).upper()
            if len(seq_str) >= min_length:
                out_handle.write(f">{record.id}\n{seq_str}\n")
    return out_path


def fasta_fetcher(input_fasta: Path, wanted_ids: set[str], output_fasta: Path) -> int:
    """Extract specific sequences by ID into a new FASTA file."""
    count = 0
    with open(output_fasta, "w", encoding="utf-8") as out:
        for record in SeqIO.parse(str(input_fasta), "fasta"):
            rec_id = record.id.strip()
            if rec_id in wanted_ids or any(w in rec_id for w in wanted_ids):
                SeqIO.write(record, out, "fasta-2line")
                count += 1
    return count


def annotate_codingseq(aa_fasta: Path, codingseq_fasta: Path, output_file: Path, logger: logging.Logger) -> Path:
    """Propagate updated protein headers onto the coding sequence (.codingseq) transcripts."""
    logger.info("Transferring annotations to transcript coding sequences: %s", output_file.name)
    headers = {}
    with open(aa_fasta, "r", encoding="utf-8", errors="replace") as f:
        for line in f:
            if line.startswith(">"):
                full = line.strip()[1:]
                qid = full.split("|")[0].strip()
                headers[qid] = full

    with open(output_file, "w", encoding="utf-8") as out:
        for record in SeqIO.parse(str(codingseq_fasta), "fasta"):
            matched_header = None
            for qid, hdr in headers.items():
                if qid in record.id:
                    matched_header = hdr
                    break
            if matched_header:
                record.id = matched_header
                record.description = ""
            record.seq = record.seq.upper()
            SeqIO.write(record, out, "fasta-2line")

    return output_file


def modify_comet_params(comet_params_path: Path, logger: logging.Logger) -> None:
    """Ensure required Comet configuration flags are set."""
    if not comet_params_path.exists():
        logger.warning("Comet params file not found at %s", comet_params_path)
        return

    text = comet_params_path.read_text(encoding="utf-8", errors="replace")
    text = re.sub(r"decoy_search\s*=\s*[0-9]", "decoy_search = 1", text)
    text = re.sub(r"output_pepxmlfile\s*=\s*[0-9]", "output_pepxmlfile = 0", text)
    text = re.sub(r"output_percolatorfile\s*=\s*[0-9]", "output_percolatorfile = 1", text)
    text = re.sub(r"decoy_prefix\s*=\s*\S+", "decoy_prefix = DECOY_", text)
    comet_params_path.write_text(text, encoding="utf-8")


def quantitative_proteomics(parsed_dir: Path, basename: str, out_dir: Path) -> Path:
    """Aggregate parsed Percolator results to calculate unique/total peptides and spectrums."""
    parsed_files = list(parsed_dir.glob("*_parsed.tsv"))
    dfs = []
    for pf in parsed_files:
        try:
            df = pd.read_csv(pf, sep="\t")
            if not df.empty and "ProteinID" in df.columns:
                dfs.append(df)
        except Exception:
            continue

    total_out = out_dir / f"{basename}_Total_Proteomics_Quantification.tsv"
    if not dfs:
        # Create empty quantification table
        with open(total_out, "w", encoding="utf-8") as f:
            f.write("ProteinID\tUnique Peptide\tTotal Peptide\tUnique Spectrum\tTotal Spectrum\n")
        return total_out

    data = pd.concat(dfs, ignore_index=True)

    def compute_counts(feature: str):
        df_feat = data[["ProteinID", feature]].dropna()
        # Count totals
        total_counts = df_feat.groupby("ProteinID").size().rename(f"Total {feature}")
        # Count unique occurrences across entire dataset
        val_counts = df_feat[feature].value_counts()
        unique_vals = set(val_counts[val_counts == 1].index)
        unique_counts = df_feat[df_feat[feature].isin(unique_vals)].groupby("ProteinID").size().rename(f"Unique {feature}")
        return pd.concat([unique_counts, total_counts], axis=1).fillna(0).astype(int)

    pep_df = compute_counts("Peptide")
    spec_df = compute_counts("Spectrum")
    combined = pd.concat([pep_df, spec_df], axis=1).fillna(0).astype(int).reset_index()

    combined.sort_values(by="ProteinID", inplace=True)
    combined.to_csv(total_out, sep="\t", index=False)
    return total_out


class AnnotaPipeline:
    """Orchestrator for AnnotaPipeline steps."""

    def __init__(self, config_path: Path, genome_file: Path | None, protein_file: Path | None, gff_file: Path | None):
        self.config_path = config_path.resolve()
        with open(self.config_path, "r", encoding="utf-8") as stream:
            self.config = yaml.safe_load(stream)

        self.genome_file = genome_file.resolve() if genome_file else None
        self.protein_file = protein_file.resolve() if protein_file else None
        self.gff_file = gff_file.resolve() if gff_file else None

        pipeline_cfg = self.config.get("pipeline", {})
        self.basename = pipeline_cfg.get("basename", "Annota")
        self.organism = pipeline_cfg.get("organism", "Unknown")
        self.threads = int(pipeline_cfg.get("threads", 8))
        self.keywords = pipeline_cfg.get("keywords", ["hypothetical", "unspecified", "fragment", "partial"])

        # Directory structure
        self.root_dir = Path.cwd() / f"AnnotaPipeline_{self.basename}"
        self.root_dir.mkdir(parents=True, exist_ok=True)

        self.pred_dir = self.root_dir / f"1_GenePrediction_{self.basename}"
        self.sim_dir = self.root_dir / f"2_SimilarityAnalysis_{self.basename}"
        self.func_dir = self.root_dir / f"3_FunctionalAnnotation_{self.basename}"
        self.tx_dir = self.root_dir / f"4_TranscriptQuantification_{self.basename}"
        self.prot_dir = self.root_dir / f"5_PeptideIdentification_{self.basename}"

        self.logger = setup_logger(self.root_dir / "AnnotaPipeline_LogFile.log")

    def run_step_prediction(self, resume: bool = False, force: bool = False) -> Path:
        """Step 1: Gene prediction via AUGUSTUS or ingest provided protein sequences."""
        self.pred_dir.mkdir(parents=True, exist_ok=True)
        self.logger.info("===============================================================")
        self.logger.info("STEP 1: Gene Prediction / Protein Sequence Preparation")
        self.logger.info("===============================================================")

        min_size = int(self.config.get("seq-cleaner", {}).get("minsize-seq", 0))
        cleaned_aa = self.pred_dir / f"Clear_{self.basename}.aa"

        if resume and not force and cleaned_aa.exists() and cleaned_aa.stat().st_size > 0:
            self.logger.info("Cleaned protein sequences already exist at %s. Skipping Step 1.", cleaned_aa)
            return cleaned_aa

        if self.protein_file:
            self.logger.info("Using provided protein FASTA: %s", self.protein_file)
            target_aa = self.pred_dir / f"{self.basename}.aa"
            shutil.copy2(self.protein_file, target_aa)
            sequence_cleaner(target_aa, min_length=min_size, output_file=cleaned_aa)
        else:
            self.logger.info("Running AUGUSTUS prediction on genome: %s", self.genome_file)
            aug_cfg = self.config.get("augustus", {})
            aug_opt = self.config.get("augustus-optional", {})

            aug_path = aug_cfg.get("augustus-path", "conda")
            if str(aug_path).lower() == "conda":
                aug_bin = "augustus"
                aug_script = "getAnnoFasta.pl"
            else:
                aug_bin = str(Path(aug_path) / "bin" / "augustus")
                aug_script = str(Path(aug_path) / "scripts" / "getAnnoFasta.pl")

            aug_gff = self.pred_dir / f"AUGUSTUS_{self.basename}.gff"
            cmd = [aug_bin]
            for param, val in aug_cfg.items():
                if param == "augustus-path":
                    continue
                if str(val).lower() == "flag":
                    cmd.append(f"-{param}")
                else:
                    cmd.append(f"--{param}={val}")
            for param, val in aug_opt.items():
                if val is not None:
                    if str(val).lower() == "flag":
                        cmd.append(f"-{param}")
                    else:
                        cmd.append(f"--{param}={val}")
            cmd.append(str(self.genome_file))

            self.logger.info("Executing AUGUSTUS: %s", " ".join(cmd))
            with open(aug_gff, "w", encoding="utf-8") as f_gff:
                subprocess.run(cmd, stdout=f_gff, check=True)

            # Extract FASTA files using getAnnoFasta.pl
            cmd_extract = [aug_script, str(aug_gff), f"--seqfile={self.genome_file}"]
            self.logger.info("Extracting predicted proteins with getAnnoFasta.pl...")
            subprocess.run(cmd_extract, cwd=str(self.pred_dir), check=True)

            aug_aa = self.pred_dir / f"AUGUSTUS_{self.basename}.aa"
            sequence_cleaner(aug_aa, min_length=min_size, output_file=cleaned_aa)

        self.logger.info("Gene prediction / sequence cleaning completed: %s", cleaned_aa)
        return cleaned_aa

    def run_step_similarity(self, resume: bool = False, force: bool = False) -> dict[str, Path]:
        """Step 2: Similarity search against SwissProt and secondary DB."""
        self.sim_dir.mkdir(parents=True, exist_ok=True)
        self.logger.info("===============================================================")
        self.logger.info("STEP 2: Similarity Analysis (BLASTp)")
        self.logger.info("===============================================================")

        annot_file = self.sim_dir / f"{self.basename}_annotated_products.txt"
        hypo_file = self.sim_dir / f"{self.basename}_hypothetical_products.txt"
        nohit_file = self.sim_dir / f"{self.basename}_no_hit_products.txt"

        if resume and not force and annot_file.exists() and hypo_file.exists() and nohit_file.exists():
            self.logger.info("Similarity analysis outputs exist. Skipping Step 2.")
            return {"annotated": annot_file, "hypothetical": hypo_file, "no_hit": nohit_file}

        cleaned_aa = self.pred_dir / f"Clear_{self.basename}.aa"
        if not cleaned_aa.exists():
            cleaned_aa = self.run_step_prediction(resume=resume, force=force)

        db_cfg = self.config.get("databases", {})
        blast_cfg = self.config.get("local-aligner", {})

        secondary_format = str(db_cfg.get("secondary-format", "nrdb")).lower()
        if secondary_format == "nr":
            secondary_format = "nrdb"
        elif secondary_format in ("eupathdb", "custom"):
            secondary_format = "specificdb"

        results = run_blastp_similarity(
            query_fasta=cleaned_aa,
            swissprot_db=db_cfg.get("swissprot-db"),
            secondary_db=db_cfg.get("secondary-db"),
            secondary_type=secondary_format,
            basename=self.basename,
            blastp_bin=blast_cfg.get("blastp", "blastp"),
            threads=self.threads,
            evalue=float(blast_cfg.get("evalue", 1e-5)),
            max_targets=int(blast_cfg.get("max_target_seqs", 10)),
            min_identity=float(blast_cfg.get("identity", 40.0)),
            min_positivity=float(blast_cfg.get("positivity", 60.0)),
            min_coverage=float(blast_cfg.get("coverage", 30.0)),
            keywords=self.keywords,
            customsep=db_cfg.get("customsep", "|") or "|",
            customcolumn=int(db_cfg.get("customcolumn", 5) or 5),
            output_dir=self.sim_dir,
            logger=self.logger,
        )
        return results

    def run_step_functional(self, resume: bool = False, force: bool = False) -> Path:
        """Step 3: Functional annotation (InterProScan, HMMER, RPS-BLAST) and file assembly."""
        self.func_dir.mkdir(parents=True, exist_ok=True)
        self.logger.info("===============================================================")
        self.logger.info("STEP 3: Functional Annotation & Final Sequence Assembly")
        self.logger.info("===============================================================")

        final_annot = self.root_dir / "All_Annotated_Products.txt"
        final_prot_fasta = self.root_dir / f"AnnotaPipeline_{self.basename}_proteins.fasta"

        if resume and not force and final_annot.exists() and final_prot_fasta.exists():
            self.logger.info("Functional annotation products already exist. Skipping Step 3.")
            return final_annot

        # Ensure prerequisites from Step 1 and 2
        cleaned_aa = self.pred_dir / f"Clear_{self.basename}.aa"
        if not cleaned_aa.exists():
            cleaned_aa = self.run_step_prediction(resume=resume, force=force)

        sim_annot = self.sim_dir / f"{self.basename}_annotated_products.txt"
        sim_hypo = self.sim_dir / f"{self.basename}_hypothetical_products.txt"
        sim_nohit = self.sim_dir / f"{self.basename}_no_hit_products.txt"
        if not (sim_annot.exists() and sim_hypo.exists() and sim_nohit.exists()):
            self.run_step_similarity(resume=resume, force=force)

        # 1. Prepare Annotated & Hypothetical FASTA subsets
        annot_ids = {line.split("\t")[0].strip() for line in open(sim_annot) if line.strip()}
        hypo_ids = {line.strip() for line in open(sim_hypo) if line.strip()}
        nohit_ids = {line.strip() for line in open(sim_nohit) if line.strip()}
        all_hypo_ids = hypo_ids | nohit_ids

        annot_fasta = self.func_dir / "Annotated_Products.fasta"
        hypo_fasta = self.func_dir / "Hypothetical_Products.fasta"

        fasta_fetcher(cleaned_aa, annot_ids, annot_fasta)
        fasta_fetcher(cleaned_aa, all_hypo_ids, hypo_fasta)

        # 2. InterProScan
        ipr_cfg = self.config.get("interproscan", {})
        ipr_annot_gff = self.func_dir / f"{self.basename}_interproscan_annotated_output.gff3"
        ipr_hypo_gff = self.func_dir / f"{self.basename}_interproscan_hypothetical_output.gff3"

        def _run_ipr(in_fasta: Path, out_gff: Path):
            if in_fasta.exists() and in_fasta.stat().st_size > 0:
                cmd = ["interproscan.sh", "-i", str(in_fasta), "-o", str(out_gff), "-f", "GFF3", "-t", "p", "-goterms", "-iprlookup"]
                for k, v in ipr_cfg.items():
                    if v is not None:
                        if str(v).lower() == "flag":
                            cmd.append(f"-{k}")
                        else:
                            cmd.extend([f"-{k}", str(v)])
                self.logger.info("Executing InterProScan: %s", " ".join(cmd))
                subprocess.run(cmd, check=False)
            if not out_gff.exists():
                out_gff.touch()

        if not (resume and not force and ipr_annot_gff.exists()):
            _run_ipr(annot_fasta, ipr_annot_gff)
        if not (resume and not force and ipr_hypo_gff.exists()):
            _run_ipr(hypo_fasta, ipr_hypo_gff)

        # 3. HMMER
        hmm_out = self.func_dir / f"{self.basename}_hmmscan_output.txt"
        pfam_db = self.config.get("databases", {}).get("pfam")
        hmm_cfg = self.config.get("hmmer", {})
        if not (resume and not force and hmm_out.exists()) and hypo_fasta.exists() and hypo_fasta.stat().st_size > 0 and pfam_db:
            cmd_hmm = ["hmmscan", "--cpu", str(self.threads), "--tblout", str(hmm_out), "--noali"]
            for k, v in hmm_cfg.items():
                if v is not None:
                    if k in ("E", "Z", "T"):
                        cmd_hmm.extend([f"-{k}", str(v)])
                    else:
                        cmd_hmm.extend([f"--{k}", str(v)])
            cmd_hmm.extend([str(pfam_db), str(hypo_fasta)])
            self.logger.info("Executing HMMER hmmscan: %s", " ".join(cmd_hmm))
            subprocess.run(cmd_hmm, check=False)
        if not hmm_out.exists():
            hmm_out.touch()

        # 4. RPS-BLAST
        rps_out = self.func_dir / f"{self.basename}_rpsblast_output.outfmt6"
        cdd_db = self.config.get("databases", {}).get("cdd-db")
        rps_cfg = self.config.get("rpsblast", {})
        if not (resume and not force and rps_out.exists()) and hypo_fasta.exists() and hypo_fasta.stat().st_size > 0 and cdd_db:
            cmd_rps = [
                "rpsblast", "-query", str(hypo_fasta), "-out", str(rps_out), "-db", str(cdd_db),
                "-outfmt", "6 qseqid sseqid sacc bitscore evalue ppos pident qcovs stitle",
                "-num_threads", str(self.threads),
            ]
            for k, v in rps_cfg.items():
                if v is not None:
                    cmd_rps.extend([f"-{k}", str(v)])
            self.logger.info("Executing RPS-BLAST: %s", " ".join(cmd_rps))
            subprocess.run(cmd_rps, check=False)
        if not rps_out.exists():
            rps_out.touch()

        # 5. Parse Functional Annotations
        self.logger.info("Parsing InterProScan, HMMER, and RPS-BLAST outputs...")
        parse_functional_annotations(
            ipr_annot=ipr_annot_gff,
            ipr_hyp=ipr_hypo_gff,
            hmm_file=hmm_out,
            rpsblast_file=rps_out,
            basename=self.basename,
            output_dir=self.func_dir,
            logger=self.logger,
        )

        # 6. Build All_Annotated_Products.txt in root_dir
        self.logger.info("Integrating terms into All_Annotated_Products.txt...")
        build_annotated_products(
            ipr1=ipr_annot_gff,
            ipr2=ipr_hypo_gff,
            annot_file=sim_annot,
            hypo_file=sim_hypo,
            nohit_file=sim_nohit,
            output_file=final_annot,
            output_dir=self.root_dir,
        )

        # 7. Generate Annotated FASTA & GFF outputs
        gff_source = self.gff_file or (self.pred_dir / f"AUGUSTUS_{self.basename}.gff")
        if gff_source and gff_source.exists():
            self.logger.info("Annotating FASTA with genomic coordinates from GFF: %s", gff_source)
            annotate_fasta_with_gff(
                gff_file=gff_source,
                annot_file=final_annot,
                fasta_file=cleaned_aa,
                basename=self.basename,
                organism=self.organism,
                output_dir=self.root_dir,
            )
            self.logger.info("Annotating GFF file...")
            annotate_gff(
                annot_file=final_annot,
                gff_file=gff_source,
                basename=self.basename,
                output_dir=self.root_dir,
            )
        else:
            self.logger.info("No GFF provided; generating simple annotated FASTA.")
            write_simple_fasta(
                annot_file=final_annot,
                fasta_file=cleaned_aa,
                basename=self.basename,
                organism=self.organism,
                output_dir=self.root_dir,
            )

        # 8. If AUGUSTUS codingseq exists, annotate transcripts
        codingseq_source = self.pred_dir / f"AUGUSTUS_{self.basename}.codingseq"
        if codingseq_source.exists():
            transcripts_fasta = self.root_dir / f"AnnotaPipeline_{self.basename}_transcripts.fasta"
            annotate_codingseq(final_prot_fasta, codingseq_source, transcripts_fasta, self.logger)

        self.logger.info("Functional annotation step completed.")
        return final_annot

    def run_step_transcriptomics(self, resume: bool = False, force: bool = False) -> Path | None:
        """Step 4: Transcript quantification using Kallisto."""
        kt_cfg = self.config.get("kallisto", {})
        if not kt_cfg.get("run-kallisto", False):
            self.logger.info("Transcriptomics (Kallisto) disabled in config. Skipping Step 4.")
            return None

        self.tx_dir.mkdir(parents=True, exist_ok=True)
        self.logger.info("===============================================================")
        self.logger.info("STEP 4: Transcript Quantification (Kallisto)")
        self.logger.info("===============================================================")

        tx_out = self.tx_dir / f"{self.basename}_Transcript_Quantification.tsv"
        if resume and not force and tx_out.exists():
            self.logger.info("Transcript quantification output exists at %s. Skipping Step 4.", tx_out)
            return tx_out

        transcripts_fasta = self.root_dir / f"AnnotaPipeline_{self.basename}_transcripts.fasta"
        if not transcripts_fasta.exists():
            self.logger.warning("Transcripts FASTA not found (%s). Cannot run Kallisto.", transcripts_fasta)
            return None

        idx_file = self.tx_dir / f"{self.basename}_kallisto_index.idx"
        if not idx_file.exists():
            self.logger.info("Generating Kallisto index: %s", idx_file)
            subprocess.run(["kallisto", "index", "-i", str(idx_file), str(transcripts_fasta)], check=True)

        reads = kt_cfg.get("rna-seq", [])
        if not isinstance(reads, list):
            reads = [reads]
        reads_str = [str(r) for r in reads]

        kallisto_out_folder = self.tx_dir / f"{self.basename}_kallisto_output"
        cmd_quant = ["kallisto", "quant", "-i", str(idx_file), "-o", str(kallisto_out_folder), "-b", str(kt_cfg.get("bootstrap", 100))]

        if len(reads) == 1:
            cmd_quant.extend(["--single", "-l", str(kt_cfg.get("l", 200)), "-s", str(kt_cfg.get("s", 20))])
        else:
            if kt_cfg.get("l"):
                cmd_quant.extend(["-l", str(kt_cfg.get("l"))])
            if kt_cfg.get("s"):
                cmd_quant.extend(["-s", str(kt_cfg.get("s"))])

        cmd_quant.extend(reads_str)
        self.logger.info("Executing Kallisto quant: %s", " ".join(cmd_quant))
        subprocess.run(cmd_quant, check=True)

        abundance_file = kallisto_out_folder / "abundance.tsv"
        th_val = kt_cfg.get("threshold", "median")
        if th_val in ("mean", "median"):
            method = th_val
            num_val = None
        else:
            method = "value"
            num_val = float(th_val)

        parse_kallisto(
            kallisto_file=abundance_file,
            basename=self.basename,
            method=method,
            threshold_value=num_val,
            output_dir=self.tx_dir,
        )
        self.logger.info("Transcript quantification completed: %s", tx_out)
        return tx_out

    def run_step_proteomics(self, resume: bool = False, force: bool = False) -> Path | None:
        """Step 5: Peptide identification via Comet MS/MS and Percolator."""
        prot_cfg = self.config.get("proteomics", {})
        if not prot_cfg.get("run-comet", False):
            self.logger.info("Proteomics (Comet) disabled in config. Skipping Step 5.")
            return None

        self.prot_dir.mkdir(parents=True, exist_ok=True)
        self.logger.info("===============================================================")
        self.logger.info("STEP 5: Peptide Identification (Comet MS/MS & Percolator)")
        self.logger.info("===============================================================")

        total_quant = self.prot_dir / f"{self.basename}_Total_Proteomics_Quantification.tsv"
        if resume and not force and total_quant.exists():
            self.logger.info("Proteomics quantification output exists. Skipping Step 5.")
            return total_quant

        prot_fasta = self.root_dir / f"AnnotaPipeline_{self.basename}_proteins.fasta"
        comet_params = Path(prot_cfg.get("comet-params"))
        modify_comet_params(comet_params, self.logger)

        spec_dir = Path(prot_cfg.get("comet-spectrometry"))
        spec_ext = prot_cfg.get("comet-ext", "mzXML")
        spec_files = list(spec_dir.glob(f"*.{spec_ext}"))

        if not spec_files:
            self.logger.warning("No mass spectrometry files found with extension .%s in %s", spec_ext, spec_dir)
            return None

        # Execute Comet
        cmd_comet = ["comet", f"-P{comet_params}", f"-D{prot_fasta}"]
        if prot_cfg.get("first") is not None and prot_cfg.get("last") is not None:
            cmd_comet.extend([f"-F{prot_cfg.get('first')}", f"-L{prot_cfg.get('last')}"])
        cmd_comet.extend([str(f) for f in spec_files])

        self.logger.info("Running Comet MS/MS...")
        subprocess.run(cmd_comet, cwd=str(self.prot_dir), check=False)

        # Process .pin files with Percolator
        raw_dir = self.prot_dir / "PERCOLATOR_Raw"
        parsed_dir = self.prot_dir / "PERCOLATOR_Parsed"
        raw_dir.mkdir(exist_ok=True)
        parsed_dir.mkdir(exist_ok=True)

        pin_files = list(spec_dir.glob("*.pin")) + list(self.prot_dir.glob("*.pin"))
        qval = float(prot_cfg.get("percolator-qvalue", 0.05))

        for pin in pin_files:
            pin_stem = re.sub(r"\.[0-9].*", "", pin.stem)
            perc_out = raw_dir / f"{pin_stem}_percolator_output.tsv"

            cmd_perc = ["percolator", "-r", str(raw_dir / f"{pin_stem}_peptide.tsv"), "-m", str(perc_out), str(pin)]
            self.logger.info("Running Percolator on %s...", pin.name)
            subprocess.run(cmd_perc, check=False)

            if perc_out.exists():
                parse_percolator(
                    percolator_file=perc_out,
                    basename=f"{self.basename}_{pin_stem}",
                    qvalue=qval,
                    output_dir=parsed_dir,
                )

        # Aggregate unique and total counts
        self.logger.info("Aggregating quantitative proteomics counts...")
        quantitative_proteomics(parsed_dir, self.basename, self.prot_dir)
        return total_quant

    def run_step_summary(self, resume: bool = False, force: bool = False) -> Path:
        """Step 6: Aggregate multi-omics results into final TSV summary table."""
        self.logger.info("===============================================================")
        self.logger.info("STEP 6: Aggregating Summary Report")
        self.logger.info("===============================================================")

        summary_file = self.root_dir / f"AnnotaPipeline_{self.basename}_Summary.tsv"
        if resume and not force and summary_file.exists():
            self.logger.info("Summary file already exists: %s", summary_file)
            return summary_file

        annot_file = self.root_dir / "All_Annotated_Products.txt"
        ipr_annot = self.func_dir / f"{self.basename}_interproscan_annotated_output.gff3"
        ipr_hyp = self.func_dir / f"{self.basename}_interproscan_hypothetical_output.gff3"

        tx_file = self.tx_dir / f"{self.basename}_Transcript_Quantification.tsv"
        prot_file = self.prot_dir / f"{self.basename}_Total_Proteomics_Quantification.tsv"

        out_summary = generate_summary(
            basename=self.basename,
            annot_file=annot_file,
            ipr_annot=ipr_annot,
            ipr_hyp=ipr_hyp,
            transcript_file=tx_file if tx_file.exists() else None,
            proteomics_file=prot_file if prot_file.exists() else None,
            output_dir=self.root_dir,
        )
        self.logger.info("Pipeline completed successfully! Summary: %s", out_summary)
        return out_summary

    def run(self, step: str = "all", start_from: str | None = None, resume: bool = False, force: bool = False) -> None:
        """Run requested steps in sequence with continuation support."""
        if start_from:
            start_idx = STEPS.index(start_from)
            steps_to_run = STEPS[start_idx:]
        elif step != "all":
            steps_to_run = [step]
        else:
            steps_to_run = STEPS

        self.logger.info("Running pipeline steps: %s", " -> ".join(steps_to_run))

        for s in steps_to_run:
            if s == "prediction":
                self.run_step_prediction(resume=resume, force=force)
            elif s == "similarity":
                self.run_step_similarity(resume=resume, force=force)
            elif s == "functional":
                self.run_step_functional(resume=resume, force=force)
            elif s == "transcriptomics":
                self.run_step_transcriptomics(resume=resume, force=force)
            elif s == "proteomics":
                self.run_step_proteomics(resume=resume, force=force)
            elif s == "summary":
                self.run_step_summary(resume=resume, force=force)


def main(argv: list[str] | None = None) -> int:
    parser = cli()
    args = parser.parse_args(argv)

    pipeline = AnnotaPipeline(
        config_path=Path(args.config),
        genome_file=Path(args.seq) if args.seq else None,
        protein_file=Path(args.protein) if args.protein else None,
        gff_file=Path(args.gff) if args.gff else None,
    )

    pipeline.run(
        step=args.single_step,
        start_from=args.start_from,
        resume=args.resume,
        force=args.force,
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
