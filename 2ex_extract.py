#!/usr/bin/env python3
"""
Extract individual exon sequences *and* optional full-length transcript
sequences from a GFF3/GTF annotation and a matching genome FASTA.

"""
from __future__ import annotations

import argparse
import logging
import os
from collections import defaultdict
from pathlib import Path
from typing import Iterable, List, Mapping, Sequence

import gffutils
from Bio import SeqIO
from Bio.Seq import Seq

###############################################################################
# Accepted attribute keys (edit here if your GFF/GTF uses other names) ########
###############################################################################

# Gene biotype / gene class
ACCEPTED_GENE_BIOTYPE_KEYS: tuple[str, ...] = (
    "featuretype",   # Sometimes MAKER/EVM promotes the feature type field
    "gene_biotype",
    "biotype",
    "gene_type",
)

# Transcript identifier keys found in exon or transcript features
ACCEPTED_TRANSCRIPT_ID_KEYS: tuple[str, ...] = (
    "transcript_id",
    "orig_transcript_id",
    "Parent"
    "ID",         # last‑ditch fallback
)

# Gene name keys
ACCEPTED_GENE_NAME_KEYS: tuple[str, ...] = (
    "Name",
    "gene_name",
    "gene_id",
    "ID",
)

# Explicit transcript featuretypes we’ll recognise
TRANSCRIPT_TYPES: tuple[str, ...] = (
    "transcript",
    "mRNA",
    "lnc_RNA",
    "ncRNA",
    "tRNA",
    "pseudogenic_transcript",
)

###############################################################################
# Utility helpers #############################################################
###############################################################################

def first_attr(feature: "gffutils.Feature", keys: Sequence[str], default: str = "") -> str:
    """Return the first non‑empty attribute value for *keys*.

    Search order:
    1. ``feature.attributes`` (standard GFF/GTF)
    2. direct property access (Maker/EVM sometimes promotes attrs)
    """
    for key in keys:
        # 1. attribute dict
        if key in feature.attributes and feature.attributes[key]:
            return feature.attributes[key][0]
        # 2. property access
        val = getattr(feature, key, None)
        if val:
            return val[0] if isinstance(val, (list, tuple)) else str(val)
    return default


def database_path(gff_file: str | os.PathLike) -> str:
    root, _ = os.path.splitext(str(gff_file))
    return f"{root}.db"

###############################################################################
# Core routines ###############################################################
###############################################################################

def build_gff_db(gff_file: str | os.PathLike, *, force: bool = False) -> str:
    """Build (if needed) and return the gffutils DB path."""
    db_path = database_path(gff_file)
    if Path(db_path).exists() and not force:
        logging.info("Using existing database %s", db_path)
        return db_path

    logging.info("Creating gffutils database %s …", db_path)
    gffutils.create_db(
        data=str(gff_file),
        dbfn=db_path,
        force=True,
        keep_order=True,
        merge_strategy="merge",
        sort_attribute_values=True,
        verbose=True,
    )
    logging.info("Database built")
    return db_path


def exon_fasta(
    gff_file: str | os.PathLike,
    fasta_file: str | os.PathLike,
    out_fasta: str | os.PathLike,
    *,
    allowed_biotypes: Iterable[str] | None = None,
) -> None:
    """Write *out_fasta* containing per‑exon FASTA records.  Each transcript –
    even alternative isoforms of the same gene – is exported separately."""
    if allowed_biotypes is None:
        allowed_biotypes = {
            "gene",
            "protein_coding",
            "pseudogene",
            "lncRNA",
            "transcribed_pseudogene",
        }

    db = gffutils.FeatureDB(build_gff_db(gff_file), keep_order=True)

    # Genome FASTA ➜ dict
    logging.info("Reading genome FASTA …")
    with open(fasta_file) as fh:
        genome = SeqIO.to_dict(SeqIO.parse(fh, "fasta"))
    genome = {rec.id.split()[0]: rec for rec in genome.values()}

    logging.info("Processing genes …")
    genes = list(db.features_of_type(["gene", "pseudogene"]))

    with open(out_fasta, "w") as out:
        for gene in genes:
            biotype = first_attr(gene, ACCEPTED_GENE_BIOTYPE_KEYS, default="unknown")
            if allowed_biotypes and biotype not in allowed_biotypes and biotype != "unknown":
                continue

            gene_name = (
                first_attr(gene, ACCEPTED_GENE_NAME_KEYS) or gene.id
            ).replace("-", "_")
            chrom = gene.seqid
            strand = gene.strand

            # ------------------------------------------------------------------
            # 1. Prefer explicit transcript children (mRNA etc.)
            # ------------------------------------------------------------------
            transcripts = list(db.children(gene, featuretype=TRANSCRIPT_TYPES, order_by="start"))
            if transcripts:
                logging.debug("%s: %d transcript children", gene_name, len(transcripts))
            else:
                logging.debug("%s: no transcript children; will infer from exons", gene_name)

            if transcripts:
                for t_idx, tx in enumerate(transcripts, start=1):
                    tx_id = first_attr(tx, ACCEPTED_TRANSCRIPT_ID_KEYS, default="")
                    if not tx_id:
                        tx_id = f"tx{t_idx}"
                    tx_id = tx_id.replace("-", "_")

                    exons = list(db.children(tx, featuretype="exon", order_by="start"))
                    if not exons:
                        # Some annotations nest exons directly under gene even when transcript exists
                        exons = list(db.children(gene, featuretype="exon", order_by="start"))

                    write_exons(out, exons, gene_name, tx_id, biotype, chrom, strand, genome)
            else:
                # ------------------------------------------------------------------
                # 2. Fallback: group exons by transcript_id / Parent on the exon itself
                # ------------------------------------------------------------------
                exons_by_tx: Mapping[str, List[gffutils.Feature]] = defaultdict(list)
                for exon in db.children(gene, featuretype="exon", order_by="start"):
                    tx_id = first_attr(exon, ACCEPTED_TRANSCRIPT_ID_KEYS, default="")
                    if not tx_id:
                        tx_id = gene_name  # single‑transcript genes collapse here
                    tx_id = tx_id.replace("-", "_")
                    exons_by_tx[tx_id].append(exon)

                for tx_id, exons in exons_by_tx.items():
                    write_exons(out, exons, gene_name, tx_id, biotype, chrom, strand, genome)

    logging.info("Wrote %s", out_fasta)


def write_exons(
    handle,
    exon_list: List[gffutils.Feature],
    gene_name: str,
    tx_id: str,
    biotype: str,
    chrom: str,
    strand: str,
    genome: Mapping[str, SeqIO.SeqRecord],
) -> None:
    """Helper to write exon FASTA for a single transcript."""
    exon_list.sort(key=lambda x: x.start)
    if strand == "-":
        exon_list = exon_list[::-1]

    for count, exon in enumerate(exon_list, start=1):
        seq_slice = genome[chrom][exon.start - 1 : exon.end].seq
        if strand == "-":
            seq_slice = seq_slice.reverse_complement()
        header = f"{gene_name}-{tx_id}-exon_{count}-{biotype}-{chrom}"
        handle.write(f">{header}\n{seq_slice}\n")


def concatenate_exons(exon_fasta_file: str | os.PathLike, out_fasta: str | os.PathLike) -> None:
    """Concatenate per‑exon FASTA into per‑transcript FASTA."""

    logging.info("Concatenating …")
    with open(exon_fasta_file) as fh:
        records = SeqIO.to_dict(SeqIO.parse(fh, "fasta"))

    transcripts: Mapping[str, List[tuple[str, Seq]]] = defaultdict(list)
    for key, rec in records.items():
        gene_name, tx_id, exon_tag, biotype, chrom = key.split("-")
        header_base = f"{gene_name}-{tx_id}-{biotype}-{chrom}"
        transcripts[header_base].append((exon_tag, rec.seq))

    with open(out_fasta, "w") as out:
        for header, seqs in transcripts.items():
            seqs.sort(key=lambda kv: int(kv[0].split("_")[1]))
            concat_seq = Seq("".join(str(seq) for _, seq in seqs))
            out.write(f">{header}\n{concat_seq}\n")

    logging.info("Wrote %s", out_fasta)

###############################################################################
# Command‑line interface ######################################################
###############################################################################

def main() -> None:
    parser = argparse.ArgumentParser(description="Extract exon/transcript FASTA from GFF3/GTF")
    sub = parser.add_subparsers(dest="cmd", required=True)

    p_build = sub.add_parser("build", help="create gffutils database")
    p_build.add_argument("gff_file")
    p_build.add_argument("--force", action="store_true")

    p_extract = sub.add_parser("extract", help="extract exon FASTA")
    p_extract.add_argument("gff_file")
    p_extract.add_argument("genome_fasta")
    p_extract.add_argument("-e", "--exons", default="exons.fasta")
    p_extract.add_argument("-c", "--concatenated", default="concatenated_exons.fasta")
    p_extract.add_argument("--skip-concat", action="store_true")
    p_extract.add_argument("-v", "--verbose", action="count", default=0)

    args = parser.parse_args()

    # logging setup -----------------------------------------------------------
    level = logging.ERROR  # default
    if args.verbose == 1:
        level = logging.INFO
    elif args.verbose >= 2:
        level = logging.DEBUG
    logging.basicConfig(level=level, format="%(levelname)s: %(message)s")

    # dispatch ---------------------------------------------------------------
    if args.cmd == "build":
        build_gff_db(args.gff_file, force=args.force)
    elif args.cmd == "extract":
        exon_fasta(args.gff_file, args.genome_fasta, args.exons)
        if not args.skip_concat:
            concatenate_exons(args.exons, args.concatenated)


if __name__ == "__main__":
    main()
