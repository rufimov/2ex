#!/usr/bin/env python3
"""blast_pipeline.py

A tidied-up, self-contained replacement for the *BLAST* stage of the exon
pipeline.  Compared with the original script, this version

* drops the deprecated *Bio.Alphabet* import,
* uses **Pathlib**, **logging**, and type hints for clarity,
* collapses repeated `list.sort()` chains into single `sorted(…, key=…)` calls
  with tuple keys so the intent is obvious and stable,
* splits work into small, testable helper functions, and
* exposes a single CLI (`python blast_pipeline.py run …`) that reproduces the
  old behaviour but is easier to tweak.

Dependencies
~~~~~~~~~~~~
* Biopython ≥ 1.75
* BLAST+ in your ``$PATH`` (``makeblastdb`` / ``blastn``)
"""
# from __future__ import annotations

import argparse
import logging
import os
import shutil
import subprocess
import tempfile
from collections import defaultdict
from pathlib import Path
from typing import Iterable, List, Sequence

from Bio import SeqIO

###############################################################################
# BLAST wrappers ##############################################################
###############################################################################

def make_blast_db(fasta: Path, db_prefix: Path) -> None:
    """Create a nucleotide BLAST database from *fasta* with basename *db_prefix*."""
    logging.info("Building BLAST DB %s …", db_prefix)
    subprocess.run(
        [
            "makeblastdb",
            "-in",
            str(fasta),
            "-dbtype",
            "nucl",
            "-parse_seqids",
            "-out",
            str(db_prefix),
        ],
        check=True,
    )


def run_blastn(
    query: Path,
    db_prefix: Path,
    out_file: Path,
    *,
    task: str = "blastn",
    outfmt: str = "6 qaccver saccver pident qcovhsp evalue bitscore",
    threads: int = 4,
) -> None:
    """Run *blastn* against the given database."""
    logging.info("BLASTing %s against %s …", query, db_prefix)
    subprocess.run(
        [
            "blastn",
            "-task",
            task,
            "-query",
            str(query),
            "-db",
            str(db_prefix),
            "-outfmt",
            outfmt,
            "-out",
            str(out_file),
            "-num_threads",
            str(threads),
        ],
        check=True,
    )

###############################################################################
# FASTA helpers ###############################################################
###############################################################################

def write_corrected_names(src_fasta: Path, dst_fasta: Path) -> None:
    """Write *dst_fasta* where each header becomes <transcript>-<biotype>."""
    with src_fasta.open() as handle:
        records = SeqIO.to_dict(SeqIO.parse(handle, "fasta"))

    with dst_fasta.open("w") as out:
        for key, rec in records.items():
            parts = key.split("-")
            transcript = parts[1]
            biotype = parts[3]
            out.write(f">{transcript}-{biotype}\n{rec.seq}\n")

###############################################################################
# Sorting helpers #############################################################
###############################################################################

def best_hits(in_table: Path, out_table: Path) -> List[str]:
    """
    Keep just one (best) BLAST hit for every probe-ID that follows the
    first dash in the query column.

    Parameters
    ----------
    in_table : Path
        6-column BLAST tabular file
        (qaccver saccver pident qcovhsp evalue bitscore)
    out_table : Path
        File to write the deduplicated table to.

    Returns
    -------
    list[str]
        The surviving lines, in the same order they are written to `out_table`.
    """
    with in_table.open() as fh:
        lines = fh.readlines()

    # ------------------------------------------------------------------
    # 1. sort so the "longest-coverage, highest-identity, lowest e-value,
    #    high bitscore" hit is *first* inside each probe group
    # ------------------------------------------------------------------
    lines.sort(
        key=lambda ln: (
            ln.split()[0].split("-")[1],   # group key  (probe token after '-')
            -float(ln.split()[3]),         # 1️⃣ larger qcovhsp  (% coverage)
            -float(ln.split()[2]),         # 2️⃣ larger pident   (% identity)
            float(ln.split()[4]),          # 3️⃣ smaller e-value
            -float(ln.split()[5]),         # 4️⃣ larger bitscore
        )
    )

    # ------------------------------------------------------------------
    # 2. keep only the first line for every probe token
    # ------------------------------------------------------------------
    seen: set[str] = set()
    best_hits: list[str] = []

    for ln in lines:
        probe_token = ln.split()[0].split("-")[1]
        if probe_token not in seen:
            best_hits.append(ln)
            seen.add(probe_token)

    # ------------------------------------------------------------------
    # 3. write the deduplicated table (so downstream steps use it)
    # ------------------------------------------------------------------
    with out_table.open("w") as out:
        out.writelines(best_hits)

    logging.info("Best-hit table written to %s (%d rows)", out_table, len(best_hits))
    return best_hits

###############################################################################
# Pipeline steps ##############################################################
###############################################################################

def split_best_hits_to_exons(
    best_hits_lines: Sequence[str],
    separated_exons: Path,
    output_fasta: Path,
) -> None:
    """Write one FASTA per exon that belongs to a best hit locus."""
    with separated_exons.open() as handle:
        all_exons = SeqIO.to_dict(SeqIO.parse(handle, "fasta"))

    with output_fasta.open("w") as out:
        for hit in best_hits_lines:
            locus = hit.split()[1].split("-")[0]
            probe = hit.split()[0]
            for key, rec in all_exons.items():
                if locus+'-' in key:
                    out.write(f">{probe}_{rec.id}\n{rec.seq}\n")


def filter_and_sort_exon_hits(lines: Sequence[str]) -> List[str]:
    """Apply the exon-level filtering + dedup sorting from the original script."""
    # Accept only lines where query and subject share probe name

    lines = [l for l in lines if l.split()[0].split("_")[0].split('-')[1] == l.split()[1].split('-')[1]]

    # Composite sort mirroring original priority
    def composite(line: str):
        cols = line.split()
        return (
            cols[1].split("-")[1],                      # probe base
            # cols[0].split("-")[2],                      # gene/tx id
            int(cols[0].split("-")[3].split("_")[1]), # exon number (ascending)
            -float(cols[3]),                            # qcov desc
            -float(cols[2]),                            # pident desc
            float(cols[4]),                             # evalue asc
            -float(cols[5]),                            # bitscore desc
        )

    sorted_lines = sorted(lines, key=composite)

    deduped: list[str] = []
    seen: set[str] = set()
    for line in sorted_lines:
        query = line.split()[0]
        if query not in seen:
            deduped.append(line)
            seen.add(query)
    return deduped


def extract_sequences_from_hits(
    hits: Sequence[str],
    probes_fa: Path,
    exons_fa: Path,
    result_loci: Path,
    result_exons: Path,
) -> None:
    """From BLAST exon hits produce two FASTAs mimicking the original outputs."""
    probes = SeqIO.to_dict(SeqIO.parse(probes_fa.open(), "fasta"))
    exons = SeqIO.to_dict(SeqIO.parse(exons_fa.open(), "fasta"))

    with result_loci.open("w") as out1, result_exons.open("w") as out2:
        for hit in hits:
            cols = hit.split()
            locus = cols[1]
            locus_to_write = cols[0].split('-')[2]+'-'+cols[1].split('-')[1]
            exon_query = cols[0]
            exon_num = exon_query.split("-")[3].split("_")[1]
            s_start, s_end = map(int, cols[6:8])  # subject positions in query–>probe blast
            q_start, q_end = map(int, cols[8:10])  # query positions in exon db

            # slice probe sequence (subject)
            probe_seq = probes[locus].seq
            if s_start > s_end:
                subseq = probe_seq[s_end - 1 : s_start].reverse_complement()
            else:
                subseq = probe_seq[s_start - 1 : s_end]
            out1.write(f">{locus_to_write}_exon_{exon_num}\n{subseq}\n")

            # slice exon sequence (query)
            exon_seq = exons[exon_query].seq
            if q_start > q_end:
                subseq_exon = exon_seq[q_end - 1 : q_start].reverse_complement()
            else:
                subseq_exon = exon_seq[q_start - 1 : q_end]
            out2.write(f">{locus_to_write}_exon_{exon_num}\n{subseq_exon}\n")


def concatenate_optimized(result_exons: Path) -> Path:
    """Concatenate exon fragments back into full transcripts (optimized probe set)."""
    concat_out = result_exons.with_name(result_exons.stem + "_concatenated" + result_exons.suffix)

    records = SeqIO.to_dict(SeqIO.parse(result_exons.open(), "fasta"))
    grouped: defaultdict[str, list[str]] = defaultdict(list)
    for key in records:
        base, locus = key.split("-")
        grouped[f"{base}-{locus.split('_')[0]}"] .append(key)

    with concat_out.open("w") as out:
        for header, keys in grouped.items():
            keys.sort(key=lambda k: int(k.split("-")[1].split("_")[2]))
            concat_seq = "".join(str(records[k].seq) for k in keys)
            out.write(f">{header}\n{concat_seq}\n")

    logging.info("Wrote %s", concat_out)
    return concat_out

###############################################################################
# Orchestration ###############################################################
###############################################################################

def run_pipeline(
    concatenated_exons: Path,
    separated_exons: Path,
    probes: Path,
    *,
    blast_task: str = "blastn",
    threads: int = 4,
) -> None:
    tmpdir = Path(tempfile.mkdtemp(prefix="blastpipe_"))
    try:
        # 1. Correct names in concatenated exons
        corrected_fa = tmpdir / "corrected_names.fasta"
        write_corrected_names(concatenated_exons, corrected_fa)

        # 2. Build DB + BLAST probes vs concatenated exons
        db_prefix = tmpdir / "concat_exons_db"
        make_blast_db(corrected_fa, db_prefix)
        blast_out1 = tmpdir / "probes_vs_concat.txt"
        blast_out1_dedup = tmpdir / "probes_vs_concat_dedup.txt"
        run_blastn(probes, db_prefix, blast_out1, task=blast_task, threads=threads)

        # 3. Choose best hit per probe
        best1 = best_hits(blast_out1, blast_out1_dedup)

        # 4. Write exon FASTA for best loci
        best_exons_fa = Path("best_hits_exons.fasta")
        split_best_hits_to_exons(best1, separated_exons, best_exons_fa)

        # 5. Build DB of probes and BLAST exons back
        probes_db_prefix = tmpdir / "probes_db"
        make_blast_db(probes, probes_db_prefix)  # BLAST+ will add .n* files next to probes
        blast_out2 = tmpdir / "exons_vs_probes.txt"
        run_blastn(best_exons_fa, probes_db_prefix, blast_out2,
                   task=blast_task,
                   threads=threads,
                   outfmt="6 qaccver saccver pident qcovhsp evalue bitscore sstart send qstart qend")

        cleaned = filter_and_sort_exon_hits(blast_out2.read_text().strip().split("\n"))
        # Write cleaned hits to a file for debugging
        cleaned_file = tmpdir / "cleaned_exon_hits.txt"
        cleaned_file.write_text("\n".join(cleaned) + "\n")
        # 6. Slice sequences and write final FASTAs
        result_loci = Path("probes_separated_to_exons.fasta")
        result_exons = Path("optimized_probes.fasta")
        extract_sequences_from_hits(cleaned, probes, best_exons_fa, result_loci, result_exons)

        # 7. Concatenate optimized probes
        concatenate_optimized(result_exons)

    finally:
        shutil.rmtree(tmpdir)

###############################################################################
# CLI #########################################################################
###############################################################################

def main() -> None:
    p = argparse.ArgumentParser(description="Run BLAST-based probe optimisation pipeline")
    p.add_argument("concatenated_exons", type=Path, help="FASTA of concatenated exons")
    p.add_argument("separated_exons", type=Path, help="FASTA of per-exon sequences")
    p.add_argument("probes", type=Path, help="FASTA of probe sequences")
    p.add_argument("-t", "--blast-task", choices=["blastn", "megablast", "dc-megablast"], default="blastn")
    p.add_argument("-j", "--threads", type=int, default=4, help="Number of BLAST threads")
    p.add_argument("-v", "--verbose", action="count", default=0, help="-v / -vv for INFO/DEBUG")

    args = p.parse_args()

    logging.basicConfig(
        level=logging.ERROR if args.verbose == 0 else logging.INFO if args.verbose == 1 else logging.DEBUG,
        format="%(levelname)s: %(message)s",
    )

    run_pipeline(
        args.concatenated_exons,
        args.separated_exons,
        args.probes,
        blast_task=args.blast_task,
        threads=args.threads,
    )


if __name__ == "__main__":
    main()
