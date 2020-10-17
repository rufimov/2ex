#!usr/local/bin/python3
import argparse
import glob
import os

from Bio import SeqIO
from Bio.Alphabet import generic_dna
from Bio.Blast.Applications import NcbimakeblastdbCommandline, NcbiblastnCommandline


def blast():
    with open(concat_exons) as concatenated:
        concatenated_exons = SeqIO.to_dict(SeqIO.parse(concatenated, 'fasta', generic_dna))
    with open(f"{concat_exons.split('.')[0]}_names_corrected.fas", 'w') as corrected:
        for key in concatenated_exons.keys():
            corrected.write(f">{key.split('-')[1]}-{key.split('-')[3]}\n{str(concatenated_exons[key].seq)}\n")
    print('Building database for %s...' % concat_exons)
    NcbimakeblastdbCommandline(dbtype='nucl', input_file=f"{concat_exons.split('.')[0]}_names_corrected.fas",
                               out=concat_exons, parse_seqids=True)()
    print('Done')
    print(f'Blasting {probes} against {concat_exons}')
    NcbiblastnCommandline(task=blast_task, query=probes, db=concat_exons,
                          out=f'{probes}_against_{concat_exons}.txt',
                          outfmt="6 qaccver saccver pident qcovhsp evalue bitscore",
                          num_threads=4)()
    print('Done')


def best_hit_search(file_with_hittable, list_with_best_hits):
    print('Choosing the best hit...')
    with open(file_with_hittable) as blast_results:
        hits = blast_results.readlines()
    hits.sort(key=lambda x: float(x.split()[5]), reverse=True)
    hits.sort(key=lambda x: float(x.split()[4]))
    hits.sort(key=lambda x: float(x.split()[2]), reverse=True)
    hits.sort(key=lambda x: float(x.split()[3]), reverse=True)
    hits.sort(key=lambda x: x.split()[0].split('-')[1])
    uniqe_hits = set()
    for hit in hits:
        if hit.split()[0].split('-')[1] not in uniqe_hits:
            list_with_best_hits.append(hit)
            uniqe_hits.add(hit.split()[0].split('-')[1])
    print('Done')


def split_to_exons():
    print('Splitting best hits to exons...')
    with open(separat_exons) as all_exons:
        all_exons_parsed = SeqIO.to_dict(SeqIO.parse(all_exons, 'fasta', generic_dna))
    with open(best_separate_exons, 'w') as best_exons:
        for besthit in best_hits:
            locus = besthit.split()[1].split('-')[0]
            probe = besthit.split()[0]
            exons = [val for key, val in all_exons_parsed.items() if locus in key]
            for exon in exons:
                name = str(exon.id)
                sequence = str(exon.seq)
                best_exons.write(f'>{probe}_{name}\n{sequence}\n')
    NcbimakeblastdbCommandline(dbtype='nucl', input_file=probes,
                               out=probes, parse_seqids=True)()
    NcbiblastnCommandline(task=blast_task, query=best_separate_exons, db=probes,
                          out=f'{best_separate_exons}_against_{probes}.txt', num_threads=4,
                          outfmt='6 qaccver saccver pident qcovhsp evalue bitscore sstart send qstart qend')()
    with open(f'{best_separate_exons}_against_{probes}.txt') as new_blast_results:
        hits = new_blast_results.readlines()
    cleaned_hits = []
    for hit in hits:
        if hit.split()[0].split('_')[0] == hit.split()[1]:
            cleaned_hits.append(hit)
    cleaned_hits.sort(key=lambda x: float(x.split()[5]), reverse=True)
    cleaned_hits.sort(key=lambda x: float(x.split()[4]))
    cleaned_hits.sort(key=lambda x: float(x.split()[2]), reverse=True)
    cleaned_hits.sort(key=lambda x: float(x.split()[3]), reverse=True)
    cleaned_hits.sort(key=lambda x: int(x.split()[0].split('-')[3].split('_')[1]))
    cleaned_hits.sort(key=lambda x: x.split()[0].split('-')[2])
    hits_exons = set()
    cleaned_dedup_hits = []
    for cleaned_hit in cleaned_hits:
        if cleaned_hit.split()[0] not in hits_exons:
            cleaned_dedup_hits.append(cleaned_hit)
            hits_exons.add(cleaned_hit.split()[0])
    cleaned_dedup_hits.sort(key=lambda x: int(x.split()[0].split('-')[3].split('_')[1]))
    cleaned_dedup_hits.sort(key=lambda x: x.split()[1].split('-')[1])
    with open(f'{best_separate_exons}_against_{probes}.txt', 'w') as new_blast_results:
        for cleaned_hit in cleaned_dedup_hits:
            new_blast_results.write(cleaned_hit)
    with open(probes) as probes_to_parse:
        probes_as_dict = SeqIO.to_dict(SeqIO.parse(probes_to_parse, 'fasta', generic_dna))
    with open(best_separate_exons) as best_exons:
        best_exons_as_dict = SeqIO.to_dict(SeqIO.parse(best_exons, 'fasta', generic_dna))
    with open(result_file, 'w') as resultfile, open(result_file2, 'w') as resultfile2:
        for cleaned_dedup_hit in cleaned_dedup_hits:
            name_of_locus = cleaned_dedup_hit.split()[1]
            name_of_exon = cleaned_dedup_hit.split()[0]
            num_exon = cleaned_dedup_hit.split()[0].split('-')[3].split('_')[1]
            if int(cleaned_dedup_hit.split()[6]) > int(cleaned_dedup_hit.split()[7]):
                start = int(cleaned_dedup_hit.split()[7])
                end = int(cleaned_dedup_hit.split()[6])
                sequence = str(probes_as_dict[name_of_locus][start - 1:end].seq.reverse_complement())
            else:
                start = int(cleaned_dedup_hit.split()[6])
                end = int(cleaned_dedup_hit.split()[7])
                sequence = str(probes_as_dict[name_of_locus][start - 1:end].seq)

            resultfile.write(f'>{name_of_locus}_exon_{num_exon}\n{sequence}\n')
            if int(cleaned_dedup_hit.split()[8]) > int(cleaned_dedup_hit.split()[9]):
                start_opt = int(cleaned_dedup_hit.split()[9])
                end_opt = int(cleaned_dedup_hit.split()[8])
                sequence_opt = str(best_exons_as_dict[name_of_exon][start_opt - 1:end_opt].seq.reverse_complement())
            else:
                start_opt = int(cleaned_dedup_hit.split()[8])
                end_opt = int(cleaned_dedup_hit.split()[9])
                sequence_opt = str(best_exons_as_dict[name_of_exon][start_opt - 1:end_opt].seq)
            resultfile2.write(f'>{name_of_locus}_exon_{num_exon}\n{sequence_opt}\n')
    print('Done')


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("concatenated_exons", help="")
    parser.add_argument("separated_exons", help="")
    parser.add_argument("probes", help="")
    parser.add_argument("-be", "--best_separated", help="", default="best_hits_as_exons.fasta")
    parser.add_argument("-pe", "--probes_separated", help="", default="probes_separated_to_exons.fasta")
    parser.add_argument("-op", "--optimized_probes", help="", default="optimized_probes.fasta")
    parser.add_argument("-t", "--blast_task", help="", default="blastn", choices=["megablast", "dc-megablast", "blastn"])
    args = parser.parse_args()
    concat_exons = args.concatenated_exons  # concatenated exons
    separat_exons = args.separated_exons  # separated exons
    probes = args.probes  # probe file
    best_separate_exons = args.best_separated
    result_file = args.probes_separated  # final output file
    result_file2 = args.optimized_probes
    blast_task = args.blast_task  # blast task ('megablast', 'dc-megablast', 'blastn')
    blast()
    best_hits = []
    best_hit_search(f'{probes}_against_{concat_exons}.txt', best_hits)
    split_to_exons()
    for file in glob.glob('*.n*'):
        os.remove(file)
