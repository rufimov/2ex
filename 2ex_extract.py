#!usr/local/bin/python3
import argparse

import gffutils
from Bio import SeqIO


def build_gff_db():
    print("Creating gff database...")
    gffutils.create_db(
        data=gff,
        dbfn=gff[:-4],
        force=True,
        keep_order=True,
        merge_strategy="merge",
        sort_attribute_values=True,
        verbose=True,
    )
    print("Done")


def extract_exons():
    print("Creating feature database...")
    feature_db = gffutils.FeatureDB(gff[:-4], keep_order=True)
    genes = list(feature_db.features_of_type("gene"))
    pseudogenes = list(feature_db.features_of_type("pseudogene"))
    genes.extend(pseudogenes)
    print("Done")
    print("Reading genome fasta...")
    with open(fasta_file) as fasta:
        fasta_parsed = SeqIO.to_dict(SeqIO.parse(fasta, "fasta"))
    print("Done")
    print("Writing fasta with exons...")
    correct_fasta_parsed = {}
    for key in fasta_parsed:
        correct_fasta_parsed[key.split()[0]] = fasta_parsed[key]
    with open(result_fasta, "w") as all_exons:
        for gene in genes:
            print(gene)
            type_of_feature = gene.attributes["gene_biotype"][0]
            if (
                type_of_feature == "protein_coding"
                or type_of_feature == "pseudogene"
                or type_of_feature == "lncRNA"
                or type_of_feature == "transcribed_pseudogene"
            ):
                strand = gene.strand
                name = gene.attributes["Name"][0].replace("-", "_")
                chromosome = gene.seqid
                exons = []
                for i in feature_db.children(
                    gene, featuretype="exon", order_by="start"
                ):
                    print(exons)
                    exons.append(i)
                try:
                    exons.sort(key=lambda x: x.start)
                    exons.sort(key=lambda x: x.attributes["transcript_id"][0])
                except KeyError:
                    try:
                        exons.sort(key=lambda x: x.start)
                        exons.sort(key=lambda x: x.attributes["orig_transcript_id"][0])
                    except KeyError:
                        exons.sort(key=lambda x: x.start)
                if strand == "-":
                    exons = exons[::-1]
                    count = 1
                    transcripts = set()
                    for exon in exons:
                        try:
                            transcript = exon.attributes["transcript_id"][0].replace(
                                "-", "_"
                            )
                        except KeyError:
                            try:
                                transcript = exon.attributes["orig_transcript_id"][
                                    0
                                ].replace("-", "_")
                            except KeyError:
                                transcript = name
                        if transcript not in transcripts:
                            count = 1
                            transcripts.add(transcript)
                        start = exon.start
                        end = exon.end
                        sequence = str(
                            correct_fasta_parsed[chromosome][
                                start - 1 : end
                            ].seq.reverse_complement()
                        )
                        all_exons.write(
                            f">{name}-{transcript}-exon_{str(count)}-{type_of_feature}-{chromosome}\n"
                            f"{sequence}\n"
                        )
                        count += 1
                elif strand == "+":
                    count = 1
                    transcripts = set()
                    for exon in exons:
                        try:
                            transcript = exon.attributes["transcript_id"][0].replace(
                                "-", "_"
                            )
                        except KeyError:
                            try:
                                transcript = exon.attributes["orig_transcript_id"][
                                    0
                                ].replace("-", "_")
                            except KeyError:
                                transcript = name
                        if transcript not in transcripts:
                            count = 1
                            transcripts.add(transcript)
                        start = exon.start
                        end = exon.end
                        sequence = str(
                            correct_fasta_parsed[chromosome][start - 1 : end].seq
                        )
                        all_exons.write(
                            f">{name}-{transcript}-exon_{str(count)}-{type_of_feature}-{chromosome}\n"
                            f"{sequence}\n"
                        )
                        count += 1


def concatenate():
    print("Done")
    print("Concatenating exons...")
    with open(result_fasta) as fasta_to_concatenate, open(
        concatenated_fasta, "w"
    ) as concat_fasta:
        fasta_parsed = SeqIO.to_dict(SeqIO.parse(fasta_to_concatenate, "fasta"))
        current_transcript = ""
        count = 1
        list_of_keys = list(fasta_parsed.keys())
        list_of_keys.sort(key=lambda x: int(x.split("-")[2].split("_")[1]))
        list_of_keys.sort(key=lambda x: x.split("-")[1])
        list_of_keys.sort(key=lambda x: x.split("-")[0])
        list_of_keys.sort(key=lambda x: x.split("-")[-1])
        for key in list_of_keys:
            transcript = key.split("-")[1]
            locus = key.split("-")[0]
            type_of_feature = key.split("-")[3]
            chromosome = key.split("-")[4]
            if count == 1:
                if transcript != current_transcript:
                    concat_fasta.write(
                        ">"
                        + locus
                        + "-"
                        + transcript
                        + "-"
                        + type_of_feature
                        + "-"
                        + chromosome
                        + "\n"
                        + str(fasta_parsed[key].seq)
                    )
                else:
                    concat_fasta.write(str(fasta_parsed[key].seq))
            else:
                if transcript != current_transcript:
                    concat_fasta.write(
                        "\n>"
                        + locus
                        + "-"
                        + transcript
                        + "-"
                        + type_of_feature
                        + "-"
                        + chromosome
                        + "\n"
                        + str(fasta_parsed[key].seq)
                    )
                else:
                    concat_fasta.write(str(fasta_parsed[key].seq))
            current_transcript = transcript
            count += 1
    print("Done")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "command",
        help="command to be executed: 'build' -- to build database with genome annotations, "
        "'extract' -- to extract exons and create a file with predicted transcripts",
        choices=["build", "extract"],
    )
    parser.add_argument("gff_file", help="file with genome annotations in gff format")
    parser.add_argument(
        "genome_fasta",
        help="file with genome sequences in fasta format, "
        "be sure that sequences in fasta and gff file are named identically",
    )
    parser.add_argument(
        "-e",
        "--exons",
        default="exons.fasta",
        help="name of the file to write extracted exons",
    )
    parser.add_argument(
        "-c",
        "--concatenated",
        default="concatenated_exons.fasta",
        help="name of the file " "with concatenated exons",
    )
    args = parser.parse_args()
    task = args.command
    gff = args.gff_file
    fasta_file = args.genome_fasta
    result_fasta = args.exons
    concatenated_fasta = args.concatenated
    if task == "build":
        build_gff_db()
    else:
        extract_exons()
        concatenate()
