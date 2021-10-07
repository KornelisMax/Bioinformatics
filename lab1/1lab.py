from Bio import SeqIO
import os
from scipy.spatial import distance_matrix
import pandas as pd

IUPAC_PROTEIN_ALPHABET = list("ACDEFGHIKLMNPQRSTVWY")
TRANSLATION_TABLE = 11 #https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi#SG11
MAX_ORF_LEN = 100
ROOT = 'data'

def get_codon_dicodon_frequencies():   
    codon_frequencies = {}
    dicodon_frequencies = {}

    for file_name in os.listdir(ROOT):
        record = import_fasta_data(file_name)
        orf_list = find_orfs(record.seq, TRANSLATION_TABLE, MAX_ORF_LEN)
        codon_frequencies[record.name] = get_codon_frequency(orf_list).values()
        dicodon_frequencies[record.name] = get_dicodon_frequency(orf_list).values()
    return codon_frequencies, dicodon_frequencies

def import_fasta_data(file_name):
    filepath = ROOT +'/' + file_name
    
    return SeqIO.read(filepath, "fasta")

def find_orfs(sequence, translation_table, min_protein_length):
    orf_list = []
    protein_sequence_len = len(sequence)
    for strand, nuc in [(+1, sequence), (-1, sequence.reverse_complement())]:
        for frame in range(3):
            protein_sequence = str(nuc[frame:].translate(translation_table))
            protein_sequence_len = len(protein_sequence)
            orf_start = 0
            orf_end = 0
            while orf_start < protein_sequence_len:
                orf_end = protein_sequence.find("*", orf_start)
                if orf_end == -1:
                    orf_end = protein_sequence_len
                if orf_end - orf_start >= min_protein_length:
                    if strand == 1:
                        start = frame + orf_start * 3
                        end = min(protein_sequence_len, frame + orf_end * 3 + 3)
                    else:
                        start = protein_sequence_len - frame - orf_end * 3 - 3
                        end = protein_sequence_len - frame - orf_start * 3
                    orf_list.append((start, end, strand, protein_sequence[orf_start:orf_end]))
                orf_start = orf_end + 1
    orf_list.sort()
    return orf_list

def get_codon_frequency(orf_list):
    codon_hash = dict.fromkeys(IUPAC_PROTEIN_ALPHABET, 0)
    total_number_of_codons = 0
    for _start, _end, _strand, protein_sequence in orf_list:
        total_number_of_codons = total_number_of_codons + len(protein_sequence)
        for codon in list(protein_sequence):
            codon_hash[codon] += 1
    for key in codon_hash:
        codon_hash[key] = codon_hash[key] / total_number_of_codons
    return codon_hash

def get_dicodon_frequency(orf_list):
    codon_hash = dict.fromkeys([a + b for a in IUPAC_PROTEIN_ALPHABET for b in IUPAC_PROTEIN_ALPHABET], 0)
    total_number_of_dicodons = 0
    for _start, _end, _strand, protein_sequence in orf_list:
        if (len(protein_sequence) % 2 != 0):
            protein_sequence = protein_sequence[:-1]
        total_number_of_dicodons = total_number_of_dicodons + (len(protein_sequence) / 2)
        dicodon_list = [protein_sequence[i:i+2] for i in range(0, len(protein_sequence), 2)]
        for dicodon in dicodon_list:
            codon_hash[dicodon] += 1
    for key in codon_hash:
        codon_hash[key] = codon_hash[key] / total_number_of_dicodons
    return codon_hash

def print_distance_matrix(frequencies):
    names = []
    data = []
    for key in frequencies:
        names.append(key)
        data.append(list(frequencies[key]))
    df = pd.DataFrame(data, index=names)
    print(pd.DataFrame(distance_matrix(df.values, df.values), index=df.index, columns=df.index))

codon_frequencies, dicodon_frequencies = get_codon_dicodon_frequencies()

print('\nCodon distance_matrix:')
print_distance_matrix(codon_frequencies)

print('\nDicodon distance_matrix:')
print_distance_matrix(dicodon_frequencies)