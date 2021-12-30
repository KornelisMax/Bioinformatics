from Bio import SeqIO
import matplotlib.pyplot as plt; plt.rcdefaults()
import numpy as np
import matplotlib as mpl
import scipy.signal as fp
import matplotlib.pyplot as plt
import numpy as np
from scipy.signal import argrelextrema
from matplotlib.pyplot import *
from Bio.Blast import NCBIWWW
import xml.etree.ElementTree as elTree
plt.rcParams.update({'font.size': 8})

def find_encoding():
    sanger = 0
    solexa = 0
    illumina1_3 = 0
    illumina1_5 = 0
    illumina1_8 = 0
    lowest_quality = 100
    highest_quality = 0
    with open("reads_for_analysis.fastq", "rU") as handle:
        for record in SeqIO.parse(handle, "fastq"):
            quality_list = record.letter_annotations['phred_quality']
            lowest = min(quality_list)
            highest = max(quality_list)
            if lowest + 33 >= ord("!"):
                if (lowest < lowest_quality):
                    lowest_quality = lowest
                if (highest_quality < highest):
                    highest_quality = highest
                if highest + 33 == ord("J"):

                    illumina1_8 = illumina1_8 + 1
                if highest + 33 < ord("J"):
                    sanger = sanger + 1

    lowest_quality = lowest_quality + 33
    highest_quality = highest_quality + 33

    if lowest_quality > ord("B") and highest_quality < ord("i"):
        illumina1_5 = illumina1_5 + 1
    if lowest_quality > ord("@") and highest_quality < ord("h"):
        illumina1_3 = illumina1_3 + 1
    if lowest_quality > ord(";") and highest_quality < ord("h"):
        solexa = solexa + 1
    if sanger > 0:
        print("Sanger Phred+33")
    else:
        print("Illumina 1.8+ Phred+33")
    if illumina1_5 > 0:
        print("illumina1_5")
    if illumina1_3 > 0:
        print("illumina1_3")
    if solexa > 0:
        print("solexa")
            
def generate_percentage_dictionary():
    dict = {}
    for i in range(0, 101, 1):
        dict[i] = 0
    return dict

def get_c_g_content():
    c_g_content = generate_percentage_dictionary()
    all_percentages = []

    with open("reads_for_analysis.fastq", "rU") as handle:
        for record in SeqIO.parse(handle, "fastq"):
            seq = record.seq
            c_g_sum = seq.count("C") + seq.count("G")
            total_sum = len(seq)
            percentage = int((c_g_sum / total_sum) * 100)
            all_percentages.append(percentage)
            c_g_content[percentage] += 1
    return all_percentages, c_g_content

def find_local_peaks(percentage_dict, x1, x2):
    local_peaks = [0] * (x2 + 1)
    values = list(percentage_dict.values())
    peaks, _ = fp.find_peaks(values)
    val = [x for x in peaks if x <= x2 and x >= x1]
    for p in val:
        local_peaks[p] = values[p]

    peak = max(local_peaks)
    max_index = local_peaks.index(peak)
    return max_index

def print_chart(percentage_dict):
    objects = percentage_dict.keys()
    y_pos = np.arange(len(objects))
    performance = percentage_dict.values()
    plt.plot(y_pos, performance, linewidth=1)
    plt.xlabel('C/G content in sequence (%)')
    plt.ylabel('Number of reads')
    plt.title('Number of reads with certain C/G content in sequence')
    plt.show()

def get_peak_seqs(peak):
    sequences_dict = {}
    nr = 5
    with open("reads_for_analysis.fastq", "rU") as handle:
        for record in SeqIO.parse(handle, "fastq"):
            seq = record.seq
            c_g_sum = seq.count("C") + seq.count("G")
            total_sum = len(seq)
            percentage = int((c_g_sum / total_sum) * 100)
            if (percentage is peak):
                sequences_dict[record.id] = seq
                # peaks_list.remove(percentage)
                nr = nr - 1
                if nr is 0:
                    return sequences_dict

def do_blast_search(seq):
    result_handle = NCBIWWW.qblast("blastn", "nt", seq, hitlist_size=1, entrez_query="bacteria[organism]")
    result_xml = result_handle.read()
    print(result_xml)
    result_title = ""
    try:
        result_title = elTree.fromstring(result_xml)\
            .find("BlastOutput_iterations")\
            .find("Iteration")\
            .find("Iteration_hits")\
            .find("Hit")\
            .find("Hit_def")\
            .text
    except AttributeError:
        result_title = "EMPTY"
    return result_title

def analyze_seqs(peaks):
    myfile = open('analysis.txt', 'w')
    myfile.write("Read ID, organism\n")
    peaks_sequences1 = get_peak_seqs(find_local_peaks(peaks, 0, 40))
    peaks_sequences2 = get_peak_seqs(find_local_peaks(peaks, 40, 70))
    peaks_sequences3 = get_peak_seqs(find_local_peaks(peaks, 70, 100))
    for k, v in peaks_sequences1.items():
        res = do_blast_search(v)
        myfile.write("%s, %s\n" % (k, res))
        print("%s, %s" % (k, v))
    for k, v in peaks_sequences2.items():
        res = do_blast_search(v)
        myfile.write("%s, %s\n" % (k, res))
        print("%s, %s" % (k, v))
    for k, v in peaks_sequences3.items():
        res = do_blast_search(v)
        myfile.write("%s, %s\n" % (k, res))
        print("%s, %s" % (k, v))
    myfile.close()

find_encoding() # outputs encoding type of given fastq file; 1 task 
all_percentages, c_g_content = get_c_g_content() # finds C/G content according to seqs; 
print("Number of all reads: %d\n" % len(all_percentages)) # prints all seqs with appropriate C/G content %
print_chart(c_g_content) # visualizes seqs with C/G percentage
# analyze_seqs(c_g_content) # Analyzes sequences: finds main peaks, 5 sequences from each peak, does blast search, writes output to txt file