# This code calculates the pobability significance of the the presence of
# "TTAGGG" sequcnes in each chromosome q of the genome
# This code assumes that hexameric sequences are present in the chromosome

from statsmodels.stats.multitest import fdrcorrection
from matplotlib.pyplot import figure
import matplotlib.pyplot as plt
from scipy.stats import norm
import numpy
import math

chr_coords = {"chr1q_start":28100001, "chr2q_start":96000000, "chr3q_start":94000000, "chr4q_start":50000001, "chr5q_start":51400000, \
           "chr6q_start":59800001, "chr7q_start":60100001, "chr8q_start":45200001, "chr9q_start":43000001, "chr10q_start":41600000, \
           "chr11q_start":53400001, "chr12q_start":35500001, "chr13q_start":18900001, "chr14q_start":17200001, "chr15q_start":19000001, \
           "chr16q_start":36800001, "chr17q_start":25100001, "chr18q_start":18500001, "chr19q_start":26200001, "chr20q_start": 28100001, \
           "chr21q_start":12000001, "chr22q_start":15000001}

chr_names = ["chr1.fa", "chr2.fa", "chr3.fa", "chr4.fa", "chr5.fa", "chr6.fa", "chr7.fa", "chr8.fa", "chr9.fa", "chr10.fa", "chr11.fa", "chr12.fa", \
             "chr13.fa", "chr14.fa", "chr15.fa", "chr16.fa", "chr17.fa", "chr18.fa", "chr19.fa", "chr20.fa", "chr21.fa", "chr22.fa"]

for name in chr_names:
    f = open(name)
    chr = ''
    for line in f:
        line = line.rstrip()
        if line[0] == ">":
            True
        else:
            chr = chr + line
    chr = chr.upper()
    chr_name = name[:-3] + "q_start"
    chrq = chr[chr_coords[chr_name]-1:] #Try chr[chr_coords[chr_name]-1:]
    chrq_len = len(chrq)
    print(name[:-3], "length is", chrq_len)
    print(" ")

    bin_size = 100000
    bin_num = math.ceil(chrq_len / bin_size)
        
    target1 = "TTAGGG"
    target2 = "CCCTAA"

    start_index = 0
    final_index = start_index + bin_size

    x_array = []
    for i in range(1, bin_num+1):
        x_array.append(i)

    target1_counts = []
    p_values = []
    sig_counts = {} #Statistically sigificant counts
    non_sig_counts = {}

    #Statistics
    standard_error = math.sqrt(((1/(4**6))*(1-(1/4**6)))/(bin_size))
    expected_value = int((1/4**6)*(bin_size)) 

    for i in range(bin_num):
        
        chrq_bin = chrq[start_index: final_index]
        
        target1_count = chrq_bin.count(target1)
        target1_counts.append(target1_count)
        
        z_score = ((target1_count/bin_size) - 1/(4**6))/standard_error
        p_value = norm.cdf(z_score) * 2
        p_values.append(p_value)
        
        start_index += bin_size
        final_index += bin_size

        if final_index > chrq_len:
            final_index = chrq_len +1

    significant, q_values = fdrcorrection(p_values)

    expected_line = []
    for i in range(bin_num):
        expected_line.append(expected_value)
        
    sig_counts = {}
    non_sig_counts = {}        
    for i in range(len(q_values)):
        if significant[i] == True:
            sig_counts[i+1] = target1_counts[i]
        else:
            non_sig_counts[i+1] = target1_counts[i]
    f_name = name[:-3] + "_results.txt"       
    f = open(f_name,'a')
    f.write("target1_counts - " + "p_values - " + "q_values - " + "significant\n")
    for i in range(len(p_values)):
        f.write(str(target1_counts[i]) + "\t" + str(p_values[i]) + "\t" + str(q_values[i]) + "\t" + str(significant[i])+ "\n")

    sig_array = sig_counts.keys()
    sig_values = sig_counts.values()
    non_sig_array = non_sig_counts.keys()
    non_sig_values = non_sig_counts.values()

    expected_line = []
    for i in range(bin_num):
        expected_line.append(expected_value)
    fig = plt.figure(figsize=(25, 14), dpi=150)
    plt.plot(x_array, expected_line)
    plt.plot(x_array, target1_counts, "g", linewidth=3.2)
    plt.scatter(sig_array, sig_values, marker='.', color='#5A5AFF', s = 170)
    plt.scatter(non_sig_array, non_sig_values, marker='^', color='#FF0000', s = 75)
    zip_object = zip(sig_array, sig_values)
    for i, j in zip_object:
        if j < expected_value:
            plt.text(i, j+0.2, "-", size = 20)
        elif j > expected_value:
            plt.text(i, j+0.2, "+")
        else:
            plt.text(i, j+0.2, "Nu")
    plt.xlabel('Window Number', fontsize=20)
    plt.ylabel('Number of TTAGGG telomeric repeats found in each window', fontsize=20)
    title = "A window distributon of TTAGGG telomeric repeats found in " + name[:-3]
    plt.title(title, fontsize=24)
    plot_name = name[:-3] + "q.png"
    fig.savefig(plot_name, dpi = 200)
    f.close()
