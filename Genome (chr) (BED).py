import numpy

try:
    f = open("chr21.fa")
except IOError:
    print("File chr21.fa does not exist!!")
    
chr21 = ''
for line in f:
    line = line.rstrip()
    if line[0] == ">":
        True
    else:
        chr21 = chr21 + line
    
print("chr21:")
target1 = "TTAGGG"
target2 = "CCCTAA"
target1_count = chr21.count(target1)
i = 0
if target1_count > 0:
    print(target1 + ": " + str(target1_count) + " instances")
    final_index = len(chr21)
    y = 0
    for i in range(target1_count):
        x = chr21.find(target1, y, final_index)
        y = x+6
        print("chr21" + "\t" + str(x) + "\t" + str(y))

else:
    print("No telomere sequences found in chr21.")
        
f.close()

#        seqs[name] = ''
#    else:
#        seqs[name] = seqs[name] + line

# int j iterates through the characters (nucleotide bases) in the line.
# i is used to add new character/nucleotide bases to the value field -
# of the chromosome key field.
       
        
