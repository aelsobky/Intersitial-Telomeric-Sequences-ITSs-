import pyranges as pr
import os

# Fetching the path of chr_all.bed (the file that has blastn query results)
path = os.path.abspath("chr_all.bed")
# converting the file into a pyrange object called (chr_all)
chr_all = pr.read_bed(path, False, nrows=1258)
# Fetching the path of chr_published.bed (the file that has the published its loci)
path = os.path.abspath("chr_published.bed")
chr_published = pr.read_bed(path, False, nrows=125)
# converting the file into a pyrange object called (chr_published)
print(chr_all.intersect(chr_published))




