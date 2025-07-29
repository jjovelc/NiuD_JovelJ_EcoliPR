from BCBio import GFF
from Bio import SeqIO

input_file = "NZ_HG941718.1.gbk"
output_file = "NZ_HG941718.1.gff"

with open(input_file) as gbk_file, open(output_file, "w") as gff_file:
    GFF.write(SeqIO.parse(gbk_file, "genbank"), gff_file)
