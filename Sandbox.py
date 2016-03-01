__author__ = 'ilap'

import os.path

gene_name="RHO"

p = "Name="+gene_name+";"

print p+"ALmaaa....:"+ p[5:5+len(gene_name)]

from Bio.Align.Applications import ClustalwCommandline

clustalw2 = "./clustalw2"
assert os.path.isfile(clustalw2), "Clustal W executable missing"
cline = ClustalwCommandline(clustalw2, infile="./BBBB")
print cline

#help(ClustalwCommandline)
stdout, stderr = cline()

print stdout
print stderr