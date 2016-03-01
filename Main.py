__author__ = 'ilap'

import os.path
from Bio import SwissProt
import gzip
from Bio.Align.Applications import ClustalwCommandline
from Bio import AlignIO
from Bio import Phylo

'''
Find gene of Rhodopsin protein in Uniprot
    = Manually download database from Uniprot.
    = Search for rho in database and download them and convert to Fasta.
Alignment
    = MSA? Clustal or TCoffeee or Muscle.
Consensus
    =
Visualise
    = Tree
    = Heatmap
'''

class Uniprot():
    def __init__(self, database_file):
        try:
            self.fasta_base = "found_genes"
            self.fasta_file = self.fasta_base+".fasta"
            self.fd = gzip.open(database_file,mode='rb')
            #self.fd = open(database_file,mode='r')
        except Exception as e:
            print "Cannot open and parse the \"{}\" file".format(database_file)
            exit(1)
                
    def __del__(self):

        self.fd.close()

    def get_genes (self,gene_name=""):
        if gene_name != "":
            print "Finding \"{}\" gene in Uniprot database...".format(gene_name)
            upper_name = gene_name.upper() # Rho --> RHO

            output_handle = open(self.fasta_file, "w")

            for record in SwissProt.parse (self.fd):

                match = record.gene_name[5:5+len(upper_name)+1].upper()
                # Name=Rhodop; --> RHOD (Length of the queried name (rho)+1)
                # For matching the two possibilities
                # 1) Name=Rho;
                # 2) Name=rho {ECO.....}
                # So, it fill compare the queried gene name and match one e.g.
                # in 1st case "RHO " == "RHO;" or "RHO;" == "RHO;"
                # in 2nd case "RHO " == "RHO " or "RHO;" == "RHO "
                # We do not consider gene names differ to "Name=...;" in swisprot file



                if (upper_name+" ") == match or (upper_name+";") == match:
                    print "Add protein to fasta file: " + record.entry_name + ", ...." + record.gene_name
                    output = ">"+record.entry_name+"\n"+record.sequence.format("fasta")+"\n"
                    #print output
                    output_handle.write(output)
            output_handle.write("")
            output_handle.close()

    def align (self):
        clustalw2 = "./clustalw2"
        assert os.path.isfile(clustalw2), "Clustal W executable missing"

        cline = ClustalwCommandline(clustalw2, infile=self.fasta_file)
        print "Aligning fasta files.."
        stdout, stderr = cline ()

        align = AlignIO.read(self.fasta_base+".aln", "clustal")
        print align

        tree = Phylo.read(self.fasta_base + ".dnd", "newick")
        Phylo.draw_ascii(tree)


#### MAIN ##########
#print "AAAAAAA", os.path.curdir
uniprot = Uniprot ("./uniprot_sprot.dat.gz")

#uniprot.get_genes("rho")

uniprot.align()

