# -*- coding: utf-8 -*-
"""
Created on Wed Nov 20 12:25:28 2019

@author: gnick
"""

from Bio import SeqIO
from Bio.Data import CodonTable
import pandas as pd

 
mito_table = CodonTable.unambiguous_dna_by_name["Vertebrate Mitochondrial"]
codon = "AGG"
aminoAcid = mito_table.forward_table[codon]
print(aminoAcid)