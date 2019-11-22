# -*- coding: utf-8 -*-
"""
Created on Wed Nov 20 12:18:49 2019

@author: gnick
"""
from Bio import SeqIO
from Bio.Data import CodonTable
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
from Bio.SeqUtils import GC
from Bio.SeqUtils.ProtParam import ProteinAnalysis
import pandas as pd


def get_sequences_from_file(fasta_fn):
    """Description: converts fasta file to dictionary, storing the description, species_name, and sequence

    Parameters:
        fasta_fn  :  the fasta file containing the sequence and sample information

    Return: a dictionary that was created from fasta_fn

    Example of usage:
        >> dictionary = get_sequences_from_file(fasta_fn)

    Output(dictionary):
        [dictionary containing that sequences and some metadata]
    """
    sequence_data_dict = {}
    for record in SeqIO.parse(fasta_fn, "fasta"):
        description = record.description.split()
        species_name = description[1] + " " + description[2]
        sequence_data_dict[species_name] = record.seq
    return (sequence_data_dict)

cytb_seqs = get_sequences_from_file("penguins_cytb.fasta")

penguins_df = pd.read_csv("penguins_mass.csv") # Includes only data for body mass
species_list = list(penguins_df.species)

#Inspecting the Data
#print(penguins_df.dtypes)
#print(penguins_df.columns)
#print(penguins_df.head(10))
#print(penguins_df.shape)


#print(cytb_seqs)

#Object has 12 rows and 11 columns

#print(pd.unique(penguins_df['species']))
#print(pd.unique(penguins_df['mass']))

#adding columns with "NaN" as the value
penguins_df["mol_weight"] = "NaN"
#print(penguins_df.head(10))
penguins_df["GC_content"] = "NaN"
print(penguins_df.head(10))

for key, value in cytb_seqs.items():
    aa_seq = quick_translate(value)
    molWeight = compute_molecular_weight(aa_seq)
    gc = GC_content(value)
    #The values will now be filled into the correct column and row
    #filling in molecular weight
    row = penguins_df[penguins_df.species == key].index[0]
    penguins_df[row, ] = molWeight
