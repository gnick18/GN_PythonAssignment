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

def translate_sequence(string_nucleotides):
    """Description: translates a string of nucleotides

    Parameters:
        string_nucleotides  :  the string of nucleotides to be translated

    Return: returns the amino acid sequence as a string

    Example of usage:
      >> aa_seq_string = translate_sequence(string_nucleotides)
      >> print (aa_seq_string)

    Output:
        [amino acid sequence]
    """
    aa_seq_string = ""  # this string will hole the translated sequence
    mito_table = CodonTable.unambiguous_dna_by_name["Vertebrate Mitochondrial"]
    codon = ""
    codonCounter = 0
    for i in range(len(string_nucleotides)):
        # this loop will make the set of three nucleotides that make up the codon
        codon = codon + string_nucleotides[i]
        codonCounter = codonCounter + 1
        if (codonCounter == 3):
            codonCounter = 0
            # making sure the codon isn't a stop codon
            if (codon == "AGA" or codon == "AGG" or codon == "TAA" or codon == "TAG"):
                break
            else:
                aminoAcid = mito_table.forward_table[codon]
                aa_seq_string = aa_seq_string + aminoAcid
                codon = ""
    return (aa_seq_string)


def quick_translate(string_nucleotides):
    """Description: translates a string of nucleotides

        Parameters:
            string_nucleotides  :  the string of nucleotides to be translated

        Return: returns the amino acid sequence as a string

        Example of usage:
          >> aa_seq_string = quick_translate(string_nucleotides)
          >> print (aa_seq_string)

        Output (String):
            [amino acid sequence]
        """
    mito_table = CodonTable.unambiguous_dna_by_name["Vertebrate Mitochondrial"]
    coding_dna = Seq(string_nucleotides, IUPAC.unambiguous_dna)
    AA_sequence: Seq = coding_dna.translate(to_stop = True, table=mito_table)
    AA_sequence_String = AA_sequence._data
    return AA_sequence_String

def compute_molecular_weight(aa_seq):
    """Description: calculates the molecular weight from an amino acid sequence

    Parameters:
        aa_seq : the amino acid sequence (must be string object)

    Return: returns the molecular weight of the amino acid sequence

    Example of usage:
      >> MolecularWeight = compute_molecular_weight(aa_seq)
      >> print(MolecularWeight)

    Output (String):
        [molecular weight of amino acid sequence]
    """
    analysed_seq = ProteinAnalysis(aa_seq)
    weight = analysed_seq.molecular_weight()
    return weight

def GC_content(dna_sequence):
    """Description: finds the proportion of "G" and "C" within a DNA sequence.

    Parameters:
        dna_sequence : the dna sequence (must be a string object)

    Return: returns the proportion of "G" and "C" as a string
    """
    seq = Seq(dna_sequence, IUPAC.unambiguous_dna)
    gc = GC(seq)
    return gc

#############################################
cytb_seqs = get_sequences_from_file("penguins_cytb.fasta")

penguins_df = pd.read_csv("penguins_mass.csv") # Includes only data for body mass
species_list = list(penguins_df.species)

#adding columns with "NaN" as the value
penguins_df["mol_weight"] = "NaN"
#print(penguins_df.head(10))
penguins_df["GC_content"] = "NaN"
print(penguins_df.head(10))

#molWeight = 2
#i = penguins_df[penguins_df.species == 'Pygoscelis adeliae'].index[0]
#penguins_df.loc[[index, 'mol_weight']] = molWeight
#print(penguins_df)

for key, value in cytb_seqs.items():
    value_string = value._data
    aa_seq = quick_translate(value_string)
    molWeight = compute_molecular_weight(aa_seq)
    gc = GC_content(value_string)
    #The values will now be filled into the correct column and row
    row = penguins_df[penguins_df.species == key].index[0]
    #filling in molecular weight
    penguins_df.iloc[row, 2] = molWeight
    #filling in the GC content
    penguins_df.iloc[row, 3] = gc

print(penguins_df)

def get_sequence_length (seq):
    """Description: finds the length of the sequence

    Parameter:
        seq: the sequence object containing the sequence stored in it
    Return:
        the length of the sequence
    """
    sequenceString = seq._data
    length = len(sequenceString)
    return length

penguins_df['number_nucleotides'] = "NaN"
for key, value in cytb_seqs.items():
    length = get_sequence_length(value)
    #The values will now be filled into the correct column and row
    row = penguins_df[penguins_df.species == key].index[0]
    #filling in the number of nucleotides column
    penguins_df.iloc[row, 4] = length
    print(penguins_df)