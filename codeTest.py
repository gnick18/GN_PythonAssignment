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



def translate_sequence(string_nucleotides):
        aa_seq_string = "" #this string will hole the translated sequence
        mito_table = CodonTable.unambiguous_dna_by_name["Vertebrate Mitochondrial"]
        end_translation = False # translation will end when the sequence runs out of nucleotides or a stop codon is found)
        beginCodon = 0 # inclusive counter for subsetting the string
        endCodon = 3 # exclusive counter for end of subset
        length = len(string_nucleotides)
        while(end_translation == False):
            codon = string_nucleotides[beginCodon:endCodon]
            #making sure the codon isn't a stop codon
            if(codon == "AGA" or codon == "AGG" or codon == "TAA" or codon == "TAG"):
                end_translation = True
                return(aa_seq_string)
            else:
                aminoAcid = mito_table.forward_table[codon]
                aa_seq_string = aa_seq_string + aminoAcid
                beginCodon = beginCodon + 3
                endCodon = endCodon + 3
                if(endCodon > length): #in this scenerio there are not enough nucleotides to continue translation, thus the loop ends
                    end_translation == True
        return(aa_seq_string)

def translate_sequence2(string_nucleotides):
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
    aa_seq_string = "" #this string will hole the translated sequence
    mito_table = CodonTable.unambiguous_dna_by_name["Vertebrate Mitochondrial"]
    codon = ""
    codonCounter = 0
    for i in range(len(string_nucleotides)): 
        #this loop will make the set of three nucleotides that make up the codon
        codon = codon + string_nucleotides[i]
        codonCounter = codonCounter + 1
        if(codonCounter == 3):
            codonCounter = 0
            #making sure the codon isn't a stop codon
            if(codon == "AGA" or codon == "AGG" or codon == "TAA" or codon == "TAG"):
                break
            else:
                aminoAcid = mito_table.forward_table[codon]
                aa_seq_string = aa_seq_string + aminoAcid
                codon = ""
    return(aa_seq_string)
    

#nucleotides = "AGTCCCCCCGGGAGCGAGC"
#AA2 = translate_sequence2(nucleotides)
#print(AA2)

##############################################################################

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
    AA_sequence = coding_dna.translate(to_stop = True, table=mito_table)
    AA_sequence_String = AA_sequence[0:len(AA_sequence)]
    return AA_sequence_String


#mito_table = CodonTable.unambiguous_dna_by_name["Vertebrate Mitochondrial"]
nucleotides = "AGTCCCCCCGGGAGCGAG"
#coding_dna = Seq(nucleotides, IUPAC.unambiguous_dna)
#AA_sequence = coding_dna.translate(table=mito_table, to_stop=True)
#AA_sequence_String = AA_sequence[0:len(AA_sequence)]
#print(AA_sequence_String)
#print(AA_sequence)
#AA = quick_translate(nucleotides)
#AA2 = translate_sequence2(nucleotides)
#print(AA)
#print(AA2)

###########################################################################################################

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

my_seq = "MAEGEITTFTALTEKFNLPPGNYKKPKLLYCSNGGHFLRILPDGTVDGTRDRSDQHIQLQLSAESVGEVYIKSTETGQYLAMDTSGLLYGSQTPSEECLFLERLEENHYNTYTSKKHAEKNWFVGLKKNGSCKRGPRTHYGQKAILFLPLPV"
print(compute_molecular_weight(my_seq))

####################################################################################################################

def GC_content(dna_sequence):
    """Description: finds the proportion of "G" and "C" within a DNA sequence.

    Parameters:
        dna_sequence : the dna sequence (must be a string object)

    Return: returns the proportion of "G" and "C" as a string
    """
    seq = Seq(dna_sequence, IUPAC.unambiguous_dna)
    gc = GC(seq)
    return gc

nucleotides = "AGAGAGAGCCCTTT"
gc = GC_content(nucleotides)
print(gc)

