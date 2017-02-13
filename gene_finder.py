# -*- coding: utf-8 -*-
"""
YOUR HEADER COMMENT HERE

@author: Christina Gee

"""

import random
from amino_acids import aa, codons, aa_table   # you may find these useful
from load import load_seq


def shuffle_string(s):
    """Shuffles the characters in the input string
        NOTE: this is a helper function, you do not
        have to modify this in any way """
    return ''.join(random.sample(s, len(s)))

# YOU WILL START YOUR IMPLEMENTATION FROM HERE DOWN ###


def get_complement(nucleotide):
    """ Returns the complementary nucleotide

        nucleotide: a nucleotide (A, C, G, or T) represented as a string
        returns: the complementary nucleotide
    >>> get_complement('A')
    'T'
    >>> get_complement('C')
    'G'
    """
    # TODO: implement this
    if (nucleotide == "A"):
        return 'T'
    if (nucleotide == 'T'):
        return 'A'
    if (nucleotide == 'C'):
        return 'G'
    if (nucleotide == 'G'):
        return 'C'


def get_reverse_complement(dna):
    """ Computes the reverse complementary sequence of DNA for the specfied DNA
        sequence

        dna: a DNA sequence represented as a string
        returns: the reverse complementary DNA sequence represented as a string
    >>> get_reverse_complement("ATGCCCGCTTT")
    'AAAGCGGGCAT'
    >>> get_reverse_complement("CCGCGTTCA")
    'TGAACGCGG'
    """
    # TODO: implement this
    newDNA = ''
    for nucleo in dna:
        newNucleo = get_complement(nucleo)
        newDNA = newDNA + newNucleo
    return newDNA [: :-1] #-1 is used to go backwards


def rest_of_ORF(dna):
    """ Takes a DNA sequence that is assumed to begin with a start
        codon and returns the sequence up to but not including the
        first in frame stop codon.  If there is no in frame stop codon,
        returns the whole string.

        dna: a DNA sequence
        returns: the open reading frame represented as a string
    >>> rest_of_ORF("ATGTGAA")
    'ATG'
    >>> rest_of_ORF("ATGAGATAGG")
    'ATGAGA'
    """

    ORF=''
    end_codon = False
    codons= [dna[i:i+3] for i in range (0,len(dna),3)]
    stop_codons = ["TAG","TAA","TGA"]
    for codon in codons:
        if codon in stop_codons:
            break
        else:
            ORF += codon

    return ORF
    # x = 1
    # for i in dna: #for every letter in the DNA
    #     ORF = ORF + i #the letter is going to be added to the variable codon
    #     if len(ORF)== (3 * x):
    #         if ORF[-3:] == "TAA" or ORF[-3:]=="TAG" or ORF[-3:]=="TGA":
    #             return dna[ : len(ORF)-3] #[Start:End:Steps]
    #         x = x + 1



def find_all_ORFs_oneframe(dna):
    """ Finds all non-nested open reading frames in the given DNA
        sequence and returns them as a list.  This function should
        only find ORFs that are in the default frame of the sequence
        (i.e. they start on indices that are multiples of 3).
        By non-nested we mean that if an ORF occurs entirely within
        another ORF, it should not be included in the returned list of ORFs.

        dna: a DNA sequence
        returns: a list of non-nested ORFs
    >>> find_all_ORFs_oneframe("ATGCATGAATGTAGATAGATGTGCCC")
    ['ATGCATGAATGTAGA', 'ATGTGCCC']
    """
    #
    # returnList = []
    # codon = []
    # index = 0
    #
    # while index <= len(dna) - 1:
    #     codon.append(dna[index])
    #
    #     if dna[index] == "T" and index + 2 <= len(dna) - 1:
    #         if dna[index + 1] == "A":
    #             if dna[index + 2] == "A" or dna[index + 2] == "G":
    #                 # print(codon[:-1])
    #                 returnList.append("".join(codon[:-1]))
    #                 codon = []
    #                 index += 2
    #
    #         elif dna[index + 1] == "G":
    #             if dna[index + 2] == "A":
    #                 # print(codon[:-1])
    #                 returnList.append("".join(codon[:-1]))
    #                 codon = []
    #                 index += 2
    #
    #     index += 1
    #
    # returnList.append("".join(codon))
    #
    # return returnList
    # ATTEMPT #2
    i=0
    newdna=[]
    while (i<len(dna)): #i is less than the entire strand, i is a number
        if dna[i:i+3] == "ATG":
            new_orf = rest_of_ORF(dna[i:])
            i = i + len(new_orf)
            newdna.append(new_orf)
        else:
            i += 3
    return newdna


def find_all_ORFs(dna):
    """ Finds all non-nested open reading frames in the given DNA sequence in
        all 3 possible frames and returns them as a list.  By non-nested we
        mean that if an ORF occurs entirely within another ORF and they are
        both in the same frame, it should not be included in the returned list
        of ORFs.

        dna: a DNA sequence
        returns: a list of non-nested ORFs

    >>> find_all_ORFs("ATGCATGAATGTAG")
    ['ATGCATGAATGTAG', 'ATGAATGTAG', 'ATG']
    """
    i=0
    newdna=[]
    while i<3:
            newdna = newdna + find_all_ORFs_oneframe (dna[i:])
            i +=1
    return newdna


#offsets

def find_all_ORFs_both_strands(dna):
    """ Finds all non-nested open reading frames in the given DNA sequence on both
        strands.

        dna: a DNA sequence
        returns: a list of non-nested ORFs
    >>> find_all_ORFs_both_strands("ATGCGAATGTAGCATCAAA")
    ['ATGCGAATG', 'ATGCTACATTCGCAT']
    """
    # TODO: implement this
    pass

    all_orfs=find_all_ORFs(dna) #finding all ORFs for DNA
    all_orfs += find_all_ORFs(get_reverse_complement(dna))
    return all_orfs

#call fnd_all_ORFs_one frame - forward and backward DNA

def longest_ORF(dna):
    """ Finds the longest ORF on both strands of the specified DNA and returns it
        as a string
    >>> longest_ORF("ATGCGAATGTAGCATCAAA")
    'ATGCTACATTCGCAT'
    """
    # TODO: implement this
    all_orfs = find_all_ORFs_both_strands(dna)
    longest = max(all_orfs,key=len)
    return longest


def longest_ORF_noncoding(dna, num_trials):
    """ Computes the maximum length of the longest ORF over num_trials shuffles
        of the specfied DNA sequence

        dna: a DNA sequence
        num_trials: the number of random shuffles
        returns: the maximum length longest ORF """
    # TODO: implement this
    i=0
    greatest_length = 0
    while i < num_trials:
        shuffle = shuffle_string(dna)
        longest = longest_ORF(shuffle)
        if greatest_length < len(longest):
            greatest_length = len(longest)
        i += 1

    return greatest_length



def coding_strand_to_AA(dna):
    """ Computes the Protein encoded by a sequence of DNA.  This function
        does not check for start and stop codons (it assumes that the input
        DNA sequence represents an protein coding region).

        dna: a DNA sequence represented as a string
        returns: a string containing the sequence of amino acids encoded by the
                 the input DNA fragment

        >>> coding_strand_to_AA("ATGCGA")
        'MR'
        >>> coding_strand_to_AA("ATGCCCGCTTT")
        'MPA'
    """
    # TODO: implement this
    translated_aa = []
    i=0
    while i < len(dna):
        codeon = dna[i:i+3]
        if len(codeon) == 3:
            aa = aa_table[codeon]
            translated_aa += aa
            i += 3
        else:
            break
    return translated_aa

    #take DNA stran, make it to triplets and translate it to amino acids
    #look up dictionary values  - data structure for programming


def gene_finder(dna):
    """ Returns the amino acid sequences that are likely coded by the specified dna

        dna: a DNA sequence
        returns: a list of all amino acid sequences coded by the sequence dna.
    """
    # TODO: implement this
    pass

if __name__ == "__main__":
    import doctest
    doctest.testmod()
