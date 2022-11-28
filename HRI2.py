
# Notes:
# Written in Python 3.8 base.

from Bio.Seq import Seq
from Bio import SeqIO

# 1: Import the yeast genome.
def import_fasta(filename):
    '''Opens fasta file and saves it as dictionary with chromosome number as a key and chromosome sequence as the
    value.'''
    file = open(filename, 'r')                      # Open file in read mode.
    name = file.readlines()                         # Read lines in file and save them under variable.
    file.close()                                    # Close file.
    dictionary = {}                                 # Create empty dictionary.
    for string in name:                             # For every line in variable.
        if string[0] == ">":                        # If line starts with >.
            string_p = string.strip('\n')           # Remove trailing whitespace containing \n.
            chr = string_p.replace('>', '')         # Replace > with nothing (returns header without >).
            dictionary[chr] = ""                    # Add chr as key in dictionary and specify value as empty string.
        else:
            string_s = string.strip('\n')           # Remove trailing whitespace containing \n.
            value = dictionary[chr]                 # Call key for chr to give the value, give value a variable name.
            endvalue = value + string_s             # Add value to the new string, update the dictionary.
            dictionary[chr] = "".join(endvalue)     # Dictionary[key] = value | add value to last entered key.
    return dictionary


# 2: Import the S3 material.
def import_txt(filename):
    '''Opens text file and filters out lines that does not map to a chromosome, if a chromosome is listed with roman
    numerals, it is converted to a number, e.g. chrIII -> chr3.'''
    file = open(filename, 'r')                      # Open file in read mode.
    s1 = file.readlines()                           # Read lines in file and save them under variable.
    file.close()                                    # Close file.
    s2 = []                                         # Create new empty list.
    chrR = {'chrI': 'chr1', 'chrII': 'chr2', 'chrIII': 'chr3', 'chrIV': 'chr4', 'chrV': 'chr5', 'chrVI': 'chr6', 'chrVII': 'chr7', 'chrVIII': 'chr8', 'chrIX': 'chr9', 'chrX': 'chr10', 'chrXI': 'chr11', 'chrXII': 'chr12',
            'chrXIII': 'chr13', 'chrXIV': 'chr14', 'chrXV': 'chr15', 'chrXVI': 'chr16', 'chrXVII': 'chr17', 'chrXVIII': 'chr18', 'chrXIX': 'chr19', 'chrXX': 'chr20', 'chrXXI': 'chr21', 'chrXXII': 'chr22',
            'chrXXIII': 'chr23'}                    # Dictionary for roman numerals.
    for line in s1:                                 # For every line in list of strings.
        temp1 = line.strip('\n')                    # Remove trailing whitespace.
        temp2 = temp1.split('\t')                   # Split strings by tabs.
        s2.append(temp2)                            # Append line to s2.
    for line in s2:                                 # For every line in list of strings.
        if line[0] in chrR:                         # If string in index 0 is a key in dict chrR.
            value = chrR[line[0]]                   # Call key for chrR to give the value, give value a variable name.
            line[0] = value                         # Replace key with value.
    return s2


def extraction(genome, coordinates, flank_nucleotide_number):
    '''Takes coordinates from start and end point, uses it to extract nucleotides from corresponding position in
    genome file, pairs them for alignment. Genome is outut from import_fasta, coordinates is output from import_txt,
    number of nucleotides is the amount of nucleotides to be included upstream and downstream of start and end.'''
    coor = coordinates                          # Rename coordinates to coor.
    num = flank_nucleotide_number               # Rename number_of_nucleotides to num.
    results = []                                # Create empty list.
    for x in coor:                              # For every row (circle) in coordinates file.
        position = x[0]                         # Takes out chromosome number, gives it variable name.
        if position in genome:                  # If chromosome name is in the genome file. call place position
            start = int(x[1])                   # Takes out start position of eccDNA, gives it variable name.
            end = int(x[2])                     # Takes out end position of eccDNA, gives it variable name.
            seq = genome.get(position)          # Takes out the sequence for the chromosome the eccDNA is from.
            circle_seq = seq[start-1-num:end-1+num]  # Take subset of string from chromosome, that contains the
            # eccDNA +- the specified number of nucleotides upstream and downstream.
            if len(circle_seq) >= num*4:        # If the length of circle sequence (including flanks) >= flanks*4.
                up = circle_seq[:num*2]         # Subset the upstream flank of the eccDNA, twice the length of flank.
                down = circle_seq[-num*2:]      # Subset the downstream flank of the eccDNA, twice the length of flank.
                results.append([position, start, end, up, down])  # Append chromosome name to list.
            else:                               # If the length of circle sequence (including flanks) < flanks*4.
                hl = int(len(circle_seq)/2)     # Takes length of sequence, divides by 2, then converts to integer.
                up = circle_seq[:hl]            # Subset the upstream flank of the eccDNA.
                down = circle_seq[hl:]          # Subset the downstream flank of the eccDNA.
                results.append([position, start, end, hl, up, down])  # Append chromosome name to list.
        else:                                   # If chromosome name is not in the genome file.
            results.append([position])          # Append chromosome name to list.
    return results

