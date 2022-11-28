
# Notes:
# Written in Python 3.8 base.
#
# Script requires:
# 1) Coordinates file: 1st column contains chromosome of origin in the form 'chr3' or 'chrIII', otherwise row is
# filtered out if key or value is not in dictionary chrR.
# 2) Coordinates file: 2nd and 3rd column contains start and end coordinate of eccDNA circles.
# 3) Organism can max have 23 chromosomes in the genome, unless more is added to the dictionary chrR.


def extract_circle_sequences(genome, coordinates, flank_nucleotide_number):
    '''Input genome.fa, return dictionary. Input coordinates.txt, return list. Input dictionary, list and
    flank_nucleotide_number, return Query.fa with upstream circle sequences, Subject.fa with downstream circle
    sequences and Remaining.txt with the circles not found in the genome.'''
    # Step 1: Input genome.fa, return dictionary.
    file = open(genome, 'r')                        # Open file in read mode.
    arg1 = file.readlines()                         # Read lines in file and save them under variable name arg1.
    file.close()                                    # Close file.
    genome = {}                                     # Create empty dictionary called genome.
    for string in arg1:                             # For every line in arg1.
        s_string = string.strip('\n')               # Remove trailing whitespace containing \n.
        if s_string[0] == ">":                      # If line starts with >, do the following.
            chr = s_string.replace('>', '')         # Replace > with nothing (title without >).
            genome[chr] = ""                        # Add chr as key in dictionary and specify value as empty string.
        else:                                       # If line does not begin with >, do the following.
            value = genome[chr]                     # Call key for chr to give the value, name it value.
            endvalue = value + s_string             # Add value to the new string, update the genome dict.
            genome[chr] = "".join(endvalue)         # Dictionary[key] = value | add value to last entered key.
    # Step 2: Input coordinates.txt, return list and file with coordinates.
    updated_coordinates = open("Coor.txt", "a")     # Open file, (will create file, if not found).
    file = open(coordinates, 'r')                   # Open file in read mode.
    s1 = file.readlines()                           # Read lines in file and save them under variable.
    file.close()                                    # Close file.
    coordinates = []                                # Create new empty list.
    chrR = {'chrI': 'chr1', 'chrII': 'chr2', 'chrIII': 'chr3', 'chrIV': 'chr4', 'chrV': 'chr5', 'chrVI': 'chr6',
            'chrVII': 'chr7', 'chrVIII': 'chr8', 'chrIX': 'chr9', 'chrX': 'chr10', 'chrXI': 'chr11', 'chrXII': 'chr12',
            'chrXIII': 'chr13', 'chrXIV': 'chr14', 'chrXV': 'chr15', 'chrXVI': 'chr16', 'chrXVII': 'chr17',
            'chrXVIII': 'chr18', 'chrXIX': 'chr19', 'chrXX': 'chr20', 'chrXXI': 'chr21', 'chrXXII': 'chr22',
            'chrXXIII': 'chr23'}                    # Dictionary for roman numerals.
    for line in s1:                                 # For every line in list of strings.
        temp1 = line.strip('\n')                    # Remove trailing whitespace.
        temp2 = temp1.split('\t')                   # Split strings by tabs.
        coordinates.append(temp2)                   # Append line to s2.
    for line in coordinates:                        # For every line in list of strings.
        if line[0] in chrR:                         # If string in index 0 is a key in dict chrR.
            value = chrR[line[0]]                   # Call key for chrR to give the value, give value a variable name.
            line[0] = value                         # Replace key with value.
    scoor = sorted(coordinates)                     # Sort coordinates, call the sorted coordiantes scoor..
    coor = []                                       # Make new empty list.
    coor.append(scoor[0])                           # Append first line from coordinates in coor.
    for x in scoor:                                 # For every line in coordinates.
        i = len(coor)-1                             # Define index for last line in coor.
        qa = (int(coor[i][1]))-25                   # Define query start range for combining circles.
        qb = (int(coor[i][1]))+25                   # Define query end range for combining circles.
        sa = (int(coor[i][2]))-25                   # Define subject start range for combining circles.
        sb = (int(coor[i][2]))+25                   # Define subject end range for combining circles.
        if qa < 0:                                  # If query circle start range is negative.
            qa = 0                                  # Reset to zero.
        if x[0] == coor[i][0] and qa <= int(x[1]) <= qb and sa <= int(x[2]) <= sb:  # If circle is in range of the last.
            a = coor[-1]                            # Set name of last line in coor to a.
            start = (int(a[1])+int(x[1]))//2        # Take the average value of start in this line and the last.
            end = (int(a[2])+int(x[2]))//2          # Take the average value of end in this line and the last.
            a[1] = str(start)                       # Replace start value in last line with average.
            a[2] = str(end)                         # Replace end value in last line with average.
        else:                                       # If circle is not in range.
            coor.append(x)                          # Append line to coor.
    coor2 = []                                      # Coordinates 2 list.
    for line in coor:                               # For line in scoor list.
        for value in line:                          # For every value in the line.
            coor2.append(value)                     # Append value to coor2 list.
            coor2.append('\t')                      # Append tab to coor2 list.
        coor2.append('\n')                          # Append newline to coor2 list.
    coor3 = "".join(coor2)                          # Join all strings together and call it coor3.
    updated_coordinates.write(coor3)                # Write coor3 into the updated coordinates file.
    # Step 3: Input dictionary, list and flank_nucleotide_number, return Query.fa with upstream circle sequences,
    # Subject.fa with downstream circle sequences and Remaining.txt with the circles not found in the genome.
    num = flank_nucleotide_number                   # Rename number_of_nucleotides to num.
    query = open("Query.fa", "a")                   # Open file, (will create file, if not found).
    subject = open("Subject.fa", "a")               # -||-
    remaining = open("Remaining.txt", "a")          # -||-
    for x in scoor:                                  # For every row (circle) in coor list.
        position = x[0]                             # Takes out chromosome, gives it variable name.
        if position in genome:                      # If chromosome name is in the genome file.
            start = int(x[1])                       # Takes out start position of eccDNA, gives it variable name.
            end = int(x[2])                         # Takes out end position of eccDNA, gives it variable name.
            seq = genome.get(position)              # Takes out the sequence for the chromosome the eccDNA is from.
            if start > num:                         # If start is bigger than the flank nucleotide number.
                circle_seq = seq[start-1-num:end-1+num]  # Take subset of string from chromosome, that includes the
                # eccDNA +- the specified number of nucleotides upstream and downstream.
            else:                                   # If start breakpoint is too close to start of chromosome.
                circle_seq = seq[0:end-1+num]       # Include up flank till start of chromosome.
            if len(circle_seq) >= num*4:            # If the length of circle sequence (including flanks) >= flanks*4.
                up = circle_seq[:num*2]             # Subset the upstream flank of the eccDNA, 2*length of flank.
                down = circle_seq[-num*2:]          # Subset the downstream flank of the eccDNA, 2*length of flank.
            else:                                   # If the length of circle sequence (including flanks) < flanks*4.
                hl = int(len(circle_seq)/2)         # Takes length of sequence, divides by 2, then converts to integer.
                up = circle_seq[:hl]                # Subset the upstream flank of the eccDNA.
                down = circle_seq[hl:]              # Subset the downstream flank of the eccDNA.
            title = [">", position, str(start), str(end)]  # Create unique title for eccDNA.
            joined_title = "_".join(title)          # Join list of strings to make a single string.
            query.write(joined_title + "\n")        # Append title and newline to query file.
            query.write(up + "\n")                  # Append upstream sequence and newline to query file.
            subject.write(joined_title + "\n")      # Append title and newline to subject file.
            subject.write(down + "\n")              # Append downstream sequence and newline to subject file.
        else:                                       # If chromosome name is not in the genome file.
            joined = "\t".join(x)                   # Join list of strings to make a single string.
            remaining.write(joined + "\n")          # Append string as line and newline to remaining file.
    updated_coordinates.close()                     # Close updated coordinates file.
    query.close()                                   # Close query file.
    subject.close()                                 # Close subject file.
    remaining.close()                               # Close remaining file.
    return


### Running python script, on Linux server, then running output through blast.
# To check if path is added to current directory, go to python shell/interpreter, execute;
# import sys
# print(sys.path)
# Then see if path is present.
# If path to current directory is not specified, in bash terminal, execute;
# To show current directory
# pwd
# export PYTHONPATH=/folder1/folder2/folder3
# Initiating python shell/interpreter, execute;
# python3
# Importing script into the shell;
# import HRI4
# Running function from script;
# HRI4.extract_circle_sequences("genome.fa", "coordinates.txt", flank_nucleotide_number)
# Exit pyhton shell/interpreter, execute;
# Ctrl + D
# Run blastn, execute;
# blastn -word_size 7 -query Query.fa -subject Subject.fa -out alignment.txt -outfmt 6

