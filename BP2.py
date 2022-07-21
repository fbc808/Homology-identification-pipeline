# Assumes 1st and 2nd column are titles of circles

def cleanup(blastn):
    '''Filters out circles not aligned with themselves.'''
    file = open(blastn, 'r')            # Open file in read mode.
    results = file.readlines()          # Read lines in file and save them under variable.
    file.close()                        # Close file.
    clean = []                          # Make new empty list.
    align = []                          # Make new empty list.
    for line in results:                # For every line in list of strings.
        temp1 = line.strip('\n')        # Remove trailing whitespace.
        temp2 = temp1.split('\t')       # Split strings by tabs.
        clean.append(temp2)             # Append line to clean list.
    for line in clean:                  # For every line in cleaned list.
        if line[0] == line[1]:          # If circle is aligned with itself.
            align.append(line)          # Append line to align list.
    return align                        # Return list of circles aligned with themselves.


def pathway_sort(align7, align9, query, flank_nucleotide_number, allowed_breakpoint_shift):
    '''Sort results into HR, MMEJ and NHEJ.'''
    abs = allowed_breakpoint_shift                  # Rename.
    rep_path = open("Repair_pathways.txt", "a")     # Open file, (will create file, if not found).
    hr = open("HR.txt", "a")                        # Open file.
    mmej = open("MMEJ.txt", "a")                    # Open file.
    nhej = open("NHEJ.txt", "a")                    # Open file.
    excluded = open("Excluded.txt", "a")            # Open file.
    # Step 1: Remove alignments for circles not aligned with themselves.
    clean7 = cleanup(align7)                        # Call cleanup function for align7.
    clean9 = cleanup(align9)                        # Call cleanup function for align9.
    # Step 2: Remove results from clean7 list if found in clean9 list to remove doubles.
    new7 = []                                       # Make new empty list.
    for x in clean7:                                # For every line in clean7.
        if x not in clean9:                         # If line in clean7 is not in clean9.
            new7.append(x)                          # Append line to the new7 list.
    # Step 3: Separate direct from inverted repeats, put inverted repeats in Excluded file.
    for x in new7:                                  # For every line in new7.
        if int(x[6]) < int(x[7]) and int(x[9]) < int(x[10]):  # If start < end for query and subject.
            mmej.write(str(x) + "\n")               # Append line and newline.
            rep_path.write("MMEJ, " + str(x) + "\n")  # Append line and newline to rep path file.
        else:                                       # If start > end for query and/or subject.
            excluded.write("Inverted repeat, Align7, " + str(x) + "\n")  # Append category, source and string.
    dr_9 = []                                       # Direct repeat from align9.
    for x in clean9:                                # For every line in clean9.
        if int(x[6]) < int(x[7]) and int(x[9]) < int(x[10]):  # If start > end for query and subject.
            dr_9.append(x)                          # Append to hr direct repeat list.
        else:                                       # If start alignment < than end alignment for query and/or subject.
            excluded.write("Inverted repeat, Align9, " + str(x) + "\n")  # Append category, source and string.
    # Step 4: In list with direct repeats, account for breakpoint shift.
    for x in dr_9:                                  # For every line in direct repeats 9 list.
        bq = flank_nucleotide_number                # Set query break point to flank_nucleotide_number.
        bs = int(x[8]) - bq                         # Set subject breakpoint to length - flank_nucleotide_number.
        alq = bq - int(x[6])                        # Find difference of breakpoint - start of alignment for query.
        als = bs - int(x[9])                        # Find difference of breakpoint - start of alignment for subject.
        b_diff = alq-als                            # Find difference between breakpoint alignment.
        if abs >= b_diff >= -abs:                   # If breakpoint difference is > -abs and < abs.
            hr.write(str(x) + "\n")                 # Append line and newline to HR file.
            rep_path.write('HR, ' + str(x) + "\n")  # Append line and newline to rep path file.
        else:                                       # If alignment is shifted.
            excluded.write("Shifted direct repeat, Align9, " + str(x) + "\n")  # Append category, source and string.
    # Step 5: Making the NHEJ file.
    temp_nhej = []                                  # Make new empty list.
    file = open(query, 'r')                         # Open file in read mode.
    query = file.readlines()                        # Read lines in file and save them under variable.
    file.close()                                    # Close file.
    for x in query:                                 # For every line in query.
        if x[0] == ">":                             # If line starts with >.
            string = x.strip('\n')                  # Remove trailing whitespace containing \n.
            title = string.replace('>', '')         # Replace > with nothing (returns header without >).
            temp_nhej.append(title)                 # Append title to temp_nhej list.
    titles7 = []                                    # Make new empty list.
    titles9 = []                                    # Make new empty list.
    for line in new7:                               # For every line in list.
        titles7.append(line[0])                     # Append first value to titles7.
    for line in clean9:                             # For every line in list.
        titles9.append(line[0])                     # Append first value to titles9.
    new = [x for x in temp_nhej if x not in titles7]  # For every x in temp_nhej, if it is not in titles7.
    newnew = [x for x in new if x not in titles9]   # For every x in new, if it is not in titles9.
    for x in newnew:                                # For every value in newnew.
        nhej.write(x + "\n")                        # Append title and newline to nhej file.
        rep_path.write('NHEJ, ' + str(x) + "\n")    # Append line and newline to rep path file.
    # Step 6: Outputting the predicted pathways
    rep_path.close()                                # Close file.
    hr.close()                                      # Close file.
    mmej.close()                                    # Close file.
    nhej.close()                                    # Close file.
    excluded.close()                                # Close file.
    return


#test = pathway_sort("Align7.txt", "Align9.txt", "Query.fa", 250, 50)
#print(test)

