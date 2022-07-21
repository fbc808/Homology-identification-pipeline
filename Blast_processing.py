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


def pathway_sort(align_mmej, align_hr, query_nhej, flank_nucleotide_number):
    '''Sort results into HR, MMEJ and NHEJ.'''
    # Step 1: Remove alignments for circles not aligned with themselves.
    clean_mmej = cleanup(align_mmej)            # Remove circles not aligned with themselves. print("mmej", len(clean_mmej))
    clean_hr = cleanup(align_hr)                # Remove circles not aligned with themselves. print("hr", len(clean_hr))
    # Step 2: Remove results from mmej list if found in hr list to avoid doubles.
    new_mmej = []                               # Make new empty list.
    for x in clean_mmej:                        # For every line in clean_mmej.
        if x not in clean_hr:                   # If line in clean_mmej is not in clean_hr.
            new_mmej.append(x)                  # Append line to the new_mmej list. print("MMEJ", len(new_mmej))
    # Step 3: Separate direct from inverted repeats, where homology is found, put inverted repeats in Excluded file.
    dr_mmej = []                                # Direct repeat microhomology-mediated end joining.
    dr_hr = []                                  # Direct repeat homologous recombination.
    s_hr = []                                   # Shifted alignment homologous recombination.
    ab_hr = []                                  # Aligned breakpoints homologous recombination.
    inverted_repeats = []                       # Inverted repeats.
    for x in new_mmej:                          # For every line in new_mmej.
        if x[6] < x[7] and x[9] < x[10]:        # If start < end for query and subject.
            dr_mmej.append(x)                   # Append to mmej direct repeat list.
        else:                                   # If start > end for query and/or subject.
            inverted_repeats.append(x)          # Append to inverted repeats. print("DR MMEJ", len(dr_mmej))
    for x in clean_hr:                          # For every line in clean_hr.
        if x[6] < x[7] and x[9] < x[10]:        # If start alignment > than end alignment for query and subject.
            dr_hr.append(x)                     # Append to hr direct repeat list.
        else:                                   # If start alignment < than end alignment for query and/or subject.
            inverted_repeats.append(x)          # Append to inverted repeats.
    # Step 4: In list with direct repeats, account for breakpoint shift.
    for x in dr_hr:                             # For every line in direct repeats hr list.
        bq = flank_nucleotide_number            # Set query break point to flank_nucleotide_number.
        bs = int(x[8]) - bq                     # Set subject breakpoint to length - flank_nucleotide_number.
        alq = bq - int(x[6])                    # Find difference of breakpoint - start of alignment for query.
        als = bs - int(x[9])                    # Find difference of breakpoint - start of alignment for subject.
        b_diff = alq-als                        # Find difference between breakpoint alignment.
        if b_diff <= 50 and b_diff >= -50:      # If breakpoint difference is > -30 and < 30.
            ab_hr.append(x)                     # Append to ab hr list.
        else:                                   # If alignment is shifted
            s_hr.append(x)                      # Append line to shift hr list.
    # Making the NHEJ file.
    temp_nhej = []                              # Make new empty list.
    file = open(query_nhej, 'r')                # Open file in read mode.
    query = file.readlines()                    # Read lines in file and save them under variable.
    file.close()                                # Close file.
    for x in query:                             # For every line in query_nhej.
        if x[0] == ">":                         # If line starts with >.
            string = x.strip('\n')              # Remove trailing whitespace containing \n.
            title = string.replace('>', '')     # Replace > with nothing (returns header without >).
            temp_nhej.append(title)             # Append title to temp_nhej list. print("nhej", len(temp_nhej))
    mmej_titles = []                            # Make new empty list.
    hr_titles = []                              # Make new empty list.
    for line in new_mmej:                       # For every line in list.
        mmej_titles.append(line[0])             # Append first value to mmej_titles.
    for line in clean_hr:                       # For every line in list.
        hr_titles.append(line[0])               # Append first value to hr_titles.
    new = [x for x in temp_nhej if x not in mmej_titles]    # For every x in temp_nhej, if it is not in mmej_titles, append to new list.
    newnew = [x for x in new if x not in hr_titles]         # For every x in new, if it is not in hr_titles, append to newnew list.
    # print("NHEJ", len(newnew)), print(ab_hr), print(dr_mmej), print(inverted_repeats)
    return

