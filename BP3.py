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


def pathway_sort(align7, align9, query, coordinates, flank_nucleotide_number, allowed_breakpoint_shift):
    '''Sort results into HR, MMEJ and NHEJ.'''
    abs = allowed_breakpoint_shift                  # Rename allowed breakpoint shift.
    rep_path = open("Pathways_blast+.txt", "a")     # Open file, (will create file, if not found).
    path_coor = open("Pathway_coor.txt", "a")       # Open file for HR circles.
    hr = open("HR.txt", "a")                        # Open file for HR circles.
    mmej = open("MMEJ.txt", "a")                    # Open file for MMEJ circles.
    nhej = open("NHEJ.txt", "a")                    # Open file for NHEJ circles.
    excluded = open("Excluded.txt", "a")            # Open file for O excl. from HR, due to inverted/shifted repeats.
    ei7 = 0                                         # Set variable ei7 to 0, for excluded inverted repeats from new7.
    ei9 = 0                                         # Set variable ei9 to 0, for excluded inverted repeats from clean9.
    sb = 0                                          # Set variable ehs to 0, for shifted breakpoint.
# Step 1: Remove alignments for circles not aligned with themselves.
    clean7 = cleanup(align7)                        # Call cleanup function for align7.
    clean9 = cleanup(align9)                        # Call cleanup function for align9.
    for line in clean9:
        print(line)
# Step 2: Remove results from clean7 list if found in clean9 list to remove doubles.
    new7 = []                                       # New list for clean7 with doubles from clean9 removed.
    for x in clean7:                                # For every line in clean7.
        if x not in clean9:                         # If line in clean7 is not in clean9.
            new7.append(x)                          # Append line to the new7 list.
# Step : Remove alignments found inside other alignments.
    mclean7 = []                                    # New list with circles only represented once.
    m7 = {}
    for element in new7:                            # For every element in new7 list.
        if element[0] not in m7:                  # If the element is not already in the list.
            m7[element[0]] = ""
            mclean7.append(element)                 # Append element to mclean7 list.
    print(len(clean9))
    mclean9 = []                                    # New list with circles only represented once.
    m9 = {}
    for element in clean9:                          # For every element in clean9 list.
        if element[0] not in m9:                  # If the element is not already in the list.
            m9[element[0]] = ""
            mclean9.append(element)                 # Append element to mclean9 list.
    print(len(mclean9))
# Step 3: Separate direct from inverted repeats in new7 and clean9, put inverted repeats in Excluded file, and put
# direct repeats for mmej in the MMEJ file and repair pathway file, and put direct repeats in dr_9.
    u_titles = []                                   # New list for unsorted titles.
    for x in mclean7:                               # For every line in mclean7 list.
        if int(x[6]) <= int(x[7]) and int(x[9]) <= int(x[10]):  # If start <= end for query and subject.
            rep_path.write('MMEJ' + "\t")           # Write category to repair pathway file.
            for value in x:                         # For every value in the line.
                mmej.write(str(value) + "\t")       # Write line and tab to HR file.
                rep_path.write(str(value) + "\t")   # Write line and tab to repair pathway file.
            mmej.write("\n")                        # Write newline to MMEJ file.
            rep_path.write("\n")                    # Write newline to repair pathway file.
            y = 'MMEJ' + x[0]                       # Join pathway to string.
            u_titles.append(y)                      # Append line to unfiltered titles list.
        else:                                       # If start > end for query and/or subject.
            temp = []                               # Make temporary list to keep line in.
            temp.append(x)                          # Append line to temporary list.
            if temp[0] != x[0]:                     # If the index 0 in temp does not match index 0 in the line.
                excluded.write("Inverted repeat\tAlign7")  # Write category and source to excluded file.
                for value in x:                     # For every value in x.
                    excluded.write("\t")            # Write tab in excluded file.
                    excluded.write(value)           # Write value in excluded file.
                ei7 += 1                            # Add one to variable ei7.
                excluded.write("\n")                # Write newline in excluded file.
    dr_9 = []                                       # Direct repeat from align9.
    for x in mclean9:                               # For every line in mclean9.
        if int(x[6]) <= int(x[7]) and int(x[9]) <= int(x[10]):  # If start <= end for query and subject.
            dr_9.append(x)                          # Append to hr direct repeat list.
        else:                                       # If start > end alignment for query and/or subject.
            temp = []                               # Make temporary list to keep line in.
            temp.append(str(1))
            #print(x[0])
            if x[0] not in temp[0]:  # If the index 0 in temp does not match index 0 in the line.
                excluded.write("Inverted repeat\tAlign9")  # Write category and source to excluded file.
                for value in x:                     # For every value in x.
                    excluded.write("\t")            # Write tab in excluded file.
                    excluded.write(value)           # Write value in excluded file.
                ei9 += 1                            # Add one to variable ei9.
                excluded.write("\n")                # Write newline in excluded file.
            for value in x:
                temp.append(value)                          # Append line to temporary list.
# Step 4: In list with direct repeats, account for breakpoint shift.
    for x in dr_9:                                  # For every line in direct repeats 9 list.
        bq = flank_nucleotide_number                # Set query break point to flank_nucleotide_number.
        bs = int(x[8]) - bq                         # Set subject breakpoint to length - flank_nucleotide_number.
        alq = bq - int(x[6])                        # Find difference of breakpoint - start of alignment for query.
        als = bs - int(x[9])                        # Find difference of breakpoint - start of alignment for subject.
        b_diff = alq-als                            # Find difference between breakpoint alignment.
        if abs >= b_diff >= -abs:                   # If breakpoint difference is >= -abs and <= abs.
            rep_path.write('HR' + "\t")             # Write category to repair pathway file.
            for value in x:                         # For every value in the line.
                hr.write(str(value) + "\t")         # Write line and tab to HR file.
                rep_path.write(str(value) + "\t")   # Write line and tab to repair pathway file.
            hr.write("\n")                          # Write newline in HR file.
            rep_path.write("\n")                    # Write newline in repair pathway file.
            y = 'HR' + x[0]                         # Join pathway to string.
            u_titles.append(y)                      # Append line to unfiltered titles list.
        else:                                       # If alignment is shifted.
            excluded.write("Shifted direct repeat, Align9, " + str(x) + "\n")  # Append category, source and string.
            sb += 1                                 # Add one to variable sb.
# Step 5: Making the NHEJ file.
    t_nhej = []                                     # Make new temporary nhej list.
    file = open(query, 'r')                         # Open query file in read mode.
    query = file.readlines()                        # Read lines in file and save them under variable.
    file.close()                                    # Close file.
    for x in query:                                 # For every line in query.
        if x[0] == ">":                             # If line starts with >.
            string = x.strip('\n')                  # Remove trailing whitespace containing \n.
            title = string.replace('>', '')         # Replace > with nothing (returns header without >).
            t_nhej.append(title)                    # Append title to temp_nhej list.
    titles = []                                     # Titles of HR and MMEJ including Excluded.
    for line in mclean7:                            # For every line in list.
        titles.append(line[0])                      # Append first value to titles.
    for line in mclean9:                            # For every line in list.
        titles.append(line[0])                      # Append first value to titles.
    new = [x for x in t_nhej if x not in titles]    # For every x in t_nhej, if it is not in titles, append to new.
    for x in new:                                   # For every value in new.
        nhej.write(x + "\n")                        # Write title and newline to nhej file.
        rep_path.write('NHEJ' + "\t" + str(x) + "\n")  # Write line and newline to repair pathway file.
        y = 'NHEJ' + x                              # Join pathway to string.
        u_titles.append(y)                          # Append line to filtered titles list.
# Step 6: Removing repeating values in unfiltered list.
    f_titles = []                                   # New list for titles without multiples.
    for element in u_titles:                        # For every element in unsorted titles.
        if element not in f_titles:                 # If the element is not already in the list.
            f_titles.append(element)                # Append element to filtered titles.
# Step 7: Classifying circles by repair pathway.
    all_titles = []                                 # Titles for all circles except Excluded and Remaining.
    for x in f_titles:                              # For every line in filtered titles list.
        split = x.split("_")                        # Split string into chrom, start, end.
        all_titles.append(split)                    # Append titles to all titles list.
    file = open(coordinates, 'r')                   # Open file in read mode.
    s1 = file.readlines()                           # Read lines in file and save them under variable.
    file.close()                                    # Close file.
    co = []                                         # Create new empty list.
    chrR = {'chrI': 'chr1', 'chrII': 'chr2', 'chrIII': 'chr3', 'chrIV': 'chr4', 'chrV': 'chr5', 'chrVI': 'chr6',
            'chrVII': 'chr7', 'chrVIII': 'chr8', 'chrIX': 'chr9', 'chrX': 'chr10', 'chrXI': 'chr11', 'chrXII': 'chr12',
            'chrXIII': 'chr13', 'chrXIV': 'chr14', 'chrXV': 'chr15', 'chrXVI': 'chr16', 'chrXVII': 'chr17',
            'chrXVIII': 'chr18', 'chrXIX': 'chr19', 'chrXX': 'chr20', 'chrXXI': 'chr21', 'chrXXII': 'chr22',
            'chrXXIII': 'chr23'}                    # Dictionary for roman numerals.
    for line in s1:                                 # For every line in list of strings.
        temp1 = line.strip('\n')                    # Remove trailing whitespace.
        temp2 = temp1.split('\t')                   # Split strings by tabs.
        co.append(temp2)                            # Append line to s2.
    for line in co:                                 # For every line in list of strings.
        if line[0] in chrR:                         # If string in index 0 is a key in dict chrR.
            value = chrR[line[0]]                   # Call key for chrR to give the value, give value a variable name.
            line[0] = value                         # Replace key with value.
    sco = []                                        # Make new list for sorted coordinates.
    for lc in co:                                   # For every line in coordinates.
        for lt in all_titles:                       # For every line in all titles.
            if lc[0] == lt[1] and lc[1] == lt[2] and lc[2] == lt[3]:  # If lc index 0, 1 & 2 match lt index 1, 2 & 3.
                sco.append(lt[0])                   # Append pathway for circle to sorted coordinates list.
                for value in lc:                    # For every value in lc.
                    sco.append("\t")                # Append tab before value in sorted coordinates.
                    sco.append(value)               # Append value to sorted coordinates.
                sco.append("\n")                    # Append newline at the end of every line in sorted coordinates.
    s = "".join(sco)                                # Join all strings into one string and rename list to s.
    path_coor.write(s)                              # Write s list to path coor file.
# Step 8: Counting number of circles belonging to each pathway as well as percentages.
    count = []                                      # Make new empty list.
    x = 0                                           # Set variable x to 0.
    y = 0                                           # Set variable y to 0.
    z = 0                                           # Set variable z to 0.
    for line in all_titles:                         # For every line in all titles list.
        if line[0] == "MMEJ":                       # If circle is mmej.
            x += 1                                  # Add one to the count.
        elif line[0] == "HR":                       # If circle is hr.
            y += 1                                  # Add one to the count.
        elif line[0] == "NHEJ":                     # If circle is nhej.
            z += 1                                  # Add one to the count.
    print(ei7+ei9+sb)
    t = x+y+z                                       # Add all circles to count total amount.
    count.append("Total number of circles: " + str(t) + "\n")
    count.append("Pathway count:\n")
    count.append("MMEJ; " + str(x) + ". HR; " + str(y) + ". NHEJ; " + str(z) + ".\n")
    count.append("Percentages:\n")
    count.append("MMEJ; " + str(round(x/t*100, 2)) + ". HR; " + str(round(y/t*100, 2)) + ". NHEJ; " + str(round(z/t*100, 2)) + ".\n")
    rep_path.close()                                # Close file.
    path_coor.close()                               # Close file.
    hr.close()                                      # Close file.
    mmej.close()                                    # Close file.
    nhej.close()                                    # Close file.
    excluded.close()                                # Close file.
    return

