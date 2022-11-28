# Assumes 1st and 2nd column are titles of circles

def cleanup(blastn):
    '''Filters out circles not aligned with themselves.'''
    file = open(blastn, 'r')                        # Open file in read mode.
    results = file.readlines()                      # Read lines in file and save them under variable.
    file.close()                                    # Close file.
    clean = []                                      # Make new empty list.
    align = []                                      # Make new empty list.
    for line in results:                            # For every line in list of strings.
        temp1 = line.strip('\n')                    # Remove trailing whitespace.
        temp2 = temp1.split('\t')                   # Split strings by tabs.
        clean.append(temp2)                         # Append line to clean list.
    for line in clean:                              # For every line in cleaned list.
        if line[0] == line[1]:                      # If circle is aligned with itself.
            align.append(line)                      # Append line to align list.
    return align                                    # Return list of circles aligned with themselves.


def no_multiples(list1, i):
    '''Removes lines with multiples from list of lists, based on specified index of column'''
    nl = []                                         # New list.
    titles = []                                     # Titles list.
    nl.append(list1[0])                             # Append first line from list1 with multiples to new list.
    titles.append(list1[0][i])                      # Append specified column from first line in list1.
    for line in list1:                              # For every line in list1.
        if line[i] not in titles:                   # If column in line is not in the column in the new list.
            nl.append(line)                         # Append line to new list.
            titles.append(line[i])                  # Append the indexed column to titles list.
    return nl                                       # Return new list without multiples.


def pathway_sort(align7, align9, query, coordinates, circle_coor, flank_nucleotide_number, allowed_breakpoint_shift,
                 gene_column):
    '''Sort results into HR, MMEJ and NHEJ.'''
    abs = allowed_breakpoint_shift                  # Rename allowed breakpoint shift.
    path_blast = open("Pathways_blast+.txt", "a")   # Open file, (will create file, if not found).
    path_coor = open("Pathway_coor.txt", "a")       # Open file for pathway coordinates for circles.
    excluded = open("Excluded.txt", "a")            # Open file for O excl. from HR, due to inverted/shifted repeats.
# Step 1: Remove alignments for circles not aligned with themselves.
    clean7 = cleanup(align7)                        # Call cleanup function for align7.
    clean9 = cleanup(align9)                        # Call cleanup function for align9.
# Step 2: Remove results from clean7 list if found in clean9 list to remove doubles.
    new7 = []                                       # New list for clean7 with doubles from clean9 removed.
    clean9_title = []                               # New list for clean9 titles.
    for x in clean9:                                # For every line in clean9.
        clean9_title.append(x[0])                   # Append title to clean9 title list.
    for x in clean7:                                # For every line in clean7.
        if x[0] not in clean9_title:                # If the title is not in clean9 title list.
            new7.append(x)                          # Append line to new7 list.
# Step 3: Separate direct from inverted repeats in new7 and clean9, put inverted repeats in Excluded file, and put
# direct repeats for mmej in the MMEJ file and repair pathway file, and put direct repeats in dr_9.
    u_titles = []                                   # New list for unsorted titles.
    inverted = []                                   # List for inverted repeats.
    for x in new7:                                  # For every line in new7 list.
        if int(x[6]) <= int(x[7]) and int(x[9]) <= int(x[10]):  # If start <= end for query and subject.
            path_blast.write('MMEJ' + "\t")         # Write category to repair pathway file.
            for value in x:                         # For every value in the line.
                path_blast.write(str(value) + "\t")  # Write line and tab to repair pathway file.
            path_blast.write("\n")                  # Write newline to repair pathway file.
            y = 'MMEJ' + x[0]                       # Join pathway to string.
            u_titles.append(y)                      # Append line to unfiltered titles list.
        else:                                       # If start > end for query and/or subject.
            temp = []                               # Make temporary list to keep line in.
            temp.append(x)                          # Append line to temporary list.
            inv = []                                # Make temp inverted list.
            if temp[0] != x[0]:                     # If the index 0 in temp does not match index 0 in the line.
                inv.append("Inverted repeat")       # Append category to temp inverted list.
                inv.append("Align7")                # Append source to temp inverted list.
                for value in x:                     # For every value in x.
                    inv.append(value)               # Append value to inverted temp list.
                inverted.append(inv)                # Append inverted temp list to inverted list.
    dr_9 = []                                       # Direct repeat from align9.
    for x in clean9:                                # For every line in clean9.
        if int(x[6]) <= int(x[7]) and int(x[9]) <= int(x[10]):  # If start <= end for query and subject.
            dr_9.append(x)                          # Append to hr direct repeat list.
        else:                                       # If start > end alignment for query and/or subject.
            temp = []                               # Make temporary list to keep line in.
            temp.append(x)                          # Append line to temporary list.
            inv = []                                # Make temp inverted list.
            if temp[0] != x[0]:                     # If the index 0 in temp does not match index 0 in the line.
                inv.append("Inverted repeat")       # Append category to temp inverted list.
                inv.append("Align9")                # Append source to temp inverted list.
                for value in x:                     # For every value in x.
                    inv.append(value)               # Append value to inverted temp list.
                inverted.append(inv)                # Append inverted temp list to inverted list.
    # Step 4: In list with direct repeats, account for breakpoint shift.
    fnn = flank_nucleotide_number                   # Rename flank nucleotide number to fnn.
    pb = []                                         # Path blast+ list.
    ex = []                                         # Excluded list.
    for x in dr_9:                                  # For every line in direct repeats 9 list.
        split = x[0].split("_")                     # Split title based on underscore, call list split.
        if int(split[2]) < fnn:                     # If start value is smaller than the flank nucleotide number.
            bq = int(split[2])                      # Set query break point to start.
        else:                                       # If start value is bigger than the flank nucleotide number.
            bq = fnn                                # Set query break point to flank nucleotide number.
        bs = int(x[8]) - bq                         # Set subject breakpoint to length - flank nucleotide number.
        alq = bq - int(x[6])                        # Find difference of breakpoint - start of alignment for query.
        als = bs - int(x[9])                        # Find difference of breakpoint - start of alignment for subject.
        b_diff = alq-als                            # Find difference between breakpoint alignment.
        if abs >= b_diff >= -abs:                   # If breakpoint difference is >= -abs and <= abs.
            temp = []                               # Temporary list.
            temp.append('HR.DR')                    # Append category to temporary list.
            for value in x:                         # For every value in the line.
                temp.append(value)                  # Append it to the temporary list.
            pb.append(temp)                         # Append temporary list containing line to path blast+ list.
            y = 'HR.DR' + x[0]                      # Join pathway to string.
            u_titles.append(y)                      # Append line to unfiltered titles list.
        else:                                       # If alignment is shifted.
            temp = []                               # Temporary list.
            temp.append('HR.SDR')                   # Append category to temporary list.
            for value in x:                         # For every value in the line.
                temp.append(value)                  # Append it to the temporary list.
            ex.append(temp)                         # Append temporary list containing line to excluded list.
# Step : Circles that has SDR, but not DR, are moved from excluded list to path blast+ list.
    pbt = []                                        # path blast+ title list.
    for line in pb:                                 # For every in line in path blast+.
        pbt.append(line[1])                         # Append title to path blast+ title list.
    for line in ex:                                 # For every line in excluded list.
        if line[1] not in pbt:                          # If line title is in the list containing direct repeats.
            pb.append(line)                         # Append line for shifted direct repeat in path blast+ list.
            x = line[0] + line[1]                   # Join the first value with the second value.
            u_titles.append(x)                      # Append line to unfiltered titles list.
# Step : HR circles are written in path blast file, and excluded are written in excluded file.
    nm_pb = no_multiples(pb, 1)                     # Remove repeating elements in pb list based on column 2.
    pbw = []                                        # Path blast+ write list.
    for line in nm_pb:                              # For line in no multiples path blast+ list.
        for value in line:                          # For every value in the line.
            pbw.append(value)                       # Append value to pbw list.
            pbw.append('\t')                        # Append tab to pbw list.
        pbw.append('\n')                            # Append newline to pbw list.
    pbw2 = "".join(pbw)                             # Join all strings together and call it pbw2.
    path_blast.write(pbw2)                          # Write pbw2 into the path blast+ file.
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
    for line in new7:                               # For every line in list.
        titles.append(line[0])                      # Append first value to titles.
    for line in clean9:                             # For every line in list.
        titles.append(line[0])                      # Append first value to titles.
    new = [x for x in t_nhej if x not in titles]    # For every x in t_nhej, if it is not in titles, append to new.
    for x in new:                                   # For every value in new.
        path_blast.write('NHEJ' + "\t" + str(x) + "\n")  # Write line and newline to repair pathway file.
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
    c = sorted(co)                                  # Sort coordinates.
    coor = []                                       # Coordinates list for combined circles found in genome.
    file2 = open(circle_coor, "r")                  # Open file in read mode.
    s2 = file2.readlines()                          # Read lines in file and save them under variable.
    file2.close()                                   # Close file.
    for line in s2:                                 # For every line in list of strings.
        temp1 = line.strip('\n')                    # Remove trailing whitespace.
        temp2 = temp1.split('\t')                   # Split strings by tabs.
        coor.append(temp2)                          # Append line to coor list.
    coor2 = sorted(coor)                            # Sort coordinates.
    sco = []                                        # Make new list for sorted coordinates.
    for lc in coor2:                                # For every line in coordinates.
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
    pathway_count = {}                              # Make dictionary for counting circles in each pathway.
    for line in all_titles:                         # For every line in all titles list.
        if line[0] not in pathway_count:            # If title is not in pathway count dict.
            pathway_count[line[0]] = 1              # Add title as key to pathway count dict with value 1.
        else:                                       # If title is in the pathway count dict.
            count = pathway_count[line[0]]          # Find value for the title, call it count.
            count += 1                              # Take value for count and add 1.
            pathway_count[line[0]] = count          # Add updated value to title in pathway count dict.
    o = len(all_titles)                             # Count the number of circles in all titles list.
    pathway_count['Circles with pathway'] = o # Add circles with pathway as key, value pair in pathway count.
    tnc = len(c)                                    # Total number of circles.
    pathway_count['Number of circles'] = tnc        # Add number of circles as key, value pair in pathway count.
    pathway_count['Circles without pathway'] = tnc-o  # Add circles - pathway as key, value pair in pathway count.
    if 'MMEJ' in pathway_count:                     # If pathway is in pathway count dict.
        num = pathway_count['MMEJ']                 # Take value for pathway, call it num.
        percent = round(((num/o)*100), 2)           # Calculate percent of circles for pathway, rounded to 2 decimals.
        pathway_count['MMEJ%'] = percent            # Add percent for pathway as key, value pair to dict.
    if 'HR.DR' and 'HR.SDR' in pathway_count:       # If both pathways are in pathway count dict.
        num1 = pathway_count['HR.DR']               # Take value for pathway, call it num1.
        num2 = pathway_count['HR.SDR']              # Take value for pathway, call it num2.
        percent = round((((num1+num2)/o)*100), 2)   # Calculate percent of circles for pathway, rounded to 2 decimals.
        pathway_count['HR%'] = percent              # Add percent for pathway as key, value pair to dict.
    if 'HR.DR' in pathway_count:                    # If pathway is in pathway count dict.
        num = pathway_count['HR.DR']                # Take value for pathway, call it num.
        percent = round(((num/o)*100), 2)           # Calculate percent of circles for pathway, rounded to 2 decimals.
        pathway_count['HR.DR%'] = percent           # Add percent for pathway as key, value pair to dict.
    if 'NHEJ' in pathway_count:                     # If pathway is in pathway count dict.
        num = pathway_count['NHEJ']                 # Take value for pathway, call it num.
        percent = round(((num/o)*100), 2)           # Calculate percent of circles for pathway, rounded to 2 decimals.
        pathway_count['NHEJ%'] = percent            # Add percent for pathway as key, value pair to dict.
# Step 9: Percentages of inverted repeats, and overlap between other pathways and inverted repeats.
    if 'HR.SDR' in pathway_count:                   # If pathway is in pathway count dict.
        num = pathway_count['HR.SDR']               # Take value for pathway, call it num.
        percent = round(((num/o)*100), 2)           # Calculate percent of circles for pathway, rounded to 2 decimals.
        pathway_count['HR.SDR%'] = percent          # Add percent for pathway as key, value pair to dict.
    nm_inv = no_multiples(inverted, 2)              # Remove repeating values from inverted list.
    s_inv = []                                      # Make list for split inverted.
    for x in nm_inv:                                # For every line in no multiples list.
        split = x[2].split("_")                     # Split string into chrom, start, end.
        s_inv.append(split)                         # Append titles to split inverted list.
    overlap = []                                    # Overlap between inverted and other pathways list.
    overlap_index = []                              # No overlap between inverted and pathway list.
    for x in s_inv:                                 # For every line in split inverted list.
        for y in all_titles:                        # For every line in all titles list.
            if x[1] == y[1] and x[2] == y[2] and x[3] == y[3]:  # If chrom, start and end, match.
                temp = []                           # Temporary list.
                temp.append('Inverted repeat')      # Append inverted repeat title to temp list.
                temp.append(y[0])                   # Append first value in all titles list to temp list.
                for value in x[1:]:                 # For every value in split inverted line, starting with value 2.
                    temp.append(value)              # Append value to temp list.
                overlap.append(temp)                # Append temp to overlap.
                index = s_inv.index(x)              # Take index of x in split inverted.
                overlap_index.append(index)         # Append index to overlap index list.
    nm_index = []                                   # No multiples index list.
    for value in overlap_index:                     # For every value in overlap index list.
        if value not in nm_index:                   # If value is not in the no multiples index list.
            nm_index.append(value)                  # Append value to no multiples index list.
    no_overlap = []                                 # No overlap between pathway and inverted.
    for i in s_inv:                                 # For line in split inverted list.
        temp = []                                   # Temporary list.
        if s_inv.index(i) not in overlap_index:     # If s inv index is not in overlap index.
            temp.append("Inverted repeat")          # Append inverted repeat to temp list.
            for value in i:                         # For every value in index.
                temp.append(value)                  # Append value to temp list.
            no_overlap.append(temp)                 # Append temp to no overlap.
    wex = []                                        # Write excluded list.
    for line in overlap:                            # For every line in overlap list.
        for value in line:                          # For every value in the line.
            wex.append(value)                       # Append value to wex list.
            wex.append('\t')                        # Append tab to wex list.
        wex.append('\n')                            # Append newline to wex list.
    for line in no_overlap:                         # For every line in no overlap list.
        for value in line:                          # For every value in the line.
            wex.append(value)                       # Append value to wex list.
            wex.append('\t')                        # Append tab to wex list.
        wex.append('\n')                            # Append newline to wex list.
    jex = "".join(wex)                              # Join list of strings, call it jex.
    excluded.write(jex)                             # Write jex to excluded file.
    lo = len(overlap)                               # Number of inverted repeat overlap with other pathways.
    lno = len(no_overlap)                           # Number of inverted repeats with a pathway.
    to_tno = lo + lno                               # Number of inverted with pathway and number of inverted without.
    pathway_count['IR'] = to_tno                    # Add IR and number of IR as key, value pair.
    pathway_count['IR, pathway'] = lo               # Add IR pathway and number of IR as key, value pair.
    pathway_count['IR, no pathway'] = lno           # Add IR no pathway and number of IR as key, value pair.
    if 'IR' in pathway_count:                       # If IR is in pathway count dict.
        num = pathway_count['IR']                   # Take value for IR, call it num.
        percent = round(((num / o) * 100), 2)       # Calculate percent of circles for IR, rounded to 2 decimals.
        pathway_count['IR%'] = percent              # Add percent for IR as key, value pair to dict.
    if 'IR, pathway' in pathway_count:              # If IR pathway is in pathway count dict.
        num = pathway_count['IR, pathway']          # Take value for IR pathway, call it num.
        percent = round(((num / o) * 100), 2)       # Calculate % of circles for IR pathway, rounded to 2 decimals.
        pathway_count['IR, pathway%'] = percent     # Add % for IR pathway as key, value pair to dict.
    if 'IR, no pathway' in pathway_count:           # If IR no pathway is in pathway count dict.
        num = pathway_count['IR, no pathway']       # Take value for pathway, call it num.
        percent = round(((num / o) * 100), 2)       # Calculate % of circles for IR no pathway, rounded to 2 decimals.
        pathway_count['IR, no pathway%'] = percent  # Add % for pathway as key, value pair to dict.
    path_blast.close()                                # Close file.
    path_coor.close()                               # Close file.
    excluded.close()                                # Close file.
    return

