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


def average(list):
    '''Finds smallest value of index 4 and 5, and the smallest value of index 6 and 7, then calculates the average
    of all, and individually '''
    l = []                                          # Make list for total values.
    l1 = []                                         # Make list for start values.
    l2 = []                                         # Make list for end values.
    for line in list:                               # For every line in list.
        if line[4] < line[5]:                       # If qstart is smaller than qend.
            l.append(line[4])                       # Append qstart to total list.
            l1.append(line[4])                      # Append qstart to start list.
        elif line[4] > line[5]:                     # If qstart is bigger than qend.
            l.append(line[5])                       # Append qend to total list.
            l1.append(line[5])                      # Append qend to start list.
        if line[6] < line[7]:                       # If sstart is smaller than send.
            l.append(line[6])                       # Append sstart to total list.
            l2.append(line[6])                      # Append sstart to end list.
        elif line[6] > line[7]:                     # If sstart is bigger than send.
            l.append(line[7])                       # Append send to total list.
            l2.append(line[7])                      # Append send to end list.
    a = len(l)                                      # Take length of total list.
    b = len(l1)                                     # Take length of start list.
    c = len(l2)                                     # Take length of end list.
    if a > 0:                                       # If a length is larger than 0.
        la = round(sum(l)/a)                        # Take average of a list, and round to whole number.
    else:                                           # If a length is 0.
        la = "NA"                                   # Call la Not Available.
    if b > 0:                                       # If b length is larger than 0.
        l1b = round(sum(l1)/b)                      # Take average of b list, and round to whole number.
    else:                                           # If b length is 0.
        l1b = "NA"                                  # Call l1c Not Available.
    if c > 0:                                       # If c length is larger than 0.
        l2c = round(sum(l2)/c)                      # Take average of c list, and round to whole number.
    else:                                           # If c length is 0.
        l2c = "NA"                                  # Call l2c Not Available.
    return ["{} circles".format(round(a/2)), "Total: {}".format(la), "Start: {}".format(l1b), "End: {}".format(l2c)]


def pathway_sort(align7, align9, query, remaining, circle_coor, flank_nucleotide_number, allowed_breakpoint_shift,
                 gene_column):
    '''Sort results into HR, MMEJ and NHEJ.'''
    abs = allowed_breakpoint_shift                  # Rename allowed breakpoint shift.
    path_blast = open("Pathways_blast+.txt", "a")   # Open file, (will create file, if not found).
    path_coor = open("Pathway_coor.txt", "a")       # Open file for pathway coordinates for circles.
    stat = open("Statistics.txt", "a")              # Open file for statistics.
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
        else:                                       # If start > end for subject.
            inv = []                                # Make temp inverted list.
            inv.append("Inverted repeat")       # Append category to temp inverted list.
            inv.append("Align7")                # Append source to temp inverted list.
            for value in x:                     # For every value in x.
                inv.append(value)               # Append value to inverted temp list.
            inverted.append(inv)                # Append inverted temp list to inverted list.
    dr_9 = []                                       # Direct repeat from align9.
    for x in clean9:                                # For every line in clean9.
        if int(x[6]) <= int(x[7]) and int(x[9]) <= int(x[10]):  # If start <= end for query and subject.
            dr_9.append(x)                          # Append to hr direct repeat list.
        else:                                       # If start > end alignment for subject.
            inv = []                                # Make temp inverted list.
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
        bs = int(x[8]) - fnn                         # Set subject breakpoint to length - flank nucleotide number.
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
# Step 5: Circles that has SDR, but not DR, are moved from excluded list to path blast+ list.
    pbt = []                                        # path blast+ title list.
    for line in pb:                                 # For every in line in path blast+.
        pbt.append(line[1])                         # Append title to path blast+ title list.
    for line in ex:                                 # For every line in excluded list.
        if line[1] not in pbt:                          # If line title is in the list containing direct repeats.
            pb.append(line)                         # Append line for shifted direct repeat in path blast+ list.
            x = line[0] + line[1]                   # Join the first value with the second value.
            u_titles.append(x)                      # Append line to unfiltered titles list.
# Step 6: HR circles are written in path blast file, and excluded are written in excluded file.
    nm_pb = no_multiples(pb, 1)                     # Remove repeating elements in pb list based on column 2.
    pbw = []                                        # Path blast+ write list.
    for line in nm_pb:                              # For line in no multiples path blast+ list.
        for value in line:                          # For every value in the line.
            pbw.append(value)                       # Append value to pbw list.
            pbw.append('\t')                        # Append tab to pbw list.
        pbw.append('\n')                            # Append newline to pbw list.
    pbw2 = "".join(pbw)                             # Join all strings together and call it pbw2.
    path_blast.write(pbw2)                          # Write pbw2 into the path blast+ file.
# Step 7: Making the NHEJ file.
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
# Step 8: Removing repeating values in unfiltered list and formatting to all titles list.
    f_titles = []                                   # New list for titles without multiples.
    for element in u_titles:                        # For every element in unsorted titles.
        if element not in f_titles:                 # If the element is not already in the list.
            f_titles.append(element)                # Append element to filtered titles.
# Step 9: Classifying circles by repair pathway.
    all_titles = []                                 # Titles for all circles except Excluded and Remaining.
    for x in f_titles:                              # For every line in filtered titles list.
        split = x.split("_")                        # Split string into chrom, start, end.
        all_titles.append(split)                    # Append titles to all titles list.
    file2 = open(circle_coor, "r")                  # Open file in read mode.
    s1 = file2.readlines()                          # Read lines in file and save them under variable.
    file2.close()                                   # Close file.
    coor = []                                       # Coordinates list for merged circles found in genome.
    for line in s1:                                 # For every line in list of strings.
        temp1 = line.strip('\n')                    # Remove trailing whitespace.
        temp2 = temp1.split('\t')                   # Split strings by tabs.
        coor.append(temp2)                          # Append line to coor list.
    coor1 = sorted(coor)                            # Sort coordinates.
    sco = []                                        # Make new list for sorted coordinates.
    for lc in coor1:                                # For every line in coordinates.
        for lt in all_titles:                       # For every line in all titles.
            if lc[0] == lt[1] and lc[1] == lt[2] and lc[2] == lt[3]:  # If lc index 0, 1 & 2 match lt index 1, 2 & 3.
                temp = []                           # Make temporary list.
                temp.append(lt[0])                  # Append pathway for circle to temp list.
                for value in lc:                    # For every value in lc.
                    temp.append(value)              # Append value to temp list.
                sco.append(temp)                    # Append temp to sco list.
    ss = sorted(sco)                                # Sort lines by pathway.
    ss1 = []                                        # Make new empty list for the formatted sorted sco.
    for line in ss:                                 # For every line in sorted sco.
        for value in line:                          # For every value in the line.
            ss1.append(value)                       # Append value to ss1.
            ss1.append("\t")                        # Append tab to ss1.
        ss1.append("\n")                            # Append newline at the end of the line to ss1.
    s = "".join(ss1)                                # Join all strings into one string and rename list to s.
    path_coor.write(s)                              # Write s list to path coor file.
# Step 10: Counting number of circles belonging to each pathway as well as percentages.
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
    file = open(remaining, 'r')                   # Open file in read mode.
    s2 = file.readlines()                           # Read lines in file and save them under variable.
    file.close()                                    # Close file.
    tnc = len(s2) + o                               # Total number of circles.
    pathway_count['Number of circles'] = tnc        # Add number of circles as key, value pair in pathway count.
    pathway_count['Circles without pathway'] = tnc-o  # Add circles - pathway as key, value pair in pathway count.
    if 'MMEJ' in pathway_count:                     # If pathway is in pathway count dict.
        num = pathway_count['MMEJ']                 # Take value for pathway, call it num.
        percent = round(((num/o)*100), 2)           # Calculate percent of circles for pathway, rounded to 2 decimals.
        pathway_count['MMEJ%'] = percent            # Add percent for pathway as key, value pair to dict.
    if 'HR.DR' in pathway_count:                    # If both pathways are in pathway count dict.
        if 'HR.SDR' in pathway_count:
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
    if 'HR.SDR' in pathway_count:                   # If pathway is in pathway count dict.
        num = pathway_count['HR.SDR']               # Take value for pathway, call it num.
        percent = round(((num/o)*100), 2)           # Calculate percent of circles for pathway, rounded to 2 decimals.
        pathway_count['HR.SDR%'] = percent          # Add percent for pathway as key, value pair to dict.
# Step 11: Percentages of inverted repeats, and overlap between other pathways and inverted repeats.
    nm_inv = no_multiples(inverted, 2)              # Remove repeating values from inverted list.
    s_inv = []                                      # Make list for split inverted.
    for x in nm_inv:                                # For every line in no multiples list.
        split = x[2].split("_")                     # Split string into chrom, start, end.
        s_inv.append(split)                         # Append titles to split inverted list.
    overlap = []                                    # Overlap between inverted and other pathways list.
    overlap_index = []                              # No overlap between inverted and pathway list.
    for x in s_inv:                                 # For every line in split inverted list.
        for y in all_titles:                        # For every line in all titles list.
            if x[1] == y[1] and x[2] == y[2] and x[3] == y[3]: # If chrom, start and end, match.
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
# Step 12: List genes in circles by pathway.
    gc = gene_column - 1                            # Take gene column minus 1 to get python index for column.
    hrdr = []                                       # Make list to contain genes found in direct repeats.
    hrsdr = []                                      # Make list to contain genes found in shifted direct repeats.
    mmej = []                                       # Make list to contain genes found in mmej.
    for l in sco:                                   # For every line in sco.
        if l[0] == "HR.DR":                         # If the first value is HR.DR.
            t1 = l[gc].split(",")                   # Split gene column and insert comma and space, call it t1.
            for value in t1:                        # For every gene in t1.
                hrdr.append(value)                  # Append gene to hrdr list.
        elif l[0] == "HR.SDR":                      # If the first value is HR.SDR.
            t1 = l[gc].split(",")                   # Split gene column and insert comma and space, call it t1.
            for value in t1:                        # For every gene in t1.
                hrsdr.append(value)                 # Append gene to hrsdr list.
        elif l[0] == "MMEJ":                        # If the first value is MMEJ.
            t1 = l[gc].split(",")                   # Split gene column and insert comma and space, call it t1.
            for value in t1:                        # For every gene in t1.
                mmej.append(value)                  # Append gene to mmej list.
    d = []                                          # Create list for hrdr values, without repeating values.
    for value in hrdr:                              # For every value in hrdr list.
        if value not in d:                          # If value is not in d.
            d.append(value)                         # Append value to d.
    s = []                                          # Create list for hrsdr values, without repeating values.
    for value in hrsdr:                             # For every value in hrsdr list.
        if value not in s:                          # If value is not in s.
            s.append(value)                         # Append value to s.
    m = []                                          # Create list for mmej values, without repeating values.
    for value in mmej:                              # For every value in mmej list.
        if value not in m:                          # If value is not in m.
            m.append(value)                         # Append value to s.
    jd = ", ".join(d)                               # Join all values in d list into one string.
    js = ", ".join(s)                               # Join all values in s list into one string.
    jm = ", ".join(m)                               # Join all values in m list into one string.
    pathway_count['HR.DR_genes'] = jd               # Add HR.DR genes as key, value pair in pathway count dict.
    pathway_count['HR.SDR_genes'] = js              # Add HR.SDR genes as key, value pair in pathway count dict.
    pathway_count['MMEJ_genes'] = jm                # Add MMEJ genes as key, value pair in pathway count dict.
# Step 13: Average length of homologous alignments.
    hr = []                                         # Make list for HR average length of alignments.
    dr = []                                         # Make list for HR.DR average length of alignments.
    sdr = []                                        # Make list for HR.SDR average length of alignments.
    ir = []                                         # Make list for IR average length of alignments.
    for line in nm_pb:                              # For every line in nm_pb list.
        if line[0] == "HR.DR":                      # If first value in line is HR.DR.
            hr.append(int(line[4]))                 # Append alignment length to hr list.
            dr.append(int(line[4]))                 # Append alignment length to dr list.
        elif line[0] == "HR.SDR":                   # If first value in line is HR.SDR.
            hr.append(int(line[4]))                 # Append alignment length to hr list.
            sdr.append(int(line[4]))                # Append alignment length to sdr list.
    for line in nm_inv:                             # For every line in nm_inv.
        ir.append(int(line[5]))                     # Append alignment length to ir list.
    hr = round(sum(hr)/len(hr))                     # Take average of hr list, and round to whole number.
    if len(dr) > 0:                                 # If the length of dr is larger than 0.
        dr = round(sum(dr)/len(dr))                 # Take average of dr list, and round to whole number.
    else:                                           # If the length of dr is 0.
        dr = 0                                      # Replace empty list with 0.
    sdr = round(sum(sdr)/len(sdr))                  # Take average of sdr list, and round to whole number.
    ir = round(sum(ir)/len(ir))                     # Take average of ir list, and round to whole number.
    pathway_count['Av. len. align. HR'] = hr        # Add HR average length as key, value pair in pathway count dict.
    pathway_count['Av. len. align. HR.DR'] = dr     # Add HR.DR ave. len. as key, value pair in pathway count dict.
    pathway_count['Av. len. align. HR.SDR'] = sdr   # Add HR.SDR ave. len. as key, value pair in pathway count dict.
    pathway_count['Av. len. align. IR'] = ir        # Add IR ave. len. as key, value pair in pathway count dict.
# Step 14: Average length from breakpoint to homologous region.
    hl1 = []                                        # Create homology length 1 list.
    for x in nm_pb:                                 # For every line in nm_pb.
        temp = []                                   # Create temporary list.
        temp.append(x[0])                           # Append pathway to temp.
        temp.append(x[2])                           # Append title to temp.
        split = x[2].split("_")  # Split title based on underscore, call list split.
        if int(split[2]) < fnn:                     # If start value is smaller than the flank nucleotide number.
            bq = int(split[2])                      # Set query break point to start.
        else:                                       # If start value is bigger than the flank nucleotide number.
            bq = fnn                                # Set query break point to flank nucleotide number.
        bs = int(x[9]) - fnn                        # Set subject breakpoint to length - flank nucleotide number.
        qsd = bq - int(x[7])                        # query-start-difference; start breakpoint - qstart alignment
        qed = bq - int(x[8])                        # query-end-difference; start breakpoint - qend alignment
        ssd = bs - int(x[10])                       # subject-start-difference; end breakpoint - sstart alignment
        sed = bs - int(x[11])                       # subject-end-difference; end breakpoint - send alignment
        temp.append(qsd)                            # Append qsd to temp.
        temp.append(qed)                            # Append qed to temp.
        temp.append(ssd)                            # Append ssd to temp.
        temp.append(sed)                            # Append sed to temp.
        temp.append(x[4])                           # Append length of alignment to temp.
        temp.append(x[9])                           # Append length of circle to temp.
        hl1.append(temp)                            # Append temp to homology length 1 list.
# Step 15: Where does the homology reside, outside breakpoints, over, or inside.
    hl2 = []                                        # Create homology length 2 list.
    for line in hl1:                                # For every line in hl1.
        temp = []                                   # Create temporary list.
        split = line[1].split("_")                  # Split title based on underscore, call list split.
        spp = line[2] >= 0 and line[3] >= 0         # Query: Start align +, end align +.
        spn = line[2] >= 0 and line[3] < 0          # Query: Start align +, end align -.
        snp = line[2] < 0 and line[3] >= 0          # Query: Start align -, end align +.
        snn = line[2] < 0 and line[3] < 0           # Query: Start align -, end align -.
        epp = line[4] >= 0 and line[5] >= 0         # Subject: Start align +, end align +.
        epn = line[4] >= 0 and line[5] < 0          # Subject: Start align +, end align -.
        enp = line[4] < 0 and line[5] >= 0          # Subject: Start align -, end align +.
        enn = line[4] < 0 and line[5] < 0           # Subject: Start align -, end align -.
        if (spn or snp) and (epn or enp):           # If homology found across breakpoints.
            temp.append(line[0] + ".AA")                # Append AA to pathway, bc alignments are Across.
        elif (spn or snp) and (epp or enn):         # If homology found across/shifted compared to breakpoints.
            if (spn or snp) and epp:                    # If query is +- or -+ and subject is ++.
                temp.append(line[0] + ".AS.B")          # Append AS.B, bc alignments are Across/Shifted.Before.
            else:                                       # If query is +- or -+ and subject is --.
                temp.append(line[0] + ".AS.A")          # Append AS.A, bc alignments are Across/Shifted.After.
        elif (spp or snn) and (epn or enp):         # If homology found shifted compared to/across breakpoints.
            if spp and (epn or enp):                    # If query is ++ and subject is +- or -+.
                temp.append(line[0] + ".SA.B")          # Append SA.B, bc alignments are Shifted/Across.Before.
            else:                                       # If query is -- and subject is +- or -+.
                temp.append(line[0] + ".SA.A")          # Append SA.A, bc alignments are Shifted/Across.After.
        else:                                       # If homology found shifted compared to breakpoints.
            if spp and epp:                             # If query is ++ and subject is ++.
                temp.append(line[0] + ".SS.B")          # Append SS.B, bc alignments are Shifted/Shifted.Before.
            elif spp and enn:                           # If query is ++ and subject is --.
                temp.append(line[0] + ".SS.O")          # Append SS.O, bc alignments are Shifted/Shifted.Outside.
            elif snn and epp:                           # If query is -- and subject is ++.
                temp.append(line[0] + ".SS.I")          # Append SS.I, bc alignments are Shifted/Shifted.Inside.
            else:                                       # If query is -- and subject is --.
                temp.append(line[0] + ".SS.A")          # Append SS.A, bc alignments are Shifted/Shifted.After.
        for value in split[1:]:                     # For every value in split, starting with the second.
            temp.append(value)                      # Append value to temp.
        for value in line[2:]:                      # For every value in line, starting with the third.
            temp.append(value)                      # Append value to temp.
        hl2.append(temp)                            # Append temp list to hl2.
    for line in hl2:                                # For every line in hl2.
        if line[4] < 0:                             # If the number in index 4 is negative.
            line[4] = -line[4]                      # Make the number positive.
        if line[5] < 0:                             # If the number in index 5 is negative.
            line[5] = -line[5]                      # Make the number positive.
        if line[6] < 0:                             # If the number in index 6 is negative.
            line[6] = -line[6]                      # Make the number positive.
        if line[7] < 0:                             # If the number in index 7 is negative.
            line[7] = -line[7]                      # Make the number positive.
    daa = []                                        # Make list for DR.AA homology.
    dsa = []                                        # Make list for DR.SA homology.
    das = []                                        # Make list for DR.AS homology.
    dss = []                                        # Make list for DR.SS homology.
    saa = []                                        # Make list for SDR.AA homology.
    ssa = []                                        # Make list for SDR.SA homology.
    sas = []                                        # Make list for SDR.AS homology.
    sss = []                                        # Make list for SDR.SS homology.
    for line in hl2:                                # For every line in hl2.
        split = line[0].split(".")                  # Split pathway based on punctuation, call list split.
        if split[1] == "DR" and split[2] == "AA":   # If circle is DR.AA.
            daa.append(line)                        # Append to corresponding daa list.
        elif split[1] == "DR" and split[2] == "SA":  # If circle is DR.SA.
            dsa.append(line)                        # Append to corresponding dsa list.
        elif split[1] == "DR" and split[2] == "AS":  # If circle is DR.AS.
            das.append(line)                        # Append to corresponding das list.
        elif split[1] == "DR" and split[2] == "SS":  # If circle is DR.SS.
            dss.append(line)                        # Append to corresponding dss list.
        elif split[1] == "SDR" and split[2] == "AA":  # If circle is SDR.AA.
            saa.append(line)                        # Append to corresponding saa list.
        elif split[1] == "SDR" and split[2] == "SA":  # If circle is SDR.SA.
            ssa.append(line)                        # Append to corresponding ssa list.
        elif split[1] == "SDR" and split[2] == "AS":  # If circle is SDR.AS.
            sas.append(line)                        # Append to corresponding sas list.
        elif split[1] == "SDR" and split[2] == "SS":  # If circle is SDR.SS.
            sss.append(line)                        # Append to corresponding sss list.
    daa_val = average(daa)                          # Call average function, get list with values.
    pathway_count['Av. len. to homology, DR.AA'] = daa_val  # Add daa values to pathway count dict.
    dsa_val = average(dsa)                          # Call average function, get list with values.
    pathway_count['Av. len. to homology, DR.SA'] = dsa_val  # Add dsa values to pathway count dict.
    das_val = average(das)                          # Call average function, get list with values.
    pathway_count['Av. len. to homology, DR.AS'] = das_val  # Add das values to pathway count dict.
    dss_val = average(dss)                          # Call average function, get list with values.
    pathway_count['Av. len. to homology, DR.SS'] = dss_val  # Add dss values to pathway count dict.
    saa_val = average(saa)                          # Call average function, get list with values.
    pathway_count['Av. len. to homology, SDR.AA'] = saa_val  # Add saa values to pathway count dict.
    ssa_val = average(ssa)                          # Call average function, get list with values.
    pathway_count['Av. len. to homology, SDR.SA'] = ssa_val  # Add ssa values to pathway count dict.
    sas_val = average(sas)                          # Call average function, get list with values.
    pathway_count['Av. len. to homology, SDR.AS'] = sas_val  # Add sas values to pathway count dict.
    sss_val = average(sss)                          # Call average function, get list with values.
    pathway_count['Av. len. to homology, SDR.SS'] = sss_val  # Add sss values to pathway count dict.
# Step 16: Write pathway count dictionary to statistics file.
    for key, value in pathway_count.items():        # For every key, value pair in pathway count dict.
        stat.write("%s: %s\n" % (key, value))       # Write to stat file as "key: value" newline.
# Step 17: Closing output files.
    path_blast.close()                              # Close file.
    path_coor.close()                               # Close file.
    stat.close()                                    # Close file.
    excluded.close()                                # Close file.
    return

