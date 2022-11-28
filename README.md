# Homology-identification-pipeline
A method for processing output from the Circle-Map method, to identify homology and produce a suggested pathway of origin for circular DNA.

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

HRI is run once, blast is run twice with a shorter and longer word length, BP is run once.
