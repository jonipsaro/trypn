#Jonathan Ipsaro
#Last reviewed: January 12, 2019
#
#This script find overlapping proteins (by protein accession number)
#between two input files.
#
#Input files should be CSV output from MASCOT search results with protein
#accession numbers in column 4 (python index 3)
#
#One file will be output:
#1. The summary table with the number of proteins from each file and
#number that overlap

#The user should adjust the following four variables to run this script
file1 = open("CSV/Trypsin_proteins.csv", "r")
file2 = open("CSV/TrypN_proteins.csv", "r")
file1_protein = "Trypsin"			#These are names to include in headers of output
file2_protein = "Tryp-N"			#Make sure these correspond to the correct file above

from datetime import datetime

#Function to clean up MASCOT CSV format
def CleanLine(line):
	quote = False

	for i in range(len(line)):
		if quote == False and line[i] == "\"": quote = True
		elif quote == True and line[i] == "\"": quote = False

		if ((quote == True) and (line[i] == ",")): line = line[:i]+":"+line[i+1:]

	return line


proteins = {}   #dictionary that will contain proteins as keys
				#and a list [file1, file2] as entries

file1_uniq_proteins = set()
file2_uniq_proteins = set()


#I. Build a set from the first dataset
build = False
for line in file1:      								#read in each line
	if line.find("Family") == 0: 						#ignore header lines
		build = True
		line = file1.readline()

	if build == False: continue

	protein = CleanLine(line.strip()).split(",")[3]		#extract the protein accession number
	file1_uniq_proteins.add(protein)					#add the protein to the file1 set
file1.close()


#II. Build a set from the second dataset
build = False
for line in file2:      #read in each peptide and make a unique listing
	if line.find("Family") == 0: 						#read in each line
		build = True									#ignore header lines
		line = file2.readline()

	if build == False: continue

	protein = CleanLine(line.strip()).split(",")[3]		#extract the protein accession number
	file2_uniq_proteins.add(protein)					#add the protein to the file1 set
file2.close()


#III. Build a dictionary for all proteins with proteins as keys and [file1, file2] as entries
proteins = {}

for each in file1_uniq_proteins:
	proteins[each] = [file1_protein, None]

for each in file2_uniq_proteins:
	if each not in proteins:
		proteins[each] = [None, file2_protein]
	else:
		proteins[each][1] = file2_protein


#IV. Go through the dictionary and count each Venn diagram class
file1_only = file2_only = both = other = 0

for protein in proteins:
	if (proteins[protein][0:2] == [file1_protein,file2_protein]): both += 1
	elif (proteins[protein][0:2] == [file1_protein, None]): file1_only += 1
	elif (proteins[protein][0:2] == [None,file2_protein]): file2_only += 1
	else: other += 1


#IV. File Output
fout = open("protein_comparison.log", "w+")
fout.write("protein_id\t"+file1_protein+"\t"+file2_protein+"\n")
sorted_proteins = sorted(proteins.keys())

for protein in sorted_proteins:   
	if proteins[protein][0] == file1_protein: file1_out = "1"
	else: file1_out = "0"

	if proteins[protein][1] == file2_protein: file2_out = "1"
	else: file2_out = "0"

	fout.write(protein+"\t"+file1_out+"\t"+file2_out+"\n")
fout.close()


#V. Log Output
logfile = open("protein_comparison_summary.log", "w+")
print("compare_proteins_2.py", file = logfile)
print("Run on", str(datetime.now()), file = logfile)

print("================================", file = logfile)
print("INPUT FILE SUMMARY", file = logfile)
print("Number of "+file1_protein+":",len(file1_uniq_proteins), file = logfile)
print("Number of "+file2_protein+":",len(file2_uniq_proteins), file = logfile)

print("================================", file = logfile)
print("OUTPUT STATISTICS", file = logfile)
print("Number of proteins unique to the "+file1_protein+" digest:",file1_only, file = logfile)
print("Number of proteins unique to the "+file2_protein+" digest:",file2_only, file = logfile)
print("Number of overlapping proteins:",both, file = logfile)

print("================================", file = logfile)
print("SANITY CHECK", file = logfile)
if len(file1_uniq_proteins) == file1_only + both: print("All "+file1_protein+" proteins accounted for. PASS.", file = logfile)
else: print("Some "+file1_protein+" proteins not accounted for. FAIL.", file = logfile)
if len(file2_uniq_proteins) == file2_only + both: print("All "+file2_protein+" proteins accounted for. PASS.", file = logfile)
else: print("Some "+file2_protein+" proteins not accounted for. FAIL.", file = logfile)
if other != 0: print("Some proteins were not properly filtered. Please check the code and input.", file = logfile)

logfile.close()