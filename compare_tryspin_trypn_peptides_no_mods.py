#Jonathan Ipsaro
#Last reviewed: January 12, 2019
#
#This script finds overlapping peptides between
#Trypsin- and TrypN-digested samples.
#
#Trypic peptides should end in K/R.
#TrypN-derived peptides should begin in K/R.
#
#In this version of the script, peptide uniqueness is based only on sequence.
#Modifications are NOT considered.
#
#Input files should be MASCOT CSV files with peptides in the 24th column (python index 23).
#
#Note: This script assumes that all input peptides are validated as coming from the appropriate enzymes.
#(i.e. tryptic peptides should end in K/R OR the C-terminus of the protein).  This is the appropriate
#treatment of data from MS searches with protease parameters included.
#
#Note: There may be a small number of peptides that cannont be distinguished. For example,
#if KLEELELDEQQR and KLEELELDEQQK are both present in the same trypsin input file,
#after trimming, these will result in the same sequence (KLEELELDEQQ). Since we expect these to be
#rare, this should not significantly affect the analysis, but should be noted. This script
#will only count one of those instances.
#
#Two files will be output:
#1. The summary table with the number of peptides from each file and number that overlap
#2. An output file (tab-delimited) that lists all valid peptides, the matching Trypsin peptide(s), and the mathing TrypN peptide(s).

from datetime import datetime

#Function to clean up MASCOT CSV format
def CleanLine(line):
	quote = False

	for i in range(len(line)):
		if quote == False and line[i] == "\"": quote = True
		elif quote == True and line[i] == "\"": quote = False

		if ((quote == True) and (line[i] == ",")): line = line[:i]+":"+line[i+1:]

	return line

peptides = {}   #dictionary that will contain peptides as keys
				#and a list [trypsin, trypn] as entries

trypc = set()           #counters for peptides
trypn = set()
trypc_KR = set()        #counter for indistinguishable peptides (i.e. those found that end in both K AND R)
trypn_KR = set()

#I. Build sets for the Tryptic peptides.
#First, the unique peptide list will be generated.
#Then the list will be checked for peptides that cannot be distinguished by the other protease.

trypc_file = open("CSV/Trypsin.csv", "r")
trypc_uniq_peptides = set()

build = False
for line in trypc_file:      #read in each peptide and make a unique set
	if line.find("prot_hit") == 0: 
		build = True
		line = trypc_file.readline()

	if build == False: continue

	peptide = CleanLine(line.strip()).split(",")[23]

	trypc_uniq_peptides.add(peptide)

for peptide in trypc_uniq_peptides:			#Check for indistinguishable peptides
	if (peptide[-1] == "K") and (peptide[:-1]+"R" in trypc_uniq_peptides): trypc_KR.add(peptide)
	elif (peptide[-1] == "R") and (peptide[:-1]+"K" in trypc_uniq_peptides): trypc_KR.add(peptide)

	if (peptide[len(peptide)-1] == "K") or (peptide[len(peptide)-1] == "R"):    #trim C-terminal K/R if needed
		trimmed = peptide[:-1]
	else:
		trimmed = peptide
	
	if trimmed not in peptides:
		peptides[trimmed] = ["Trypsin", None, [], []]

	peptides[trimmed][2].append(peptide)

trypc_file.close()


#II. Build sets for the Tryp-N peptides
#First, the unique peptide list will be generated.
#Then the list will be checked for peptides that cannot be distinguished by the other protease.

trypn_file = open("CSV/TrypN.csv", "r")
trypn_uniq_peptides = set()

build = False
for line in trypn_file:      #read in each peptide and make a unique set
	if line.find("prot_hit") == 0: 
		build = True
		line = trypn_file.readline()

	if build == False: continue

	peptide = CleanLine(line.strip()).split(",")[23]

	trypn_uniq_peptides.add(peptide)

for peptide in trypn_uniq_peptides:			#Check for indistinguishable peptides
	if (peptide[0] == "K") and ("R"+peptide[1:] in trypn_uniq_peptides): trypn_KR.add(peptide)
	elif (peptide[0] == "R") and ("K"+peptide[1:] in trypn_uniq_peptides): trypn_KR.add(peptide)

	if (peptide[0] == "K") or (peptide[0] == "R"):    #trim N-terminal K/R if needed
		trimmed = peptide[1:]
	else:
		trimmed = peptide

	if trimmed not in peptides:
		peptides[trimmed] = [None, "TrypN", [], []]
	else:
		peptides[trimmed][1] = "TrypN"
	
	peptides[trimmed][3].append(peptide)
trypn_file.close()


#III. Go through the dictionary and count each Venn diagram class
trypc_only = trypn_only = both = other = 0
both2 = 0
for peptide in peptides:
    if (peptides[peptide][0:2] == ["Trypsin","TrypN"]): both += 1
    elif (peptides[peptide][0:2] == ["Trypsin", None]): trypc_only += 1
    elif (peptides[peptide][0:2] == [None,"TrypN"]): trypn_only += 1
    else: other += 1


#IV. Full Comparison Output
fout = open("peptide_comparison_no_mods.log", "w+")
sorted_peptides = sorted(peptides.keys())

for peptide in sorted_peptides:
    if peptides[peptide][0] == None: tryp_seq = "[None]"
    else: tryp_seq = str(peptides[peptide][2])
	
    if peptides[peptide][1] == None: trypn_seq = "[None]"
    else: trypn_seq = str(peptides[peptide][3])
	
    fout.write(peptide+"\tTrypsin\t"+tryp_seq+"\tTrypN\t"+trypn_seq+"\n")
fout.close()

#V. Summary Output
logfile = open("peptide_comparison_summary_no_mods.log", "w+")
print("compare_trypsin_trypn_peptides_no_mods.py", file = logfile)
print("Run on", str(datetime.now()), file = logfile)

print("================================", file = logfile)
print("INPUT FILE SUMMARY", file = logfile)
print("Number of Tryptic peptides:",len(trypc_uniq_peptides), file = logfile)
print("Number of TrypN peptides:",len(trypn_uniq_peptides), file = logfile)

print("================================", file = logfile)
print("OUTPUT STATISTICS", file = logfile)
print("Number of peptides unique to the Trypsin digest:",trypc_only, file = logfile)
print("Number of peptides unique to the TrypN digest:",trypn_only, file = logfile)
print("Number of overlapping peptides:",both, file = logfile)
print("Number of undistinguishable peptides in Trypsin and TrypN files:",len(trypn_KR)/2,len(trypc_KR)/2, file = logfile)

print("================================", file = logfile)
print("SANITY CHECK", file = logfile)
if (len(trypc_uniq_peptides) - len(trypc_KR)/2) == trypc_only + both: print("All Trypsin peptides accounted for. PASS.", file = logfile)
else: print("Some Trypsin peptides not accounted for. FAIL.", file = logfile)
if (len(trypn_uniq_peptides) - len(trypn_KR)/2) == trypn_only + both: print("All TrypN peptides accounted for. PASS.", file = logfile)
else: print("Some TrypN peptides not accounted for. FAIL.", file = logfile)
if other != 0: print("Some peptides were not properly filtered. Please check the code and input.", file = logfile)

logfile.close()