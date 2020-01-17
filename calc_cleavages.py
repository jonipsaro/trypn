#Jonathan Ipsaro
#Last reviewed: January 11, 2020
#
#This script takes csv files from MASCOT and determines the frequency of each N-terminal amino acid
#in each unique peptide.
#Unique peptides are defined as those that have unique sequences and the same modifications.
#
#N-terminal cleaveage frequencies are determined by using both the N-terminal residue of the
#mapped peptide and the residue immediately following the peptide.
#
#This script will work on several files specified in the filenames list in the specified directory.
#
#CSV files from MASCOT should have the peptide sequence in the 24th position (python index 23); 
#the previous amino acid should be in the 23rd position; the subsequent amino acid in the 25th position;
#modifications should be in the 26th position.
#
#A full output file ending with the name _cleavages.log will have the percent cleavages for each
#amino acid in each file as well as summary statistics (total unique peptides, total cleavages, 
#the number of internal peptides [these contribute 2 sites], or at the N- or C-terminus of the protein
#[these contribute only 1 cleavage each].
#
#Another output file will be generated that only displays the K and R specificity for easy graphing
#in other software.
#
#Finally, an additional ouput file ending with the name _table.log is a tabulated format that can be read
#into R for logo drawing.
#
#This script uses Python3 and depends on pandas.

import pandas as pd

#Function to clean up MASCOT CSV format
def CleanLine(line):
	quote = False

	for i in range(len(line)):
		if quote == False and line[i] == "\"": quote = True
		elif quote == True and line[i] == "\"": quote = False

		if ((quote == True) and (line[i] == ",")): line = line[:i]+":"+line[i+1:]

	return line

#Set file names and directory
filenames = ["filename1.csv", "filename2.csv"]		#Specify file names as strings in a list
directory = "directory/"									#Specify the path to the directory where the files are
													#be sure to include the trailing "/" in the path
prefix = "prefix"									#Specify a prefix to be used for all output files

#Set up summary output files
spec_file = open(prefix+"_summary.log", "w+")
print("Filename.csv\tK\tR\tOther", file=spec_file)

spec_uniq_file = open(prefix+"_summary2.log", "w+")
print("Filename.csv\tKR_percent\tunique_peptides", file=spec_uniq_file)

#Declare variables and set up output data frame
logfile = open(prefix+"_cleavages.log", "w+")
aa = ["A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y", "-"]
filenames_short = []
for each in filenames: filenames_short.append(each[0:each.find(".")])
if "-" in aa: aa_only = aa.remove("-")
else: aa_only = aa

summary_df = pd.DataFrame(index=aa_only, columns=filenames_short)

#Begin processing, send output to the logfile
print("Processing", len(filenames), "files...")
print("Filename.csv", file=logfile)
print("P1'_residue\tPercent_cleavages", file=logfile)
count = 0
#Go through each file and determine all cleaveage events
for file in filenames:
	fin = open(directory+file, "r")			#Open the file
	filename_short = file[0:file.find(".")]	#Make the short file name for the dataframe output

	all_peptides = set()	#Refresh the variables
	internal = set()
	nterminal = set()
	cterminal = set()
	bothtermini = set()

	aa_internal = {}		#Set up a dictionary for each type of peptide
	aa_nterminal = {}		#Each dictionary will have the format {amino_acid: [peptide1, peptide2, ...]}
	aa_cterminal = {}
	aa_bothtermini = {}
	quick_test = {}

	for each in aa:			#Populate the dictionary keys with the amino acids
		aa_internal[each] = []
		aa_nterminal[each] = []
		aa_cterminal[each] = []
		aa_bothtermini[each] = []

	build = False						#Iterate through the file until the "prot_hit" line is found
	for line in fin:					#Lines preceding this are in the header
		if line.find("prot_hit") == 0: 
			build = True
			line = fin.readline()

		if build == True:				#Then start reading in lines for the analysis
			line = CleanLine(line)

			prev_aa = line.split(",")[22]			#Extract the relevant parts of each line
			seq = line.split(",")[23]
			next_aa = line.split(",")[24]
			mod = line.split(",")[25].strip('\"')

			if mod == "": mod = "None"

			peptide = (seq, mod)		#Define each unique peptide

 
			#Determine if this is a new, unique internal peptide
			#If so, append that peptide to the appropriate {amino_acid:[peptides]}
			if (peptide not in all_peptides) and (prev_aa != "-") and (next_aa != "-"):
				internal.add(peptide)
				aa_internal[peptide[0][0]].append(peptide)
				aa_internal[next_aa].append(peptide)

			#Or if it is an N-terminal peptide
			elif (peptide not in all_peptides) and (prev_aa == "-") and (next_aa != "-"):
				nterminal.add(peptide)
				aa_nterminal[next_aa].append(peptide)

			#Or if it is a C-terminal peptide
			elif (peptide not in all_peptides) and (prev_aa != "-") and (next_aa == "-"):
				cterminal.add(peptide)
				aa_cterminal[peptide[0][0]].append(peptide)

			#Or if it is both N- and C-terminal
			elif (peptide not in all_peptides) and (prev_aa == "-") and (next_aa == "-"):
				bothtermini.add(peptide)

			all_peptides.add(peptide)

	#Determine the total number of cleavages:
	#Internal peptides provide 2, N- and C-terminal peptides each provide 1
	total_cleavages = len(internal)*2 + len(nterminal) + len(cterminal)

	print("==============================================", file=logfile)
	print(file, file=logfile)

	#Calculate the fraction of cleavages for each amino acid
	K_percent = 0
	R_percent = 0
	Other_percent = 0

	for each in aa:
		if each != "-":
			aa_total = len(aa_internal[each]) + len(aa_nterminal[each]) + len(aa_cterminal[each])
			aa_percent = (aa_total/total_cleavages)*100
			summary_df.at[each, filename_short] = aa_percent
			print(each, aa_percent, sep="\t", file=logfile)

			if each == "K": K_percent = aa_percent
			elif each == "R": R_percent = aa_percent
			else: Other_percent += aa_percent

	#Output to the main log file
	print(file, "TOTAL UNIQUE PEPTIDES:", len(all_peptides), sep="\t", file=logfile)
	print(file, "TOTAL CLEAVAGES:", total_cleavages, sep="\t", file=logfile)
	print(file, "N-TERMINAL PEPTIDES:", len(nterminal), sep="\t", file=logfile)
	print(file, "INTERNAL PEPTIDES:", len(internal), sep="\t", file=logfile)
	print(file, "C-TERMINAL PEPTIDES:", len(cterminal), sep="\t", file=logfile)
	print("==============================================", file=logfile)
	fin.close()

	#Output to the KR specificity file
	print(file, K_percent, R_percent, Other_percent, sep="\t", file=spec_file)

	#Output to the specificity_unique file
	print(file, K_percent+R_percent, len(all_peptides), sep="\t", file=spec_uniq_file)

	print(file, "... DONE")

logfile.close()
spec_file.close()
spec_uniq_file.close()


#Output the data frame for logo plotting in R
summary_file = open(prefix+"_table.log", "w+")
summary_file.write(summary_df.to_csv(sep="\t"))
summary_file.close()