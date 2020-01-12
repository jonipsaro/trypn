#Jonathan Ipsaro
#Last reviewed: January 11, 2019
#
#This script takes csv files from MASCOT and compares three replicate sample preparations.
#
#This script will analyze PROTEINS based on accession number.
#
#CSV files from MASCOT should have the peptide sequence in the 4th position (python index 3).
#
#This script uses Python3.

#Function to clean up MASCOT CSV format
def CleanLine(line):
	quote = False

	for i in range(len(line)):
		if quote == False and line[i] == "\"": quote = True
		elif quote == True and line[i] == "\"": quote = False

		if ((quote == True) and (line[i] == ",")): line = line[:i]+":"+line[i+1:]

	return line


fin1 = open("CSV/TrypN_proteins1.csv", "r")
fin2 = open("CSV/TrypN_proteins2.csv", "r")
fin3 = open("CSV/TrypN_proteins3.csv", "r")

proteins1 = set()
proteins2 = set()
proteins3 = set()

#Build a set of all proteins in the first input file
#Each element in the set is a tuple of two strings (sequence, modifications)
build = False
for line in fin1:
	if line.find("Family") == 0: 
		build = True
		line = fin1.readline()

	if build == True:
		line = CleanLine(line)
		protein = line.split(",")[3]
		proteins1.add(protein)

fin1.close()


#Build a set of all proteins in the second input file
build = False
for line in fin2:
	if line.find("Family") == 0: 
		build = True
		line = fin2.readline()

	if build == True:
		line = CleanLine(line)
		protein = line.split(",")[3]
		proteins2.add(protein)
		
fin2.close()


#Build a set of all proteins in the third input file
build = False
for line in fin3:
	if line.find("Family") == 0: 
		build = True
		line = fin3.readline()

	if build == True:
		line = CleanLine(line)
		protein = line.split(",")[3]
		proteins3.add(protein)
		
fin3.close()


#Calculate overlaps or Venn Diagram
intersect12 = proteins1.intersection(proteins2)
intersect13 = proteins1.intersection(proteins3)
intersect23 = proteins2.intersection(proteins3)
intersect123 = intersect12.intersection(proteins3)

unique1 = proteins1.difference(proteins2.union(proteins3))
unique2 = proteins2.difference(proteins1.union(proteins3))
unique3 = proteins3.difference(proteins1.union(proteins2))


#Output
logfile = open("protein_comparison.log", "w+")

print("TOTALS:", file = logfile)
print("Number of proteins in replicate 1 (total):", len(proteins1), file = logfile)
print("Number of proteins in replicate 2 (total):", len(proteins2), file = logfile)
print("Number of proteins in replicate 3 (total):", len(proteins3), file = logfile)

print("\nPRESENT IN TWO REPLICATES:", file = logfile)
print("Number of proteins shared in replicates 1 and 2:", len(intersect12), file = logfile)
print("Number of proteins shared in replicates 1 and 3:", len(intersect13), file = logfile)
print("Number of proteins shared in replicates 2 and 3:", len(intersect23), file = logfile)

print("\nPRESENT IN THREE REPLICATES:", file = logfile)
print("Number of proteins shared in replicates 1, 2, and 3:", len(intersect123), file = logfile)

print("\nPRESENT ONLY IN ONE REPLICATE:", file = logfile)
print("Number of proteins unique to replicate 1:", len(unique1), file = logfile)
print("Number of proteins unique to replicate 2:", len(unique2), file = logfile)
print("Number of proteins unique to replicate 3:", len(unique3), file = logfile)

print("\nPRESENT ONLY IN TWO REPLICATES:", file = logfile)
print("Number of proteins in ONLY replicates 1 and 2:", len(intersect12)-len(intersect123), file = logfile)
print("Number of proteins in ONLY replicates 1 and 3:", len(intersect13)-len(intersect123), file = logfile)
print("Number of proteins in ONLY replicates 2 and 3:", len(intersect23)-len(intersect123), file = logfile)

print("Vector for R script:")
print("A =", len(unique1), ", B = ", len(unique2), ", C = ", len(unique3), ", \"A&B\" =", len(intersect12)-len(intersect123), ", \"A&C\" =", len(intersect13)-len(intersect123), ", \"B&C\" =", len(intersect23)-len(intersect123), ", \"A&B&C\" =", len(intersect123))