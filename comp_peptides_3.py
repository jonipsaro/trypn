#Jonathan Ipsaro
#Last reviewed: January 11, 2020
#
#This script takes csv files from MASCOT and compares three replicate sample preparations.
#
#This script will analyze PEPTIDES and counts peptides as unique if they have different modifications
#
#CSV files from MASCOT should have the peptide sequence in the 24th position (python index 23) and 
#modifications in the 26th position (python index 25).
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

#Specify the three files to be analyzed
fin1 = open("CSV/TrypN_replicate1.csv", "r")
fin2 = open("CSV/TrypN_replicate2.csv", "r")
fin3 = open("CSV/TrypN_replicate3.csv", "r")

peptides1 = set()
peptides2 = set()
peptides3 = set()

#Build a set of all peptides in the first input file
#Each element in the set is a tuple of two strings (sequence, modifications)
build = False
for line in fin1:
	if line.find("prot_hit") == 0: 
		build = True
		line = fin1.readline()

	if build == True:
		line = CleanLine(line)
		seq = line.split(",")[23]
		mod = line.split(",")[25].strip('\"')
		if mod == "": mod = "None"

		peptide = (seq, mod)
		peptides1.add(peptide)

fin1.close()


#Build a set of all peptides in the second input file
build = False
for line in fin2:
	if line.find("prot_hit") == 0: 
		build = True
		line = fin2.readline()

	if build == True:
		line = CleanLine(line)
		seq = line.split(",")[23]
		mod = line.split(",")[25].strip('\"')
		if mod == "": mod = "None"

		peptide = (seq, mod)
		peptides2.add(peptide)
		
fin2.close()


#Build a set of all peptides in the third input file
build = False
for line in fin3:
	if line.find("prot_hit") == 0: 
		build = True
		line = fin3.readline()

	if build == True:
		line = CleanLine(line)
		seq = line.split(",")[23]
		mod = line.split(",")[25].strip('\"')
		if mod == "": mod = "None"

		peptide = (seq, mod)
		peptides3.add(peptide)
		
fin3.close()


#Calculate overlaps or Venn Diagram
intersect12 = peptides1.intersection(peptides2)
intersect13 = peptides1.intersection(peptides3)
intersect23 = peptides2.intersection(peptides3)
intersect123 = intersect12.intersection(peptides3)

unique1 = peptides1.difference(peptides2.union(peptides3))
unique2 = peptides2.difference(peptides1.union(peptides3))
unique3 = peptides3.difference(peptides1.union(peptides2))


#Output
logfile = open("peptide_comparison.log", "w+")
print("TOTALS:", file = logfile)
print("Number of peptides in replicate 1 (total):", len(peptides1), file = logfile)
print("Number of peptides in replicate 2 (total):", len(peptides2), file = logfile)
print("Number of peptides in replicate 3 (total):", len(peptides3), file = logfile)

print("\nPRESENT IN TWO REPLICATES:", file = logfile)
print("Number of peptides shared in replicates 1 and 2:", len(intersect12), file = logfile)
print("Number of peptides shared in replicates 1 and 3:", len(intersect13), file = logfile)
print("Number of peptides shared in replicates 2 and 3:", len(intersect23), file = logfile)

print("\nPRESENT IN THREE REPLICATES:", file = logfile)
print("Number of peptides shared in replicates 1, 2, and 3:", len(intersect123), file = logfile)

print("\nPRESENT ONLY IN ONE REPLICATE:", file = logfile)
print("Number of peptides unique to replicate 1:", len(unique1), file = logfile)
print("Number of peptides unique to replicate 2:", len(unique2), file = logfile)
print("Number of peptides unique to replicate 3:", len(unique3), file = logfile)

print("\nPRESENT ONLY IN TWO REPLICATES:", file = logfile)
print("Number of peptides in ONLY replicates 1 and 2:", len(intersect12)-len(intersect123), file = logfile)
print("Number of peptides in ONLY replicates 1 and 3:", len(intersect13)-len(intersect123), file = logfile)
print("Number of peptides in ONLY replicates 2 and 3:", len(intersect23)-len(intersect123), file = logfile)

#Output for Venn Diagram Plot
print("A =", len(unique1), ", B = ", len(unique2), ", C = ", len(unique3), ", \"A&B\" =", len(intersect12)-len(intersect123), ", \"A&C\" =", len(intersect13)-len(intersect123), ", \"B&C\" =", len(intersect23)-len(intersect123), ", \"A&B&C\" =", len(intersect123))