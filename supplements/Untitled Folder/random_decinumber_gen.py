import random

outFile = open("phenotype_values.txt","w")
for i in range(50):
	outFile.write(str(random.randint(1,110)/100))
	outFile.write("\n")

outFile.close();
