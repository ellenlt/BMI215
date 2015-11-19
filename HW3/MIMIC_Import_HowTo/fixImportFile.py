import sys
import string

def fixImportFile(infileName,outfileName):
	infile = open(infileName,'rb')
	outfile = open(outfileName,'wb')
	lines = infile.readlines()
	newLines = []
        PathToPatientFiles = '~/' # REPLACE THIS WITH YOUR OWN PATH TO THE PATIENT FILES FOLDER
	for line in lines:
		if line.find("/Users/estrandb") >= 0:
			newLine = line.split("PATH/")
			newLine = string.join([newLine[0], PathToPatientFiles, newLine[1]],sep="")
			newLines.append(newLine)
		else:
			newLines.append(line)
	infile.close()
	for line in newLines:
		outfile.write(line)
	outfile.close()

if len(sys.argv) >=2:
	fixImportFile(sys.argv[1],sys.argv[2])

