#PHASESPERM
A program to phase sperm. 

##Before you start
This python script requires scipy and numpy. If you do not have these installed you can type:

	sudo apt-get install python-scipy
	sudo apt-get install python-numpy

phasesperm.py expects data to be in a tab delimited file with 3 columns containing the scaffold, position and consensus information, and then one column for each sperm sequenced.
The symbol '-' is used for missing data. For example, if 4 sperm were sequenced, then a few lines of this file might look like:

	scaffold\_49     1043    T       -       -       A	T
	scaffold\_49     2240    T       T       -       T	A
	scaffold\_49     2363    G       G       -       -	C

If you want to generate test data in this format, type:

	simulatesperm.py

Pipeline:
 * trimsnp.py test.txt 10 > test-min10.txt
 * phasesperm.py test-min10.txt
 * cd MST/
 * make
 * ./MST.exe ../MST-data.txt ../MST-out.txt
 * cd ../
 * addmarkers.py MST-out.txt > MST-out-wmakers.txt
 * trimtopairs.py outR.csv > pairs.txt
 * mltrack.py pairs.txt
 * Rscipt graphoutR.R
