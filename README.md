#PHASESPERM
A program to phase sperm. 

##Before you start
sudo apt-get install python-scipy

Data should be formated :

A simple test script generates data :

Pipeline:
 * simulatesperm.py
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
