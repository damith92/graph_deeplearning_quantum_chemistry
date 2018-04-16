import csv
 
tot_smiles = []
tot_solubility = []

rownum = 0

rfile = open('delaney.csv', newline='')
csv_reader = csv.reader(rfile)
 
for row in csv_reader:
    if rownum != 0 :
        tot_smiles.append(row[9])
        tot_solubility.append(float(row[8]))
        
    rownum += 1
    
rfile.close()
