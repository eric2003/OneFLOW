import numpy as np
import csv
import sys
import os

nvar = len(sys.argv)
print('nvar=',nvar)
print('sys.argv=',sys.argv)

nt = 2000
if nvar >= 2:
    mt = sys.argv[1]
    nt = int(mt)
    print('nt=',nt)

with open('field_final0.csv', newline='') as csvfile:
    readCSV = csv.reader(csvfile, delimiter= ' ')
    icount = 0
    for row in readCSV:
        icount += 1

ni = icount
print("ni=",ni)

nm = 3
q = np.zeros( (ni, nm ) )
x = np.zeros( (ni) )

filename = 'field_final'+str(nt)+'.csv'

with open(filename, newline='') as csvfile:
    readCSV = csv.reader(csvfile, delimiter= ' ')
    i = 0
    for row in readCSV:
        x[i]    = float(row[0])
        q[i][0] = float(row[1])
        q[i][1] = float(row[2])
        q[i][2] = float(row[3])
        i += 1

#sort
sorted_indices = np.argsort(x)
xs=x[sorted_indices]
qs=q[sorted_indices,:]

data = []

for i in range(ni):
    ll = []
    ll.append("{:.25f}".format(xs[i]))
    ll.append("{:.25f}".format(qs[i,0]))
    ll.append("{:.25f}".format(qs[i,1]))
    ll.append("{:.25f}".format(qs[i,2]))
    data.append(ll)
    
basename = os.path.basename(filename)
file_name, _ = os.path.splitext(basename)

outfilename = file_name + '.bak'

print("outfilename=",outfilename)

with open(outfilename, 'w', newline='', encoding='utf-8') as csvfile:  
    writer = csv.writer(csvfile, delimiter=' ')  
    writer.writerows(data)  

