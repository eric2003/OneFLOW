import numpy as np
import csv
import sys
import os

def get_sorted_indices(arr1, arr2):
    # 创建一个字典来存储原始数组的索引  
    index_dict = {value: index for index, value in enumerate(arr1)}  
    # 使用字典来构建 index_map
    index_map = [index_dict[item] for item in arr2]  
    return index_map

nvar = len(sys.argv)
print('nvar=',nvar)
print('sys.argv=',sys.argv)

nt = 2000
if nvar >= 2:
    mt = sys.argv[1]
    nt = int(mt)
    print('nt=',nt)
    
filename_src = 'field_final'+str(nt)+'.csv'
filename_tgt = 'field_final_x_tgt.csv'
print("filename_src=",filename_src)
print("filename_tgt=",filename_tgt)
    
with open(filename_src, newline='') as csvfile:
    readCSV = csv.reader(csvfile, delimiter= ' ')
    icount = 0
    for row in readCSV:
        icount += 1

ni = icount
print("ni=",ni)

nm = 3
q = np.zeros( (ni, nm ) )
x = np.zeros( (ni) )
xt = np.zeros( (ni) )

with open(filename_src, newline='') as csvfile:
    readCSV = csv.reader(csvfile, delimiter= ' ')
    i = 0
    for row in readCSV:
        x[i]    = float(row[0])
        q[i][0] = float(row[1])
        q[i][1] = float(row[2])
        q[i][2] = float(row[3])
        i += 1
        
with open(filename_tgt, newline='') as csvfile:
    readCSV = csv.reader(csvfile, delimiter= ' ')
    i = 0
    for row in readCSV:
        xt[i]    = float(row[0])
        i += 1
        
#exit()        
        
#sort
sorted_indices = get_sorted_indices(x, xt)
qt=q[sorted_indices,:]

data = []

for i in range(ni):
    ll = []
    ll.append("{:.25f}".format(xt[i]))
    ll.append("{:.25f}".format(qt[i,0]))
    ll.append("{:.25f}".format(qt[i,1]))
    ll.append("{:.25f}".format(qt[i,2]))
    data.append(ll)
    
basename = os.path.basename(filename_src)
file_name, _ = os.path.splitext(basename)

outfilename = file_name + '.tgt'

print("outfilename=",outfilename)

with open(outfilename, 'w', newline='', encoding='utf-8') as csvfile:  
    writer = csv.writer(csvfile, delimiter=' ')  
    writer.writerows(data)  

