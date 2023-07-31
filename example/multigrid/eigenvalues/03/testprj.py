import numpy as np
import matplotlib.pyplot as plt
    
n = 12
h = 1.0/n
dx = h

print("dx=",dx)

x = np.zeros( 2*n+1 )
ww=np.array([1.0/3.0,0.5,2.0/3.0,1.0])
mylabel_list=np.array(["1/3","1/2","2/3","1"])
mycolor_list=np.array(['black','blue','red','green'])
mymarker_list=np.array(['o','s','v','^'])
num = len(ww)
print("num=",num)
eigen = np.zeros( (num,2*n+1) )

for i in range(0, 2*n+1):
    x[i] = i * dx

for j in range(0, num):
    w=ww[j]
    for k in range(0, 2*n+1):
      value = np.sin( k * np.pi / (2*n) )
      v2 = value * value
      eigen[j][k] = 1.0 - 2 * w * v2
    
id = np.zeros(2*n+1)
for i in range(0, 2*n+1):
   id[i] = i
   
fig = plt.figure("OneFLOW-CFD Solver Eigenvalues of the iteration matrix", figsize=(6, 4), dpi=100)  
ax = fig.add_subplot()

#ax.spines['left'].set_position('center')
ax.spines['bottom'].set_position('center')
ax.spines['right'].set_color('none')
ax.spines['top'].set_color('none')
ax.set_xlabel("k")
ax.set_ylabel(r'$\lambda_{k}(P_{\omega})$',rotation='horizontal')
ax.yaxis.set_label_coords(0.0, 1.0)
ax.xaxis.set_label_coords(1.01, 0.51)

for j in range(0, num):
    ym = eigen[j][-1]
    print("ym=",ym)
    mylabel=r'$\omega$ = '+mylabel_list[j]
    ax.plot(id,eigen[j],color=mycolor_list[j],linewidth=1,marker=mymarker_list[j],markerfacecolor='none',label=mylabel)
    x = id[-1]
    y = eigen[j][-1]
    ax.annotate(mylabel, xy=(0,0), xytext=(x-1,y+0.1), textcoords='data' )  
    #ax.legend()

fig.tight_layout()
plt.show()
