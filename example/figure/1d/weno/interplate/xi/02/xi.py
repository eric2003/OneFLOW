from fractions import Fraction

def calxi(x,j):
    return (x**(j+1))/(j+1)

jst = 0
jed = 4
x = Fraction(1,2)
for j in range(0, jed+1):
    v = calxi(x,j)
    print(f'Intergral[({x})^{j}]={v}')
    
