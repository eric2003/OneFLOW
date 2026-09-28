from fractions import Fraction

def intregral_xi(x,j):
    return (x**(j+1))/(j+1)
    
def calxi(j):
    half = Fraction(1,2)
    xia = j - half
    xib = j + half
    return xia, xib

jst = -2
jed = 2

for j in range(jst, jed+1):
    xia, xib = calxi(j)
    print(f'j={j:>2} xia,b=[{str(xia):>4},{str(xib):>4}]')
    ist = 0
    ied = 4
    for i in range(ist, ied+1):
        v1 = intregral_xi(xia,i)
        v2 = intregral_xi(xib,i)
        print(f'    v1={v1} v2={v2} diff ={v2-v1}')

    
