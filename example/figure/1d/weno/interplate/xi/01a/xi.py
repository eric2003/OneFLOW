from fractions import Fraction

def calxi(j):
    half = Fraction(1,2)
    xia = j - half
    xib = j + half
    return xia, xib

jst = -2
jed = 2

for j in range(jst, jed+1):
    xia, xib = calxi(j)
    #print(f'j={j:>2} xia,b=[{xia},{xib}]')
    print(f'j={j:>2} xia,b=[{str(xia):>4},{str(xib):>4}]')

