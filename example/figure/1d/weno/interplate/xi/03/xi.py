from fractions import Fraction

def integral_xi(x, j):
    return (x ** (j + 1)) / (j + 1)

print(f"{'j':>3}  {'interval':>14}  " + "  ".join(f"i={i:^8}" for i in range(5)))
print("-" * 72)

vv = [-1,0,1]
for j in vv:
    xia = Fraction(j) - Fraction(1,2)
    xib = Fraction(j) + Fraction(1,2)
    print(f"{j:>3}  [{xia:>4},{xib:>4}]  ", end="")
    for i in range(5):
        val = integral_xi(xib, i) - integral_xi(xia, i)
        s = "0" if val == 0 else str(val)
        print(f"{s:^11}", end="")
    print()