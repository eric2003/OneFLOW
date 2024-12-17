def thomas_algorithm(a, b, c, d):
    n = len(d)
    c_prime = [0] * n
    d_prime = [0] * n
    x = [0] * n
    c_prime[0] = c[0] / b[0]
    d_prime[0] = d[0] / b[0]
    
    for i in range(1, n):
        coef = 1.0 / ( b[i] - a[i] * c_prime[i-1] )
        c_prime[i] = c[i] * coef
        d_prime[i] = ( d[i] - a[i] * d_prime[i-1] ) * coef

    x[n-1] = d_prime[n-1]

    for i in range(n-2, -1, -1):
        x[i] = d_prime[i] - c_prime[i] * x[i+1]

    return x

a = [0, -1, -1, -1, -1]
b = [2, 2, 2, 2, 2]
c = [-1, -1, -1, -1, 0]
d = [1, 1, 1, 1, 1]

x = thomas_algorithm(a, b, c, d)
print(x)
