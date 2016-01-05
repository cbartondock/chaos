import fractions
def simplest_rationals_between(a, b,m):
    pairs = list(set([(p,q) for q in range(1,m) for p in  range(1,m)]))
    pairs =[e for e in pairs if float(e[0])/float(e[1]) >= a and float(e[0])/  float(e[1]) <=b and fractions.gcd(e[0],e[1])==1]
    pairs.sort(key = lambda e: e[0]+e[1])
    return pairs

