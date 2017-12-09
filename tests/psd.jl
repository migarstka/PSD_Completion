# N = 10;
# F =svdfact(rand(N,N))
# C=F[:U]*diagm(rand(N))*F[:U]'

N = 8
A = rand(-5:5,N,N)
A = 0.5 *(A+A')
F =eigfact(A)
Λ = diagm(F[:values])
Q = F[:vectors]
A = Q*max(Λ,0)*Q'

