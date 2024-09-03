import numpy as np
import random
from sage.stats.distributions.discrete_gaussian_integer import DiscreteGaussianDistributionIntegerSampler
import time


####### basic subroutine algorithms 

## These algorithms samples vector following uniform/Discrete Gaussian distributions
def rand_vector(n,q): return vector([randint(-q,q) for i in range(n)])
def noise_vector(n): return vector([D() for i in range(n)])


## This algorithm outputs a coefficientwise modulo reduced ring element.
def mod(a,q):
   temp = a%q
   if temp > (q-1)/2:
      temp -= q
   return temp


def ringmod(s,q):
   temp = 0
   for i in range(n):
      temp += zeta^i * mod(s[i],q)
   return temp



## This algorithm finds the smallest prime p such that p = 1 mod 2n.
def smallest_prime(n):
    temp = 2*n 
    while (temp +1).is_prime() !=1:
        temp += 2*n
    return temp+1


## These algorithms are used for computing the size of a ring element.
def szf(f):
   return RR(ringtovector(f,n,zeta).norm(infinity))

def ringtovector(f,n,zeta):
   tempvector = vector(QQ,ZZ(n))
   for i in range(n):
      tempvector[i] = f[i]
   return tempvector


## This algorithm outputs a coefficientwise rounded ring element.
def ringround(b,zeta,n):
   temp = 0
   for i in range(ZZ(n)):
      temp += (b[i].round())*zeta^i
   return temp





######## Parameter Setup for BGV
n = 2^5     # Secret dimension
logq = 100   # bit of underlying modulus
q = next_prime(2^logq)
sigma = 3.2    # standard deviation of the noise distribution 
t = 2^20
p = smallest_prime(n)



####### Sample a ciphertext 
R = CyclotomicField(2*n)
zeta = R.gen()

D = DiscreteGaussianDistributionIntegerSampler(sigma=sigma)
coef_s = noise_vector(n)
s = R(list(coef_s))
e = R(list(noise_vector(n)))
a = R(list(rand_vector(n,(q-1)/2)))
b = ringmod(a*s +t*e, q )


####### Construct a set of ghat used in the Algorithm 2 and 3. 
Rp = PolynomialRing(GF(p),'x')
x = Rp.gen()
f = x^n +1
factor = f.factor()
factor_list = []
for i in range(n):
    factor_list += [R(factor[i][0])]
F1 = prod(factor_list)
ghat_list = []
for i in range(n):
    ghat_list += [ringmod(F1/factor_list[i],p)]





###### construct a solution for easy comparing
sp = Rp(list(s))
sol = []

for i in range(n):
    sol += [ZZ(sp%factor[i][0])]


######## Algorithm 2

sol_list = []


for i in range(n):
    u =0 
    while true == 1:
        if szf(ringmod(  b+ t*ringround (ghat_list[i]*u*q/p ,zeta, n)  - (a+t*ringround(ghat_list[i]*q/p , zeta, n))*s ,q)) < q/4 :
            sol_list += [u]
            break
        u += 1


if sol_list == sol:
    print('Secret key is correctly recovered')
