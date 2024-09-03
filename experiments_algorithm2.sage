import numpy as np
import random
from sage.stats.distributions.discrete_gaussian_integer import DiscreteGaussianDistributionIntegerSampler
import time


####### basic subroutine algorithms 

## These algorithms samples vector following uniform/Discrete Gaussian distributions
def rand_vector(n,q): return vector([randint(-q,q) for i in range(n)])
def noise_vector(n,D): return vector([D() for i in range(n)])


## This algorithm outputs a coefficientwise modulo reduced ring element.
def mod(a,q):
   temp = a%q
   if temp > (q-1)/2:
      temp -= q
   return temp


def ringmod(s,q,n,zeta):
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
def szf(f,n,zeta):
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





######## Default Parameter Setup for BGV
sigma = 3.2    # standard deviation of the noise distribution 
t = 2^20       # message modulus



def Reaction_attack(n,logq):
    #n  :  Secret dimension, it should be a power of two
    #logq : bit of underlying modulus
    if is_power_of_two(n) == false:
        print('You should take an input n as a power of two!!')
        return 
    ## underlying setup
    q = next_prime(2^logq)
    p = smallest_prime(n)
    R = CyclotomicField(2*n)
    zeta = R.gen()
    D = DiscreteGaussianDistributionIntegerSampler(sigma=sigma)
    

    ####### Construct a set of ghat used in the Algorithm 2.
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
        ghat_list += [ringmod(F1/factor_list[i],p,n,zeta)]

    ###### Algorithm 2 test
    Iteration_number = 100
    Success_number = 0
    for tes in range(Iteration_number):
        print('tes=',tes)
        ###### Sample a secret vector
        coef_s = noise_vector(n,D)
        s = R(list(coef_s))
        sp = Rp(list(s))
        sol = []

        for i in range(n):
            sol += [ZZ(sp%factor[i][0])]
        ######## Algorithm 2 main
        sol_list = []
        for i in range(n):
            ####### Sample a ciphertext         
            e = R(list(noise_vector(n,D)))
            a = R(list(rand_vector(n,(q-1)/2)))
            b = ringmod(a*s +t*e, q, n,zeta )
            u =0 
            while true == 1:
                if szf(ringmod(  b+ t*ringround (ghat_list[i]*u*q/p ,zeta, n)  - (a+t*ringround(ghat_list[i]*q/p , zeta, n))*s ,q,n,zeta),n,zeta) < q/4 :
                    sol_list += [u]
                    break
                u += 1


        if sol_list == sol:     ## If the secret vector is successfully recovered, it gives 1. 
            Success_number += 1
    return Success_number/ Iteration_number









