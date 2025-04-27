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




def ringmod(s,q,n,Power_zeta):
   temp = 0
   for i in range(n):
      temp += Power_zeta[i] * mod(s[i],q)
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
def ringround(b,Power_zeta,n):
   temp = 0
   for i in range(ZZ(n)):
      temp += round(b[i])*Power_zeta[i]
   return temp







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
    Power_zeta = []
    for i in range(n):
        Power_zeta += [zeta^i]
    D = DiscreteGaussianDistributionIntegerSampler(sigma=sigma)
    threshold = q/5

    ####### Construct a set of ghat used in the Algorithm 2.

    
    Rp = PolynomialRing(GF(p),'x')
    x = Rp.gen()
    f = x^n +1
    factor = f.factor()
    factor_list = []
    for i in range(n):
        factor_list += [factor[i][0]]
    F1 = prod(factor_list)
    ghat_list = []
    ghat_list2 = []
    for i in range(n):
        ghat_list += [ringmod(R(Rp(F1/factor_list[i])),p,n,Power_zeta)]
        ghat_list2 += [ringtovector(t*ghat_list[i]*q/p,n,zeta)]


    ###### Algorithm 2 test
    Iteration_number = 10
    Success_number = 0
    for tes in range(Iteration_number):
        
        ###### Sample a secret vector
        coef_s = noise_vector(n,D)
        s = R(list(coef_s))
        sp = Rp(list(s))
        sol = []

        for i in range(n):
            sol += [ZZ(sp%factor[i][0])]
        ######## Algorithm 2 main
        sol_list = []

        list_m = []
        list_m = [round(random.uniform(-1, 1) * t) for i in range(n)]
        m = R(list_m)

        for i in range(n):
            
            ####### Sample a ciphertext
            liste = list(noise_vector(n,D))          
            e = R(liste)
            a = R(list(rand_vector(n,(q-1)/2)))
            b = ringmod(a*s +e + m, q, n,Power_zeta)
            u = 0 
            
            #### We observe that b_u - a_u*s = b- a*s + ghat_list2[i]*sol[i] - ghat_list2[i] * u. We thus compute b- a*s + ghat_list2[i]*sol[i] as common  
            common_vec = vector(liste) - ghat_list2[i]*sol[i]
            
            temp_vec = ghat_list2[i]
            temp_vec2 = common_vec- temp_vec
            while u < p:
                temp_vec2 += temp_vec
                # temp_vec2 = ringmod(temp_vec2, q, n, Power_zeta)
                # print(szf(temp_vec2, n, zeta) < q/5)

                for j, entry in enumerate(temp_vec2):
                    if mod(round(entry),q) > threshold:
                        u += 1
                        break
                else:
                    sol_list += [u]
                    break
#        print(sol)
#        print(sol_list)
        if sol_list == sol:     ## If the secret vector is successfully recovered, it gives 1. 
            Success_number += 1
        print(f"Running: {tes:.1f} times and {Success_number:.1f} times success")
    return Success_number/ Iteration_number









sigma = 3.2    # standard deviation of the noise distribution 
t = 2^30    # message modulus
logq = 35
n = 2^6
start_time = time.time()
Reaction_attack(n,logq)
end_time = time.time()
print(f"Execution time: {end_time - start_time:.6f} seconds")



# ######## Default Parameter Setup for CKKS
# sigma = 3.2    # standard deviation of the noise distribution 
# t = 2^30       # Delta
# n = 2^13
# logq = 35
# start_time = time.time()
# Reaction_attack(n,logq)
# end_time = time.time()
# print(f"Execution time: {end_time - start_time:.6f} seconds")
