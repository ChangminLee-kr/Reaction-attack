# Reaction Attack Simulation

This repository contains a simple SageMath script (`experiments_algorithm2.sage`) that simulates the success probability of a reaction attack. The script is designed to be loaded into SageMath, and it provides a function `Reaction_attack(n, logq)` to compute the success probability of an attack based on given parameters.

## Prerequisites

- [SageMath](https://www.sagemath.org/) installed on your system.

## Getting Started

To use the script, you first need to load it into SageMath. Follow these steps:

1. **Open SageMath**: Start your SageMath environment.

2. **Load the script**: Use the following command to load the `experiments_algorithm2.sage` script:
   load('experiments_algorithm2.sage')

3. **Run the Reaction Attack Simulation**: After loading the script, you can call the Reaction_attack function to compute the success probability. Note that the parameter n should be a power of two.
   
    Reaction_attack(n, logq)

n: This parameter should be a power of two (e.g., 2, 4, 8, 16, etc.).
logq: This parameter represents the logarithm of the value q used in the simulation.




# Example Usage

Load the script
load('experiments_algorithm2.sage')

Run the function with n = 32 (which is a power of two) and logq = 100
success_probability = Reaction_attack(8, 5)

Output the result
print(success_probability)
