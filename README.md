# Optimization-Interior-point-method

Optimization Algorithms for Linear Programming
This repository contains three algorithms for solving linear programming problems: Central Path Fixed, Central Path Adaptive, and Mehrotra Predictor. These algorithms are implemented in MATLAB and can be used to find optimal solutions to linear programming problems.

# Algorithms Overview

# Central Path Fixed

The Central Path Fixed algorithm aims to find the optimal solution to a linear programming problem by iteratively approaching the central path. The central path is the set of feasible solutions that minimize the objective function while satisfying the constraints. The algorithm follows these steps:


Initialization: Initialize the decision variables, slack variables, and Lagrange multipliers.

Iterative Process: Iterate until convergence or a maximum number of iterations is reached.

a. Compute the residuals (rc, rb, rxs) and the affine search direction (deltax_affine, deltas_affine).

b. Compute the step lengths (alpha_prime and alpha_dual) using the affine search direction.

c. Perform the corrector step to update the decision variables.

Termination: Stop the iterations when the termination criterion is met (e.g., tolerance on the central path measure).

#Central Path Adaptive
The Central Path Adaptive algorithm is an improvement over the Central Path Fixed algorithm. It dynamically adjusts the step lengths (alpha_prime and alpha_dual) during each iteration, allowing for more efficient convergence. The algorithm follows these steps:

Initialization: Initialize the decision variables, slack variables, and Lagrange multipliers.

Iterative Process: Iterate until convergence or a maximum number of iterations is reached.

a. Compute the residuals (rc, rb, rxs) and the affine search direction (deltax_affine, deltas_affine).

b. Compute the step lengths (alpha_prime and alpha_dual) using the affine search direction.

c. Perform the corrector step to update the decision variables.

d. Update the step lengths (alpha_prime and alpha_dual) using a line search technique.

Termination: Stop the iterations when the termination criterion is met (e.g., tolerance on the central path measure).

# Mehrotra Predictor

The Mehrotra Predictor algorithm is another approach to solving linear programming problems. It combines a predictor-corrector framework with a self-adjusting parameter (sigma) to improve convergence. The algorithm follows these steps:

Initialization: Initialize the decision variables, slack variables, and Lagrange multipliers.

Iterative Process: Iterate until convergence or a maximum number of iterations is reached.

a. Compute the residuals (rc, rb, rxs) and the affine search direction (deltax_affine, deltas_affine).

b. Compute the step lengths (alpha_prime and alpha_dual) using the affine search direction.

c. Perform the corrector step to update the decision variables.

d. Update the self-adjusting parameter (sigma) based on the predictor-corrector framework.

Termination: Stop the iterations when the termination criterion is met (e.g., tolerance on the central path measure).

# Usage
You can find usage example in Test_R2.m
