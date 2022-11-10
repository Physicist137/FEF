# FEF

This program is able to calculate the Field Enhancement Factor (FEF) of a general nano-structure, using classical electrostatics. It solves laplace equation by calculating all multipole moments until truncation, building a certain kind of multipole expansion in the process. It does so by numerically calculating integrals that came from minimization of a variational responsible for enforcing dirichlet boundary conditions, and saving their values in a database, and then by accessing them to extrapolate the multipole moments.

This code was made in 2020.
