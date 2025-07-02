# # Symmetric CSE example
#
# This is a simple example of how to run the symmetric CSE case.
#
# ## Install dependencies
#
# Make sure the required packages are installed by running the following from
# the root directory of this repo:

# ```julia
# using Pkg
# Pkg.instantiate()
# ```

# ## Load required modules

using ConstrainedStrategicEquilibrium
using Plots
using Distributions
using NonlinearSolve

# ## Create the symmetric CSE problem
#
# The original fortran asym code used n=16 and the following initial guess:

nval = 16
xguess = Vector{Float64}(undef, 2 * nval - 1)
xguess[1] = -2.84827109173688
xguess[2] = -2.85209052797506
xguess[3] = -2.85688597597939
xguess[4] = -2.86243170560486
xguess[5] = -2.87005433886236
xguess[6] = -2.87972857171734
xguess[7] = -2.89195095174828
xguess[8] = -2.90986034638539
xguess[9] = -2.93411212816428
xguess[10] = -2.96990679761558
xguess[11] = -3.03212449908756
xguess[12] = -3.14051739544213
xguess[13] = -3.33100693581025
xguess[14] = -3.69184084490141
xguess[15] = -5.61683383184477
xguess[16] = -7.41860036191254
xguess[17] = 4.37486841314759
xguess[18] = 3.55315379806175
xguess[19] = 3.09860286531008
xguess[20] = 2.74440338704673
xguess[21] = 2.45415125926837
xguess[22] = 2.19413505997764
xguess[23] = 1.94448901027584
xguess[24] = 1.71125414593309
xguess[25] = 1.47745377147087
xguess[26] = 1.23478257742541
xguess[27] = 1.01951187118377
xguess[28] = 0.846654741348483
xguess[29] = 0.726435662568677
xguess[30] = 0.619443991395723
xguess[31] = 1.12697386952317

#
# Create an asymmetric CSE problem with the following:
cse_prob = AsymmetricAfrprogsCSEProblem(
    inin=nval,
    maxn=nval,
    np=4,
    mc=10000,
    solver=FastShortcutNonlinearPolyalg(; autodiff=AutoFiniteDiff()),
    solver_kwargs=(; show_trace=Val(true), trace_level=TraceMinimal()),
    solver_initial_guess=xguess,
)

# ## Compute the CSE
#
# Now compute the CSE for the problem we created
solutions = compute_cse(cse_prob)

# ## Postprocessing
#
# Select the final solution from the list
sol = solutions[end]

# Plot the final solution comparing the computed CSE to the "analytical" Bayes-Nash Equilibrium.
plot(sol, dpi=300)

# We can also save the figure to a file
savefig("afr_progs_asym.png")
