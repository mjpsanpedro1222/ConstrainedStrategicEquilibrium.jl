# # Uniform distribution
#
# **Note:** this example is not working
#
# This example is using the Uniform distribution with the asymmetric problem
# but all players have the same distribution.
#
# ## Install dependencies
#
# Enter the Julia REPL (`julia`) and run:
#
# ```julia
# using Pkg
# Pkg.add("ConstrainedStrategicEquilibrium")
# ```

# ## Load required modules

using ConstrainedStrategicEquilibrium
using Plots
using Distributions
using NonlinearSolve

# ## Create the CSE problem
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
# Create an asymmetric CSE problem with the following (not we make a small change
# from the Fortran asym.f code by running from n=16 to n=17, instead of just n=16):
cse_prob = AsymmetricAfrprogsCSEProblem(
    inin=nval,
    maxn=nval,
    np=2,
    mc=10000,
    solver_kwargs=(; show_trace=Val(true)),
    solver_initial_guess=xguess,
    distributions=[Uniform(0, 1), Uniform(0, 1)],
)

# Now validate the problem we just created (this will throw an error if there is a mistake):
validate_cse_problem(cse_prob)

# ## Compute the CSE
#
# Now compute the CSE for the problem we created
solutions = compute_cse(cse_prob)

# ## Postprocessing
#
# Loop over the solutions and plot them all:
for sol in solutions
    plot(sol, dpi=300)

    # We can also save the figure to a file (comment this line out if you don't want to)
    savefig("afr-progs-asym-n$(sol.n).png")
end
