"""
$(TYPEDEF)

Base type for all CSE problems.
"""
abstract type CSEProblem end

"""
$(TYPEDEF)

Base type for all symmetric CSE problems.
"""
abstract type SymmetricCSEProblem <: CSEProblem end

"""
$(TYPEDEF)

Base type for all asymmetric CSE problems.
"""
abstract type AsymmetricCSEProblem <: CSEProblem end

"""
$(TYPEDEF)

Base type for all solver parameter objects passed to objective functions.
"""
abstract type CSESolverParams end

"""
$(TYPEDEF)

Base type for all CSE solutions.
"""
abstract type CSESolution end
