"""
SciML interface: convert MPMProblem to DiscreteProblem.
"""

"""
    to_discrete_problem(prob::MPMProblem)

Convert an MPMProblem to a SciMLBase.DiscreteProblem for use with
DifferenceEquations.jl or other SciML solvers.
"""
function to_discrete_problem(prob::MPMProblem)
    A = _get_matrix(prob)
    f = (u, p, t) -> A * u
    return SciMLBase.DiscreteProblem(f, float.(prob.n0), prob.tspan, prob.p)
end
