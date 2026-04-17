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

"""
    to_discrete_problem(prob::CoupledMPMProblem)

Convert a coupled/stateful MPM problem to a `SciMLBase.DiscreteProblem`.
The state is a `PopulationSystem`, so callbacks can inspect named components
and auxiliary state directly.
"""
function to_discrete_problem(prob::CoupledMPMProblem)
    function coupled_step(sys, p, t)
        sys_new = _copy_system(sys)
        _advance_one_day!(sys_new, Int(round(t)), p, prob.events, prob.substeps, prob.rules;
            normalize = prob.normalize)
        return sys_new
    end
    return SciMLBase.DiscreteProblem(coupled_step, _copy_system(prob.system), Float64.(prob.tspan), prob.p)
end
