"""
Plot recipes for MatrixProjectionModels using RecipesBase.
"""

# Matrix heatmap
@recipe function f(mpm::MatrixProjectionModel)
    seriestype := :heatmap
    yflip := true
    title --> "Projection Matrix"
    xlabel --> "From stage"
    ylabel --> "To stage"
    n = n_stages(mpm)
    names = string.(mpm.stage_names)
    xticks --> (1:n, names)
    yticks --> (1:n, names)
    1:n, 1:n, mpm.A
end

# Population trajectory
@recipe function f(sol::MPMSolution)
    xlabel --> "Time"
    ylabel --> "Population size"
    title --> "Population Trajectory"
    label --> "Total N"
    sol.t, [sum(u) for u in sol.u]
end
