"""
Supplement and overwrite specific matrix transitions.

Allows injecting known values (from literature, prior knowledge) into
specific matrix elements, either as fixed rates or as multipliers of
estimated values.

Reference: lefko3's supplemental() function (Shefferson & Ehrlen).
"""

"""
    TransitionSupplement

A single transition supplement specification.

# Fields
- `from`: Source stage index or name
- `to`: Target stage index or name
- `value`: The value to apply
- `type`: `:overwrite` (replace element), `:multiplier` (multiply existing),
          or `:add` (add to existing)
"""
struct TransitionSupplement{V<:Real}
    from::Int
    to::Int
    value::V
    type::Symbol
end

"""
    supplement!(A::AbstractMatrix, specs::AbstractVector{TransitionSupplement})

Apply transition supplements to a projection matrix in-place.

# Arguments
- `A`: Projection matrix to modify
- `specs`: Vector of `TransitionSupplement` specifications

# Returns
The modified matrix `A`.
"""
function supplement!(A::AbstractMatrix, specs::AbstractVector{<:TransitionSupplement})
    for s in specs
        if s.type == :overwrite
            A[s.to, s.from] = s.value
        elseif s.type == :multiplier
            A[s.to, s.from] *= s.value
        elseif s.type == :add
            A[s.to, s.from] += s.value
        else
            throw(ArgumentError("Unknown supplement type: $(s.type). Use :overwrite, :multiplier, or :add"))
        end
    end
    return A
end

"""
    supplement!(A::AbstractMatrix; overwrites=Pair[], multipliers=Pair[], additions=Pair[])

Apply supplements using keyword syntax.

# Example
```julia
A = zeros(3, 3)
supplement!(A;
    overwrites = [(2, 1) => 0.8, (3, 2) => 0.5],  # (to, from) => value
    multipliers = [(1, 3) => 1.2],
)
```
"""
function supplement!(A::AbstractMatrix;
                     overwrites::AbstractVector{<:Pair}=Pair{Tuple{Int,Int},Float64}[],
                     multipliers::AbstractVector{<:Pair}=Pair{Tuple{Int,Int},Float64}[],
                     additions::AbstractVector{<:Pair}=Pair{Tuple{Int,Int},Float64}[])
    specs = TransitionSupplement[]
    for ((to, from), val) in overwrites
        push!(specs, TransitionSupplement(from, to, Float64(val), :overwrite))
    end
    for ((to, from), val) in multipliers
        push!(specs, TransitionSupplement(from, to, Float64(val), :multiplier))
    end
    for ((to, from), val) in additions
        push!(specs, TransitionSupplement(from, to, Float64(val), :add))
    end
    return supplement!(A, specs)
end

"""
    supplement(A::AbstractMatrix; kwargs...)

Non-mutating version of `supplement!`. Returns a modified copy.
"""
function supplement(A::AbstractMatrix; kwargs...)
    B = copy(float.(A))
    supplement!(B; kwargs...)
    return B
end
