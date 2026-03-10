"""
COMPADRE/COMADRE database types for Julia.
"""

using DataFrames

struct CompadreMat{T<:Real}
    matA::Matrix{T}
    matU::Matrix{T}
    matF::Matrix{T}
    matC::Matrix{T}
end

struct CompadreDB <: MPM.AbstractCompadreDB
    data::DataFrame
    version::String
    db_type::Symbol  # :compadre or :comadre
end
