module MatrixProjectionModelsCompadreExt

import MatrixProjectionModels as MPM
import CSV
import DataFrames
import Downloads
import RData

include("CompadreExt/types.jl")
include("CompadreExt/io.jl")
include("CompadreExt/accessors.jl")
include("CompadreExt/flags.jl")
include("CompadreExt/operations.jl")
include("CompadreExt/statistics.jl")
include("CompadreExt/stage_queries.jl")
include("CompadreExt/build.jl")

end # module
