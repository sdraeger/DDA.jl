"""High-level analysis functions for the DDA binary."""
module API

using Printf
using ..DDADefaults
using ..Results
using ..Runner

export run_st, run_ct, run_de

include("API/Helpers.jl")
include("API/Conversions.jl")
include("API/FileAPI.jl")
include("API/MatrixAPI.jl")
include("API/PublicAPI.jl")

end # module API
