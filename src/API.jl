"""High-level analysis functions for the DDA binary."""
module API

using Printf
using ..DDADefaults
using ..Results: STResult, CTResult, DEResult
using ..Runner
using ..Runner: VariantResultData, DDARunner, DDARequest

include("API/Helpers.jl")
include("API/Conversions.jl")
include("API/FileAPI.jl")
include("API/MatrixAPI.jl")
include("API/PublicAPI.jl")

end # module API
