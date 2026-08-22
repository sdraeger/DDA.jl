"""High-level analysis functions for the DDA binary."""
module API

using Printf
using ..DDADefaults
using ..Results: STResult, CTResult, DEResult
using ..Runner
using ..Runner: StructuredChannelData, run_analysis_structured

include("API/Helpers.jl")
include("API/Conversions.jl")
include("API/FileAPI.jl")
include("API/MatrixAPI.jl")
include("API/PublicAPI.jl")

end # module API
