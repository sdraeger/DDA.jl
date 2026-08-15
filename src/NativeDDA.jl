"""Pure Julia DDA engine with CPU and optional CUDA backends."""
module NativeDDA

using LinearAlgebra
using Statistics
using ..ModelEncoding: generate_monomials

export NativeDDAResult, NativeFlavorResult, flavor_result, run_dda_matrix

include("NativeDDA/Types.jl")
include("NativeDDA/Preparation.jl")
include("NativeDDA/Problems.jl")
include("NativeDDA/CPU.jl")
include("NativeDDA/GPU.jl")
include("NativeDDA/Engine.jl")

end
