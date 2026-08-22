"""Lazy loading of optional heavy dependencies (Plots, DataFrames, CUDA)."""
module OptionalDeps

const _cache = Dict{Symbol, Module}()

"""
    OptionalDeps.require(:Name) -> Module

Load an optional package on first use and cache the module. Throws a helpful
error when the package is not installed.
"""
function require(name::Symbol)::Module
    return get!(_cache, name) do
        try
            Base.require(Main, name)
        catch
            error(
                "Package $(String(name)) is required for this functionality. " *
                "Install it with: using Pkg; Pkg.add(\"$(String(name))\")",
            )
        end
    end::Module
end

end # module OptionalDeps
