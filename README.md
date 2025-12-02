# DelayDifferentialAnalysis.jl

Julia package for Delay Differential Analysis (DDA).

## Installation

```julia
using Pkg
Pkg.add("DelayDifferentialAnalysis")
```

Or, for the development version:

```julia
using Pkg
Pkg.add(url="https://github.com/sdraeger/DelayDifferentialAnalysis.jl")
```

## Usage

```julia
using DelayDifferentialAnalysis

# Access variant metadata
println(ST.full_name)  # "Standard"
println(ST.abbreviation)  # "ST"

# Generate select masks
mask = generate_select_mask([ST, CT])
println(mask)  # "110000"

# Parse existing masks
variants = parse_select_mask("110000")
println(variants)  # [ST, CT]
```

## Variants

The package provides access to all DDA variants:

- `ST` - Standard
- `CT` - Continuous Time
- `CD` - Continuous Data
- `DE` - Differential Equation
- `SY` - Synchrony

## License

MIT
