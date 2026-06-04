# =============================================================================
# COMMAND BUILDING
# =============================================================================

"""Build the DDA CLI command."""
function build_command(runner::DDARunner, request::DDARequest, output_base::String)::Cmd
    args = String[]

    # File type flag
    ext = lowercase(splitext(request.file_path)[2])
    push!(args, ext == ".edf" ? "-EDF" : "-ASCII")

    # Input/output
    push!(args, "-DATA_FN", request.file_path)
    push!(args, "-OUT_FN", output_base)

    # Channels are already 1-based for Julia callers and for the binary.
    push!(args, "-CH_list")
    for ch in request.channels
        push!(args, string(ch))
    end

    # SELECT mask
    mask = request.select === nothing ? generate_select_mask(request.variants) : request.select
    push!(args, "-SELECT")
    for bit in mask
        push!(args, string(bit))
    end

    # Model encoding
    push!(args, "-MODEL")
    for p in request.model_terms
        push!(args, string(p))
    end

    # Delay values
    push!(args, "-TAU")
    for d in request.delay_params.delays
        push!(args, string(d))
    end

    # Window parameters. `-WLms` and `-WSms` are special binary flags and are
    # intentionally not emitted by this wrapper.
    wp = request.window_params
    if wp.WL !== nothing
        push!(args, "-WL", string(wp.WL))
    end
    if wp.WS !== nothing
        push!(args, "-WS", string(wp.WS))
    end

    # Model parameters
    mp = request.model_params
    push!(args, "-dm", string(mp.derivative_points))
    push!(args, "-order", string(mp.order))
    push!(args, "-nr_tau", string(mp.nr_tau))

    if wp.ct_window_length !== nothing
        push!(args, "-WL_CT", string(wp.ct_window_length))
    end
    if wp.ct_window_step !== nothing
        push!(args, "-WS_CT", string(wp.ct_window_step))
    end

    # Time bounds
    if request.time_range !== nothing
        tr = request.time_range
        push!(args, "-StartEnd", string(Int(tr.start)), string(Int(tr.stop)))
    end

    sampling_rate_args = _sampling_rate_args(request.sampling_rate)
    if !isempty(sampling_rate_args)
        push!(args, "-SR")
        append!(args, sampling_rate_args)
    end

    append!(args, request.passthrough_args)

    if REQUIRES_SHELL_WRAPPER && !Sys.iswindows()
        return `sh $(runner.binary_path) $args`
    else
        return `$(runner.binary_path) $args`
    end
end

"""Render the generated shell call as a plain string for diagnostics and tests."""
function build_command_string(
    runner::DDARunner,
    request::DDARequest,
    output_base::String,
)::String
    return join(collect(build_command(runner, request, output_base)), " ")
end
