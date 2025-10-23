using UUIDs

struct DDARunner
    binary_path::String

    function DDARunner(binary_path::String)
        if !isfile(binary_path)
            throw(BinaryNotFoundError(binary_path))
        end
        new(binary_path)
    end
end

function run(
    runner::DDARunner,
    request::DDARequest,
    start_bound::Union{Int,Nothing}=nothing,
    end_bound::Union{Int,Nothing}=nothing,
    edf_channel_names::Union{Vector{String},Nothing}=nothing
)::DDAResult
    analysis_id = string(uuid4())

    if !isfile(request.file_path)
        throw(FileNotFoundError(request.file_path))
    end

    @info "Starting DDA analysis for file: $(request.file_path)"
    @info "Channel indices (0-based from frontend): $(request.channels)"
    @info "Time range: $(request.time_range)"
    @info "Window parameters: $(request.window_parameters)"
    @info "Scale parameters: $(request.scale_parameters)"

    temp_dir = tempdir()
    output_file = joinpath(temp_dir, "dda_output_$(analysis_id).txt")

    channel_indices = if !isnothing(request.channels)
        [string(idx + 1) for idx in request.channels]
    else
        ["1"]
    end

    @info "Channel indices (1-based for DDA binary): $channel_indices"

    is_windows = Sys.iswindows()

    select_mask = if !isnothing(request.algorithm_selection.select_mask)
        request.algorithm_selection.select_mask
    else
        "1 0 0 0"
    end

    select_bits = split(select_mask)
    st_enabled = length(select_bits) >= 1 && select_bits[1] == "1"
    ct_enabled = length(select_bits) >= 2 && select_bits[2] == "1"
    has_ct_pairs = ct_enabled && !isnothing(request.ct_channel_pairs) && !isempty(request.ct_channel_pairs)

    run_st_separately = st_enabled && has_ct_pairs

    first_select_mask = if run_st_separately
        bits = copy(select_bits)
        if length(bits) > 1
            bits[2] = "0"
        end
        join(bits, " ")
    else
        select_mask
    end

    @info "First execution SELECT mask: $first_select_mask"
    if run_st_separately
        @info "Will run CT separately for $(length(request.ct_channel_pairs)) channel pairs"
    end

    file_ext = lowercase(splitext(request.file_path)[2])
    is_ascii_file = file_ext in [".ascii", ".txt"]
    file_type_flag = is_ascii_file ? "-ASCII" : "-EDF"

    @info "Using file type flag: $file_type_flag"

    actual_input_file = request.file_path
    if is_ascii_file
        temp_ascii_file = joinpath(temp_dir, "dda_input_$(analysis_id).ascii")

        content = read(request.file_path, String)
        lines = split(content, '\n')

        if isempty(lines)
            throw(ParseError("ASCII file is empty"))
        end

        first_line = lines[1]
        has_header = any(isalpha, first_line)

        data_lines = has_header ? lines[2:end] : lines

        @info "Detected header in ASCII file, stripping it for DDA binary"

        data_content = join(data_lines, '\n')
        write(temp_ascii_file, data_content)

        @info "Created temporary headerless ASCII file: $temp_ascii_file"
        actual_input_file = temp_ascii_file
    end

    cmd_args = String[]

    if !is_windows
        push!(cmd_args, "sh")
        push!(cmd_args, runner.binary_path)
    end

    append!(cmd_args, [
        "-DATA_FN", actual_input_file,
        "-OUT_FN", output_file,
        file_type_flag,
        "-CH_list"
    ])

    if run_st_separately
        append!(cmd_args, channel_indices)
        @info "Using $(length(channel_indices)) individual channels for ST variant"
    elseif ct_enabled && has_ct_pairs
        first_pair = request.ct_channel_pairs[1]
        push!(cmd_args, string(first_pair[1] + 1))
        push!(cmd_args, string(first_pair[2] + 1))
        @info "CT only - processing first pair (1-based): [$(first_pair[1] + 1), $(first_pair[2] + 1)]"
        if length(request.ct_channel_pairs) > 1
            @info "Will process $(length(request.ct_channel_pairs) - 1) additional CT pairs in separate executions"
        end
    else
        append!(cmd_args, channel_indices)
        @info "Using $(length(channel_indices)) individual channels for enabled variants"
    end

    append!(cmd_args, [
        "-dm", "4",
        "-order", "4",
        "-nr_tau", "2",
        "-WL", string(request.window_parameters.window_length),
        "-WS", string(request.window_parameters.window_step),
        "-SELECT"
    ])

    append!(cmd_args, split(first_select_mask))

    append!(cmd_args, ["-MODEL", "1", "2", "10"])

    if !isnothing(request.window_parameters.ct_window_length)
        append!(cmd_args, ["-WL_CT", string(request.window_parameters.ct_window_length)])
    end
    if !isnothing(request.window_parameters.ct_window_step)
        append!(cmd_args, ["-WS_CT", string(request.window_parameters.ct_window_step)])
    end

    delay_min = Int(request.scale_parameters.scale_min)
    delay_max = Int(request.scale_parameters.scale_max)
    push!(cmd_args, "-TAU")
    for delay in delay_min:delay_max
        push!(cmd_args, string(delay))
    end

    if !isnothing(start_bound) && !isnothing(end_bound)
        append!(cmd_args, ["-StartEnd", string(start_bound), string(end_bound)])
    end

    @info "Executing DDA command: $(join(cmd_args, ' '))"

    start_time = time()

    cmd = is_windows ? Cmd([runner.binary_path; cmd_args[5:end]]) : Cmd(cmd_args)
    try
        Base.run(cmd)
    catch e
        throw(ExecutionFailedError("Failed to execute DDA binary: $e"))
    end

    @info "DDA binary execution completed in $(round(time() - start_time, digits=2))s"

    variant_names = ["ST", "CT", "CD", "DE"]
    enabled_variants = [variant_names[i] for (i, bit) in enumerate(select_bits) if bit == "1"]

    if isempty(enabled_variants)
        throw(ExecutionFailedError("No variants enabled in SELECT mask"))
    end

    @info "Enabled variants (from original request): $enabled_variants"

    first_execution_variants = run_st_separately ? ["ST"] : enabled_variants

    output_file_stem = splitext(basename(output_file))[1]
    output_dir = dirname(output_file)

    variant_matrices = Tuple{String,Matrix{Float64}}[]

    for variant in first_execution_variants
        variant_file_path = joinpath(output_dir, "$(output_file_stem)_$(variant)")
        variant_file_with_ext = "$(output_file)_$(variant)"

        actual_output_file = if isfile(variant_file_path)
            variant_file_path
        elseif isfile(variant_file_with_ext)
            variant_file_with_ext
        else
            @warn "Output file not found for variant $variant. Tried: $variant_file_path, $variant_file_with_ext"
            continue
        end

        @info "Reading DDA output for variant $variant from: $actual_output_file"

        output_content = read(actual_output_file, String)

        @info "Output file size for $variant: $(length(output_content)) bytes"

        q_matrix = parse_dda_output(output_content)

        if !isempty(q_matrix)
            num_channels = size(q_matrix, 1)
            num_timepoints = size(q_matrix, 2)
            @info "Q matrix dimensions for $variant: $num_channels channels × $num_timepoints timepoints"

            push!(variant_matrices, (variant, q_matrix))
        else
            @warn "No data extracted from DDA output for variant $variant"
        end
    end

    if run_st_separately && has_ct_pairs
        @info "Now processing CT variant with $(length(request.ct_channel_pairs)) channel pairs"

        pairs = request.ct_channel_pairs
        combined_ct_matrix = Matrix{Float64}(undef, 0, 0)
        first_pair_processed = false

        for (pair_idx, pair) in enumerate(pairs)
            pair_output_file = joinpath(temp_dir, "dda_output_$(analysis_id)_pair$(pair_idx-1).txt")

            pair_cmd_args = String[]

            if !is_windows
                push!(pair_cmd_args, "sh")
                push!(pair_cmd_args, runner.binary_path)
            end

            append!(pair_cmd_args, [
                "-DATA_FN", actual_input_file,
                "-OUT_FN", pair_output_file,
                file_type_flag,
                "-CH_list",
                string(pair[1] + 1),
                string(pair[2] + 1),
                "-dm", "4",
                "-order", "4",
                "-nr_tau", "2",
                "-WL", string(request.window_parameters.window_length),
                "-WS", string(request.window_parameters.window_step),
                "-SELECT", "0", "1", "0", "0",
                "-MODEL", "1", "2", "10"
            ])

            if !isnothing(request.window_parameters.ct_window_length)
                append!(pair_cmd_args, ["-WL_CT", string(request.window_parameters.ct_window_length)])
            end
            if !isnothing(request.window_parameters.ct_window_step)
                append!(pair_cmd_args, ["-WS_CT", string(request.window_parameters.ct_window_step)])
            end

            push!(pair_cmd_args, "-TAU")
            for delay in delay_min:delay_max
                push!(pair_cmd_args, string(delay))
            end

            if !isnothing(start_bound) && !isnothing(end_bound)
                append!(pair_cmd_args, ["-StartEnd", string(start_bound), string(end_bound)])
            end

            @info "Executing DDA for CT pair $(pair_idx-1) (1-based): [$(pair[1] + 1), $(pair[2] + 1)]"

            pair_cmd = is_windows ? Cmd([runner.binary_path; pair_cmd_args[5:end]]) : Cmd(pair_cmd_args)
            try
                Base.run(pair_cmd)
            catch e
                @error "DDA binary failed for CT pair $(pair_idx-1): $e"
                continue
            end

            pair_output_file_stem = splitext(basename(pair_output_file))[1]
            pair_ct_file_path = joinpath(output_dir, "$(pair_output_file_stem)__CT")
            pair_ct_file_with_ext = "$(pair_output_file)_CT"

            actual_pair_ct_file = if isfile(pair_ct_file_path)
                pair_ct_file_path
            elseif isfile(pair_ct_file_with_ext)
                pair_ct_file_with_ext
            else
                @warn "CT output file not found for pair $(pair_idx-1). Skipping."
                continue
            end

            pair_content = read(actual_pair_ct_file, String)
            pair_q_matrix = parse_dda_output(pair_content)

            if !isempty(pair_q_matrix)
                @info "Adding CT pair $(pair_idx-1) results: $(size(pair_q_matrix, 1)) channels × $(size(pair_q_matrix, 2)) timepoints"

                if !first_pair_processed
                    combined_ct_matrix = pair_q_matrix
                    first_pair_processed = true
                else
                    combined_ct_matrix = vcat(combined_ct_matrix, pair_q_matrix)
                end
            end

            try
                rm(actual_pair_ct_file)
            catch
            end
            try
                rm(pair_output_file)
            catch
            end
        end

        if !isempty(combined_ct_matrix)
            num_channels = size(combined_ct_matrix, 1)
            num_timepoints = size(combined_ct_matrix, 2)
            @info "Combined CT Q matrix dimensions: $num_channels channels × $num_timepoints timepoints"
            push!(variant_matrices, ("CT", combined_ct_matrix))
        else
            @warn "No CT data extracted from any pair"
        end
    end

    try
        rm(output_file)
    catch
    end
    for variant in enabled_variants
        variant_path = joinpath(output_dir, "$(output_file_stem)_$(variant)")
        try
            rm(variant_path)
        catch
        end
    end

    if is_ascii_file && actual_input_file != request.file_path
        try
            rm(actual_input_file)
        catch
        end
    end

    if isempty(variant_matrices)
        throw(ParseError("No data extracted from any DDA variant"))
    end

    primary_variant_name, primary_q_matrix = variant_matrices[1]

    @info "Using $primary_variant_name as primary variant, $(length(variant_matrices)) total variants processed"

    channels = if !isnothing(request.channels)
        ["Channel $(idx + 1)" for idx in request.channels]
    else
        ["Channel 1"]
    end

    function variant_display_name(id::String)
        if id == "ST"
            "Single Timeseries (ST)"
        elseif id == "CT"
            "Cross-Timeseries (CT)"
        elseif id == "CD"
            "Cross-Delay (CD)"
        elseif id == "DE"
            "Delay Evolution (DE)"
        else
            id
        end
    end

    variant_results = [
        VariantResult(
            variant_id,
            variant_display_name(variant_id),
            q_matrix,
            if variant_id == "CT" && !isnothing(request.ct_channel_pairs)
                pairs = request.ct_channel_pairs
                if !isnothing(edf_channel_names)
                    ["$(edf_channel_names[pair[1] + 1]) ⟷ $(edf_channel_names[pair[2] + 1])" for pair in pairs]
                else
                    ["Ch$(pair[1] + 1) ⟷ Ch$(pair[2] + 1)" for pair in pairs]
                end
            elseif !isnothing(edf_channel_names)
                channel_indices = request.channels
                if !isnothing(channel_indices)
                    [edf_channel_names[idx+1] for idx in channel_indices if idx + 1 <= length(edf_channel_names)]
                else
                    edf_channel_names
                end
            else
                nothing
            end
        )
        for (variant_id, q_matrix) in variant_matrices
    ]

    result = DDAResult(
        analysis_id,
        request.file_path,
        channels,
        primary_q_matrix,
        variant_results,
        nothing,
        request.window_parameters,
        request.scale_parameters,
        now()
    )

    return result
end

binary_path(runner::DDARunner) = runner.binary_path
