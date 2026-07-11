using Test
using DelayDifferentialAnalysis

@testset "Claudia compatibility contract" begin
    binary_path = tempname()
    touch(binary_path)

    try
        runner = DDARunner(binary_path)
        request = DDARequest(
            "recording.custom",
            1:4,
            ["ST", "CT"];
            input_format=:ascii,
            model=[1, 2, 6],
            derivative_points=4,
            order=3,
            delays=[32, 9],
            WL=3000,
            WS=200,
            WL_CT=2,
            WS_CT=2,
        )
        command = DelayDifferentialAnalysis.Runner._logical_command_parts(
            runner,
            request,
            "RUN1/structure_selection_01_02_06",
        )

        @test command == [
            binary_path,
            "-ASCII",
            "-DATA_FN", "recording.custom",
            "-OUT_FN", "RUN1/structure_selection_01_02_06",
            "-CH_list", "1", "2", "3", "4",
            "-SELECT", "1", "1", "0", "0", "0", "0",
            "-MODEL", "1", "2", "6",
            "-TAU", "32", "9",
            "-WL", "3000",
            "-WS", "200",
            "-dm", "4",
            "-order", "3",
            "-nr_tau", "2",
            "-WL_CT", "2",
            "-WS_CT", "2",
        ]

        no_channels = DDARequest(
            "recording.custom",
            nothing,
            ["ST"];
            input_format="ascii",
        )
        no_channels_command = collect(
            DelayDifferentialAnalysis.Runner.build_command(
                runner,
                no_channels,
                "output",
            ),
        )
        @test !("-CH_list" in no_channels_command)

        MOD = make_MOD(3, 3)
        model = findall(==(1), vec(MOD[1, :]))
        P_DDA = DelayDifferentialAnalysis.StructureSelection._p_dda(3)
        model_id = DelayDifferentialAnalysis.StructureSelection._model_filename_id(
            model,
            P_DDA;
            nr_delays=2,
            order=3,
        )
        @test model_id == join(lpad.(string.(model), 2, '0'), "_")
        @test DelayDifferentialAnalysis.StructureSelection._trial_out_fn(
            "RUN1",
            model_id,
            nothing,
        ) == joinpath("RUN1", "structure_selection_$model_id")
    finally
        rm(binary_path; force=true)
    end
end
