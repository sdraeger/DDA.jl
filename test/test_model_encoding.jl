using Test
using DelayDifferentialAnalysis

@testset "ModelEncoding" begin
    @testset "generate_monomials basic" begin
        monos = generate_monomials(2, 2)
        # For num_delays=2, order=2: (0,1), (0,2), (1,1), (1,2), (2,2)
        @test length(monos) == 5
        @test monos[1] == [0, 1]
        @test monos[2] == [0, 2]
        @test monos[3] == [1, 1]
        @test monos[4] == [1, 2]
        @test monos[5] == [2, 2]
    end

    @testset "generate_monomials order=1" begin
        monos = generate_monomials(3, 1)
        @test length(monos) == 3
        @test monos == [[1], [2], [3]]
    end

    @testset "generate_monomials standard DDA (2, 4)" begin
        monos = generate_monomials(2, 4)
        @test length(monos) == 14
        # First monomials should be linear terms
        @test monos[1] == [0, 0, 0, 1]
        @test monos[2] == [0, 0, 0, 2]
    end

    @testset "monomial_to_text" begin
        @test monomial_to_text([0, 1]) == "x_1"
        @test monomial_to_text([0, 2]) == "x_2"
        @test monomial_to_text([1, 1]) == "x_1^2"
        @test monomial_to_text([1, 2]) == "x_1 * x_2"
        @test monomial_to_text([2, 2]) == "x_2^2"
    end

    @testset "monomial_to_text with tau_values" begin
        @test monomial_to_text([0, 1]; tau_values=[7, 10]) == "x(t-7)"
        @test monomial_to_text([1, 2]; tau_values=[7, 10]) == "x(t-7) * x(t-10)"
    end

    @testset "monomial_to_latex" begin
        @test monomial_to_latex([0, 1]) == "x_{1}"
        @test monomial_to_latex([1, 1]) == "x_{1}^{2}"
        @test monomial_to_latex([1, 2]) == "x_{1} \\cdot x_{2}"
    end

    @testset "decode_model_encoding text" begin
        eq = decode_model_encoding([1, 2, 10]; num_delays=2, polynomial_order=4, format="text")
        @test startswith(eq, "dx/dt = ")
        @test occursin("a_1", eq)
        @test occursin("a_2", eq)
        @test occursin("a_3", eq)
    end

    @testset "decode_model_encoding latex" begin
        eq = decode_model_encoding([1, 2]; num_delays=2, polynomial_order=2, format="latex")
        @test startswith(eq, "\\frac{dx}{dt}")
        @test occursin("a_{1}", eq)
    end

    @testset "decode_model_encoding out of range" begin
        @test_throws ErrorException decode_model_encoding([999]; num_delays=2, polynomial_order=2)
    end

    @testset "model matrix rows map to MODEL indices" begin
        pp = [
            0 0 1
            0 0 2
            1 1 1
        ]
        @test model_matrix_to_encoding(pp; num_delays=2, polynomial_order=3) == [1, 2, 6]

        order4 = [
            0 0 0 1
            0 0 0 2
            1 1 1 1
        ]
        @test model_matrix_to_encoding(order4; num_delays=2, polynomial_order=4) == [1, 2, 10]

        @test_throws ErrorException model_matrix_to_encoding(pp; num_delays=2, polynomial_order=4)
        @test_throws ErrorException model_matrix_to_encoding([0 2 1]; num_delays=2, polynomial_order=3)
        @test_throws ErrorException model_matrix_to_encoding([0 0 3]; num_delays=2, polynomial_order=3)
        @test_throws ErrorException model_matrix_to_encoding([0 0 0]; num_delays=2, polynomial_order=3)
    end

    @testset "visualize_model_space" begin
        output = visualize_model_space(2, 2)
        @test occursin("Model Space:", output)
        @test occursin("Total monomials: 5", output)
        @test occursin("x_1", output)
    end

    @testset "visualize_model_space with highlighting" begin
        output = visualize_model_space(2, 2; highlight_encoding=[1, 2])
        @test occursin("*", output)
        @test occursin("Selected terms", output)
        @test occursin("dx/dt", output)
    end
end
