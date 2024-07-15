using Test
using TwinCopulas
using Distributions
using Random
@testset "Extreme Value Copulas Tests" begin
    rng = MersenneTwister(2024)
    d = 100
    @testset "AsymGalambosCopula - sampling, pdf, cdf" begin
        for α in [0.0, rand(rng, Uniform()), rand(rng, Uniform(5.0, 9.0)), rand(rng, Uniform(10.0, 15.0))]
            for θ in [[rand(rng, Uniform(0, 1)), rand(rng, Uniform(0, 1))], [0.0, 0.0], [1.0, 1.0]]
                C = AsymGalambosCopula(α, θ)
                data = rand(rng, C, 100)
                
                # Verificar que la generación de números aleatorios tiene la forma correcta
                @test size(data) == (2, 100)
    
                # Evaluar CDF y PDF en cada columna de los datos generados
                for i in 1:d
                    u = data[:,i]
                    cdf_value = cdf(C, u)
                    pdf_value = pdf(C, u)
                    @test 0.0 <= cdf_value <= 1.0
                    @test pdf_value >= 0.0
                end    
            end
        end
    end

    @testset "AsymLogCopula - sampling, pdf, cdf" begin
        for α in [1.0, rand(rng, Uniform(1.0, 5.0)), rand(rng, Uniform(10.0, 15.0))]
            for θ in [[rand(rng, Uniform(0, 1)), rand(rng, Uniform(0, 1))], [0.0, 0.0], [1.0, 1.0]]
                C = AsymLogCopula(α, θ)
                data = rand(rng, C, 100)
                
                # Verificar que la generación de números aleatorios tiene la forma correcta
                @test size(data) == (2, 100)
    
                # Evaluar CDF y PDF en cada columna de los datos generados
                for i in 1:d
                    u = data[:,i]
                    cdf_value = cdf(C, u)
                    pdf_value = pdf(C, u)
                    @test 0.0 <= cdf_value <= 1.0
                    @test pdf_value >= 0.0
                end    
            end
        end
    end

    @testset "AsymMixedCopula - sampling, pdf, cdf" begin
        for θ in [[0.1, 0.2], [0.0, 0.0], [0.3, 0.4], [0.2, 0.4]]
            try
                C = AsymMixedCopula(θ)
                data = rand(rng, C, 100)
    
                # Verificar que la generación de números aleatorios tiene la forma correcta
                @test size(data) == (2, 100)
    
                # Evaluar CDF y PDF en cada columna de los datos generados
                for i in 1:d
                    u = data[:,i]
                    cdf_value = cdf(C, u)
                    pdf_value = pdf(C, u)
    
                    # Añadir mensaje de error detallado si el test falla
                    @test 0.0 <= cdf_value <= 1.0 || error("Fallo en la CDF: θ=$θ, u=$u, cdf_value=$cdf_value")
                    @test pdf_value >= 0.0 || error("Fallo en la PDF: θ=$θ, u=$u, pdf_value=$pdf_value")
                end
            catch e
                @test e isa ArgumentError
                println("No se pudo construir AsymMixedCopula con θ=$θ: ", e)
            end
        end
    end

    @testset "BC2Copula - sampling, pdf, cdf" begin
        for θ in [[rand(rng), rand(rng)], [1.0, 0.0], [0.5, 0.5], [-0.1, 0.2], [1.1, 0.5], [0.5, -0.2]]
            try
                C = BC2Copula(θ...)
                data = rand(rng, C, 100)
    
                # Verificar que la generación de números aleatorios tiene la forma correcta
                @test size(data) == (2, 100)
    
                # Evaluar CDF y PDF en cada columna de los datos generados
                for i in 1:d
                    u = data[:,i]
                    cdf_value = cdf(C, u)
                    pdf_value = pdf(C, u)
    
                    # Añadir mensaje de error detallado si el test falla
                    @test 0.0 <= cdf_value <= 1.0 || error("Fallo en la CDF: θ=$θ, u=$u, cdf_value=$cdf_value")
                    @test pdf_value >= 0.0 || error("Fallo en la PDF: θ=$θ, u=$u, pdf_value=$pdf_value")
                end
            catch e
                @test e isa ArgumentError
                println("No se pudo construir BC2Copula con θ=$θ: ", e)
            end
        end
    end

    @testset "CuadrasAugeCopula - sampling, pdf, cdf" begin
        for θ in [rand(rng), 0.0, 1.0, -0.1, 1.1, rand(rng)]
            try
                C = CuadrasAugeCopula(θ)
                data = rand(rng, C, 100)
    
                # Verificar que la generación de números aleatorios tiene la forma correcta
                @test size(data) == (2, 100)
    
                # Evaluar CDF y PDF en cada columna de los datos generados
                for i in 1:d
                    u = data[:,i]
                    cdf_value = cdf(C, u)
                    if θ == 1.0
                        # En el caso especial de θ = 1, se espera que PDF no esté definido
                        @test_throws ArgumentError pdf(C, u)
                    else
                        pdf_value = pdf(C, u)
                        # Añadir mensaje de error detallado si el test falla
                        @test 0.0 <= cdf_value <= 1.0 || error("Fallo en la CDF: θ=$θ, u=$u, cdf_value=$cdf_value")
                        @test pdf_value >= 0.0 || error("Fallo en la PDF: θ=$θ, u=$u, pdf_value=$pdf_value")
                    end
                end
            catch e
                @test e isa ArgumentError
                println("No se pudo construir CuadrasAugeCopula con θ=$θ: ", e)
            end
        end
    end

    @testset "GalambosCopula - sampling, pdf, cdf" begin
        for θ in [rand(rng), 0.0, Inf, -1.0, rand(rng, Uniform(1.0, 5.0)), rand(rng, Uniform(5.0, 10.0))]
            try
                C = GalambosCopula(θ)
                data = rand(rng, C, 100)
    
                # Verificar que la generación de números aleatorios tiene la forma correcta
                @test size(data) == (2, 100)
    
                # Evaluar CDF y PDF en cada columna de los datos generados
                for i in 1:d
                    u = data[:,i]
                    cdf_value = cdf(C, u)
                    if  θ == Inf
                        # En los casos especiales de θ = 0 o θ = Inf, se espera que PDF no esté definido
                        @test_throws ArgumentError pdf(C, u)
                    else
                        pdf_value = pdf(C, u)
                        # Añadir mensaje de error detallado si el test falla
                        @test 0.0 <= cdf_value <= 1.0 || error("Fallo en la CDF: θ=$θ, u=$u, cdf_value=$cdf_value")
                        @test pdf_value >= 0.0 || error("Fallo en la PDF: θ=$θ, u=$u, pdf_value=$pdf_value")
                    end
                end
            catch e
                @test e isa ArgumentError
                println("No se pudo construir GalambosCopula con θ=$θ: ", e)
            end
        end
    end

    @testset "HuslerReissCopula - sampling, pdf, cdf" begin
        for θ in [rand(rng), 0.0, Inf, -1.0, rand(rng, Uniform(1.0, 5.0)), rand(rng, Uniform(5.0, 10.0))]
            try
                C = HuslerReissCopula(θ)
                data = rand(rng, C, 100)
    
                # Verificar que la generación de números aleatorios tiene la forma correcta
                @test size(data) == (2, 100)
    
                # Evaluar CDF y PDF en cada columna de los datos generados
                for i in 1:d
                    u = data[:,i]
                    cdf_value = cdf(C, u)
                    if  θ == Inf
                        # En los casos especiales de θ = 0 o θ = Inf, se espera que PDF no esté definido
                        @test_throws ArgumentError pdf(C, u)
                    else
                        pdf_value = pdf(C, u)
                        # Añadir mensaje de error detallado si el test falla
                        @test 0.0 <= cdf_value <= 1.0 || error("Fallo en la CDF: θ=$θ, u=$u, cdf_value=$cdf_value")
                        @test pdf_value >= 0.0 || error("Fallo en la PDF: θ=$θ, u=$u, pdf_value=$pdf_value")
                    end
                end
            catch e
                @test e isa ArgumentError
                println("No se pudo construir HuslerReissCopula con θ=$θ: ", e)
            end
        end
    end

    @testset "LogCopula - sampling, pdf, cdf" begin
        for θ in [1.0, Inf, 0.5, rand(rng, Uniform(1.0, 10.0))]
            try
                C1 = LogCopula(θ)
                C2 = GumbelCopula(θ)
                data = rand(rng, C1, 100)
    
                # Verificar que la generación de números aleatorios tiene la forma correcta
                @test size(data) == (2, 100)
    
                # Evaluar CDF y PDF en cada columna de los datos generados
                for i in 1:d
                    u = data[:,i]
                    cdf_value_C1 = cdf(C1, u)
                    cdf_value_C2 = cdf(C2, u)
                    if  θ == Inf
                        # En los casos especiales de  θ = Inf, se espera que PDF no esté definido
                        @test_throws ArgumentError pdf(C1, u)
                        @test_throws ArgumentError pdf(C2, u)
                    else
                        pdf_value_C1 = pdf(C1, u)
                        pdf_value_C2 = pdf(C2, u)
                        # Añadir mensaje de error detallado si el test falla
                        @test 0.0 <= cdf_value_C1 <= 1.0 || error("Fallo en la CDF de LogCopula: θ=$θ, u=$u, cdf_value_C1=$cdf_value_C1")
                        @test 0.0 <= cdf_value_C2 <= 1.0 || error("Fallo en la CDF de GumbelCopula: θ=$θ, u=$u, cdf_value_C2=$cdf_value_C2")
                        @test pdf_value_C1 >= 0.0 || error("Fallo en la PDF de LogCopula: θ=$θ, u=$u, pdf_value_C1=$pdf_value_C1")
                        @test pdf_value_C2 >= 0.0 || error("Fallo en la PDF de GumbelCopula: θ=$θ, u=$u, pdf_value_C2=$pdf_value_C2")
    
                        # Comparar los valores de CDF y PDF entre LogCopula y GumbelCopula
                        @test isapprox(cdf_value_C1, cdf_value_C2, atol=1e-6) || error("CDF de LogCopula y GumbelCopula no coinciden: θ=$θ, u=$u, cdf_value_C1=$cdf_value_C1, cdf_value_C2=$cdf_value_C2")
                        @test isapprox(pdf_value_C1, pdf_value_C2, atol=1e-6) || error("PDF de LogCopula y GumbelCopula no coinciden: θ=$θ, u=$u, pdf_value_C1=$pdf_value_C1, pdf_value_C2=$pdf_value_C2")
                    end
                end
            catch e
                @test e isa ArgumentError
                println("No se pudo construir LogCopula con θ=$θ: ", e)
            end
        end
    end

    @testset "MixedCopula - sampling, pdf, cdf" begin
        for θ in [0.0, 1.0, 0.2, -rand(rng), 0.5]
            try
                C = MixedCopula(θ)
                data = rand(rng, C, 100)
    
                # Verificar que la generación de números aleatorios tiene la forma correcta
                @test size(data) == (2, 100)
    
                # Evaluar CDF y PDF en cada columna de los datos generados
                for i in 1:d
                    u = data[:,i]
                    cdf_value = cdf(C, u)
                    pdf_value = pdf(C, u)
                    # Añadir mensaje de error detallado si el test falla
                    @test 0.0 <= cdf_value <= 1.0 || error("Fallo en la CDF: θ=$θ, u=$u, cdf_value=$cdf_value")
                    @test pdf_value >= 0.0 || error("Fallo en la PDF: θ=$θ, u=$u, pdf_value=$pdf_value")
                end
            catch e
                @test e isa ArgumentError
                println("No se pudo construir HuslerReissCopula con θ=$θ: ", e)
            end
        end
    end

    @testset "MOCopula - sampling, pdf, cdf" begin
        # Lista de ejemplos de parámetros para probar
        param_sets = [
            (0.1, 0.2, 0.3),
            (1.0, 1.0, 1.0),
            (0.5, 0.5, 0.5),
            (-0.1, 0.2, 0.3),  # Caso inválido: parámetro negativo
            (0.1, -0.2, 0.3),  # Caso inválido: parámetro negativo
            (0.1, 0.2, -0.3),  # Caso inválido: parámetro negativo
            (rand(rng), rand(rng), rand(rng))  # Parámetros aleatorios
        ]
    
        for λ in param_sets
            try
                C = MOCopula(λ...)
                data = rand(rng, C, 100)
    
                # Verificar que la generación de números aleatorios tiene la forma correcta
                @test size(data) == (2, 100)
    
                # Evaluar CDF y PDF en cada columna de los datos generados
                for i in 1:d
                    u = data[:,i]
                    cdf_value = cdf(C, u)
                    pdf_value = pdf(C, u)
    
                    # Añadir mensaje de error detallado si el test falla
                    @test 0.0 <= cdf_value <= 1.0 || error("Fallo en la CDF: λ=$λ, u=$u, cdf_value=$cdf_value")
                    @test pdf_value >= 0.0 || error("Fallo en la PDF: λ=$λ, u=$u, pdf_value=$pdf_value")
                end
            catch e
                @test e isa ArgumentError
                println("No se pudo construir MOCopula con λ=$λ: ", e)
            end
        end
    end

    @testset "tEVCopula - sampling, pdf, cdf" begin
        # Lista de ejemplos de parámetros para probar
        param_sets = [
            (2.0, 0.5),  # ν > 0 y -1 < ρ <= 1
            (5.0, -0.5),  # ν > 0 y -1 < ρ <= 1
            (10.0, 1.0),  # ν > 0 y ρ == 1, debería crear MCopula
            (3.0, 0.0),  # ν > 0 y ρ == 0, debería crear IndependentCopula
            (-2.0, 0.5),  # Caso inválido: ν <= 0
            (3.0, 1.1),  # Caso inválido: ρ fuera de rango
            (3.0, -1.1),  # Caso inválido: ρ fuera de rango
            (rand(rng, Uniform(4.0, 10.0)), rand(rng, Uniform(-0.9, 1.0)))  # Parámetros aleatorios dentro de rango
        ]
    
        for (ν, ρ) in param_sets
            try
                C = tEVCopula(ν, ρ)
                data = rand(rng, C, 100)
    
                # Verificar que la generación de números aleatorios tiene la forma correcta
                @test size(data) == (2, 100)
    
                # Evaluar CDF y PDF en cada columna de los datos generados
                for i in 1:d
                    u = data[:,i]
                    
                    cdf_value = cdf(C, u)
                    if  ρ == 1
                        # En los casos especiales de θ = 0
                        @test_throws ArgumentError pdf(C, u)
                    else
                        pdf_value = pdf(C, u)
                        # Añadir mensaje de error detallado si el test falla
                        @test 0.0 <= cdf_value <= 1.0 || error("Fallo en la CDF: ν=$ν, ρ=$ρ, u=$u, cdf_value=$cdf_value")
                        @test pdf_value >= 0.0 || error("Fallo en la PDF: ν=$ν, ρ=$ρ, u=$u, pdf_value=$pdf_value")
                    end
                end
            catch e
                @test e isa ArgumentError
                println("No se pudo construir tEVCopula con ν=$ν, ρ=$ρ: ", e)
            end
        end
    end
end