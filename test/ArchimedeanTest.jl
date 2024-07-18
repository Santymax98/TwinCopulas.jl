using Test
using TwinCopulas
using Distributions
using Random

@testset "Archimedean Copulas Tests" begin
    rng = MersenneTwister(2024)
    d = 100
    @testset "AMHCopula - sampling, pdf, cdf" begin
        for θ in [-1.0, -rand(rng), 0.0, rand(rng)]
            C = AMHCopula(θ)
            data = rand(rng, C, 100)
                
            # Verificar que la generación de números aleatorios tiene la forma correcta
            @test size(data) == (2, 100)
            # Evaluar CDF y PDF en cada columna de los datos generados

            for i in 1:d
                u = data[:,i]
                @test 0.0 <= cdf(C, u) <= 1.0
                @test pdf(C, u) >= 0.0
            end    
        end
    end

    @testset "ClaytonCopula - sampling, pdf, cdf" begin
        for θ in [-1.0, -rand(rng), 0.0, rand(rng, Uniform(10, 20)), Inf]
            C = ClaytonCopula(θ)
            data = rand(rng, C, 100)
                
            # Verificar que la generación de números aleatorios tiene la forma correcta
            @test size(data) == (2, 100)
            # Evaluar CDF y PDF en cada columna de los datos generados

            for i in 1:d
                u = data[:,i]
                if isinf(θ) || (θ == -1)
                    # En el caso especial de θ = Inf, se espera que PDF no esté definido
                    @test_throws ArgumentError pdf(C, u)
                else
                @test 0.0 <= cdf(C, u) <= 1.0
                @test pdf(C, u) >= 0.0
                end
            end    
        end
    end

    @testset "FrankCopula - sampling, pdf, cdf" begin
        for θ in [-Inf, -rand(rng, Uniform(10, 20)), 0.0, rand(rng, Uniform(10, 20)), Inf]
            C = FrankCopula(θ)
            data = rand(rng, C, 100)
                
            # Verificar que la generación de números aleatorios tiene la forma correcta
            @test size(data) == (2, 100)
            # Evaluar CDF y PDF en cada columna de los datos generados

            for i in 1:d
                u = data[:,i]
                if isinf(θ)
                    # En el caso especial de θ = Inf, se espera que PDF no esté definido
                    @test_throws ArgumentError pdf(C, u)
                else
                @test 0.0 <= cdf(C, u) <= 1.0
                @test pdf(C, u) >= 0.0
                end
            end    
        end
    end

    @testset "GumbelBarnettCopula - sampling, pdf, cdf" begin
        for θ in [0.0, rand(rng), rand(rng), 1.0]
            C = GumbelBarnettCopula(θ)
            data = rand(rng, C, 100)
                
            # Verificar que la generación de números aleatorios tiene la forma correcta
            @test size(data) == (2, 100)
            # Evaluar CDF y PDF en cada columna de los datos generados

            for i in 1:d
                u = data[:,i]
                @test 0.0 <= cdf(C, u) <= 1.0
                @test pdf(C, u) >= 0.0
            end    
        end
    end

    @testset "GumbelCopula- sampling, pdf, cdf" begin
        for θ in [1.0, rand(rng, Uniform(1, 10)), rand(rng, Uniform(10, 20)), Inf]
            C = GumbelCopula(θ)
            data = rand(rng, C, 100)
                
            # Verificar que la generación de números aleatorios tiene la forma correcta
            @test size(data) == (2, 100)
            # Evaluar CDF y PDF en cada columna de los datos generados

            for i in 1:d
                u = data[:,i]
                if isinf(θ)
                    @test_throws ArgumentError pdf(C, u)
                else
                @test 0.0 <= cdf(C, u) <= 1.0
                @test pdf(C, u) >= 0.0
                end
            end    
        end
    end

    @testset "InvGaussianCopula - sampling, pdf, cdf" begin
        for θ in [0.01, rand(rng, Uniform(5, 10)), rand(rng, Uniform(10, 15))]
            C = InvGaussianCopula(θ)
            data = rand(rng, C, 100)
                
            # Verificar que la generación de números aleatorios tiene la forma correcta
            @test size(data) == (2, 100)
            # Evaluar CDF y PDF en cada columna de los datos generados

            for i in 1:d
                u = data[:,i]
                @test 0.0 <= cdf(C, u) <= 1.0
                @test pdf(C, u) >= 0.0
            end    
        end
    end

    @testset "JoeCopula - sampling, pdf, cdf" begin
        for θ in [1.0, rand(rng, Uniform(1, 10)), rand(rng, Uniform(5, 15)), Inf]
            C = JoeCopula(θ)
            data = rand(rng, C, 100)
                
            # Verificar que la generación de números aleatorios tiene la forma correcta
            @test size(data) == (2, 100)
            # Evaluar CDF y PDF en cada columna de los datos generados

            for i in 1:d
                u = data[:,i]
                if isinf(θ)
                    # En el caso especial de θ = Inf, se espera que PDF no esté definido
                    @test_throws ArgumentError pdf(C, u)
                else
                @test 0.0 <= cdf(C, u) <= 1.0
                @test pdf(C, u) >= 0.0
                end
            end    
        end
    end

    @testset "Nelsen2Copula - sampling, pdf, cdf" begin
        for θ in [1.0, rand(rng, Uniform(5.0, 10.0)), rand(rng, Uniform(15.0, 20.0)), Inf]
            C = Nelsen2Copula(θ)
            data = rand(rng, C, 100)
                
            # Verificar que la generación de números aleatorios tiene la forma correcta
            @test size(data) == (2, 100)
            # Evaluar CDF y PDF en cada columna de los datos generados

            for i in 1:d
                u = data[:,i]
                if isinf(θ) || θ == 1
                    # En el caso especial de θ = Inf, se espera que PDF no esté definido
                    @test_throws ArgumentError pdf(C, u)
                else
                    cdf_value = cdf(C, u)
                    if !(0.0 <= cdf_value <= 1.0)
                        error("Fallo en la CDF: θ=$θ, u=$u, cdf_value=$cdf_value")
                    end
                    @test 0.0 <= cdf_value <= 1.0
                    @test pdf(C, u) >= 0.0
                end
            end    
        end
    end
    
end
