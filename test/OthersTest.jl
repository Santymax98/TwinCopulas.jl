using Test
using TwinCopulas
using Distributions
using Random

@testset "Others Copulas Tests" begin
    rng = MersenneTwister(2024)
    d = 100
    @testset "B11Copula - sampling, pdf, cdf" begin
        for θ in [0.0, rand(rng), rand(rng), 1.0]
            C = B11Copula(θ)
            data = rand(rng, C, 100)
                
            # Verificar que la generación de números aleatorios tiene la forma correcta
            @test size(data) == (2, 100)
            # Evaluar CDF y PDF en cada columna de los datos generados

            for i in 1:d
                u = data[:,i]
                if  θ == 1
                    @test_throws ArgumentError pdf(C, u)
                else
                @test 0.0 <= cdf(C, u) <= 1.0
                @test pdf(C, u) >= 0.0
                end  
            end
        end
    end

    @testset "FrechetCopula - sampling, pdf, cdf" begin
        # Lista de ejemplos de parámetros para probar
        param_sets = [
            (0.1, 0.2),  # θ1 + θ2 <= 1
            (0.5, 0.5),  # θ1 + θ2 <= 1
            (0.7, 0.3),  # θ1 + θ2 <= 1
            (-0.1, 0.2),  # Caso inválido: θ1 < 0
            (0.1, -0.2),  # Caso inválido: θ2 < 0
            (0.5, 0.6),  # Caso inválido: θ1 + θ2 > 1
            (1.1, 0.1),  # Caso inválido: θ1 > 1
            (0.1, 1.1),  # Caso inválido: θ2 > 1
            (rand(rng), rand(rng) * (1 - rand(rng)))  # Parámetros aleatorios dentro de rango
        ]
    
        for (θ1, θ2) in param_sets
            try
                C = FrechetCopula(θ1, θ2)
                data = rand(rng, C, 100)
    
                # Verificar que la generación de números aleatorios tiene la forma correcta
                @test size(data) == (2, 100)
    
                # Evaluar CDF y PDF en cada columna de los datos generados
                for i in 1:d
                    u = data[:,i]
                    cdf_value = cdf(C, u)
                    pdf_value = pdf(C, u)
    
                    # Añadir mensaje de error detallado si el test falla
                    @test 0.0 <= cdf_value <= 1.0 || error("Fallo en la CDF: θ1=$θ1, θ2=$θ2, u=$u, cdf_value=$cdf_value")
                    @test pdf_value >= 0.0 || error("Fallo en la PDF: θ1=$θ1, θ2=$θ2, u=$u, pdf_value=$pdf_value")
                end
            catch e
                @test e isa ArgumentError
                println("No se pudo construir FrechetCopula con θ1=$θ1, θ2=$θ2: ", e)
            end
        end
    end

    @testset "MardiaCopula - sampling, pdf, cdf" begin
        # Lista de ejemplos de parámetros para probar
        param_sets = [
            -1.0,  # θ = -1, debería crear WCopula
            -0.5,  # -1 < θ < 0
            0.0,   # θ = 0, debería crear IndependentCopula
            0.5,   # 0 < θ < 1
            1.0,   # θ = 1, debería crear MCopula
            -1.1,  # Caso inválido: θ < -1
            1.1,   # Caso inválido: θ > 1
            rand(rng, Uniform(-1.0, 1.0))  # Parámetro aleatorio dentro de rango
        ]
    
        for θ in param_sets
            try
                C = MardiaCopula(θ)
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
                println("No se pudo construir MardiaCopula con θ=$θ: ", e)
            end
        end
    end

    @testset "MaresiasCopula - sampling, cdf" begin
    
        # Definición de funciones G válidas
        function G1(u)
            return u > 0.5 ? 2u - 1 : 0
        end
        
        function G2(u)
            α = rand(rng, Uniform(1.0, 2.0))
            return u^α
        end
        
        function G3(u)
            α = 5.5 # α ≥ 1
            return  u < 0.5 ? u / α : (2 - 1/α) * u - (1 - 1/α) * u
        end
        

        valid_functions = [G1, G2, G3]
    
        for G in valid_functions
            try
                C = MaresiasCopula(G)
                data = rand(rng, C, 100)
    
                # Verificar que la generación de números aleatorios tiene la forma correcta
                @test size(data) == (2, 100)
    
                # Evaluar CDF en cada columna de los datos generados
                for i in 1:d
                    u = data[:,i]
                    cdf_value = cdf(C, u)
    
                    # Añadir mensaje de error detallado si el test falla
                    @test 0.0 <= cdf_value <= 1.0 || error("Fallo en la CDF: G=$G, u=$u, cdf_value=$cdf_value")
                    @test pdf(C, u) >= 0.0
                end
            catch e
                @test e isa ArgumentError
                println("No se pudo construir MaresiasCopula con G=$G: ", e)
            end
        end
    end

    @testset "MorgensternCopula - sampling, pdf, cdf" begin
        for θ in [-1.0, -rand(rng), 0.0, rand(rng), 1.0]
            C = MorgensternCopula(θ)
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

    @testset "PlackettCopula - sampling, pdf, cdf" begin

        for θ in [0.0, 1.0, rand(rng), rand(rng, Uniform(1.5, 5.5)), rand(rng, Uniform(6.0, 10.0)), Inf]

            C = PlackettCopula(θ)
            data = rand(rng, C, 100)
    
            # Verificar que la generación de números aleatorios tiene la forma correcta
            @test size(data) == (2, 100)
    
            # Evaluar CDF y PDF en cada columna de los datos generados
            for i in 1:d
                u = data[:,i]
                cdf_value = cdf(C, u)
                if  θ == Inf || θ == 0
                    # En los casos especiales de θ = 0
                    @test_throws ArgumentError pdf(C, u)
                else
                    pdf_value = pdf(C, u)
                    # Añadir mensaje de error detallado si el test falla
                    @test 0.0 <= cdf_value <= 1.0
                    @test pdf_value >= 0.0
                end
            end
        end
    end

    @testset "RafteryCopula - sampling, pdf, cdf" begin

        for θ in [0.0, 1.0, rand(rng), rand(rng)]

            C = RafteryCopula(θ)
            data = rand(rng, C, 100)
    
            # Verificar que la generación de números aleatorios tiene la forma correcta
            @test size(data) == (2, 100)
    
            # Evaluar CDF y PDF en cada columna de los datos generados
            for i in 1:d
                u = data[:,i]
                cdf_value = cdf(C, u)
                if  θ == 1
                    # En los casos especiales de θ = 0
                    @test_throws ArgumentError pdf(C, u)
                else
                    pdf_value = pdf(C, u)
                    # Añadir mensaje de error detallado si el test falla
                    @test 0.0 <= cdf_value <= 1.0
                    @test pdf_value >= 0.0
                end
            end
        end
    end
end
