using Test
using TwinCopulas
using Distributions
using Random

@testset "Elliptical Copulas Tests" begin
    rng = MersenneTwister(2024)
    d = 100
    @testset "GaussianCopula - sampling, pdf, cdf" begin
        for θ in [-1.0, -rand(rng), 0.0, rand(rng), 1.0]
            C = GaussianCopula(θ)
            data = rand(rng, C, 100)
                
            # Verificar que la generación de números aleatorios tiene la forma correcta
            @test size(data) == (2, 100)
            # Evaluar CDF y PDF en cada columna de los datos generados

            for i in 1:d
                u = data[:,i]
                if θ == -1 || θ == 1
                    @test_throws ArgumentError pdf(C, u)
                else
                @test 0.0 <= cdf(C, u) <= 1.0
                @test pdf(C, u) >= 0.0
                end  
            end
        end
    end

    @testset "tCopula - sampling, pdf, cdf" begin
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
                C = tCopula(ν, ρ)
                data = rand(rng, C, 100)
    
                # Verificar que la generación de números aleatorios tiene la forma correcta
                @test size(data) == (2, 100)
    
                # Evaluar CDF y PDF en cada columna de los datos generados
                for i in 1:d
                    u = data[:,i]
                    cdf_value = cdf(C, u)
                    if  ρ == 1
                        # En los casos especiales de θ = 1
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
                println("No se pudo construir copula t con θ=$ν, ρ=$ρ: ", e)
            end
        end
    end
end
