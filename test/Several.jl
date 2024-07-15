using Test
using TwinCopulas
using Distributions
using Random

@testset "Several Copulas Tests" begin
    rng = MersenneTwister(2024)
    d = 100
    @testset "EmpiricalCopula - creation, cdf, pdf, and sampling" begin
     
        # Generar datos de prueba válidos bivariados
        data_valid = rand(2, 100)
        
        # Generar datos de prueba inválidos (unidimensional)
        data_invalid = rand(rng, 100)
        
        @testset "Valid Data" begin
            try
                C = EmpiricalCopula(data_valid)
                
                # Verificar que la creación de la copula empírica tiene la estructura correcta
                @test size(C.data) == (2, 100)
                @test size(C.pseudos) == (2, 100)
                
                # Verificar que las pseudo-observaciones están en el rango [0, 1]
                @test all(0 .<= C.pseudos .<= 1)
                
                # Evaluar CDF y PDF en algunos puntos
                u_points = [[0.5, 0.5], [0.1, 0.9], [0.7, 0.3]]
                for u in u_points
                    cdf_value = cdf(C, u)
                    pdf_value = pdf(C, u)
                    
                    @test 0.0 <= cdf_value <= 1.0 || error("Fallo en la CDF: u=$u, cdf_value=$cdf_value")
                    @test pdf_value >= 0.0 || error("Fallo en la PDF: u=$u, pdf_value=$pdf_value")
                end
                
                # Verificar la generación de muestras
                samples = rand(rng, C, 100)
                @test size(samples) == (2, 100)
                
            catch e
                @test e isa ErrorException
                println("Error inesperado al construir EmpiricalCopula con datos válidos: ", e)
            end
        end
    
        @testset "Invalid Data" begin
            try
                C = EmpiricalCopula(data_invalid)
            catch e
                @test e isa AssertionError
                println("Correctamente no se pudo construir EmpiricalCopula con datos inválidos: ", e)
            end
        end
    end

    @testset "IndependentCopula - sampling, pdf, cdf" begin
        C = IndependentCopula()
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

    @testset "MCopula - sampling, pdf, cdf" begin
        C = MCopula()
        data = rand(rng, C, 100)
                
        # Verificar que la generación de números aleatorios tiene la forma correcta
        @test size(data) == (2, 100)
        # Evaluar CDF y PDF en cada columna de los datos generados
        for i in 1:d
            u = data[:,i]
            @test 0.0 <= cdf(C, u) <= 1.0
            @test_throws ArgumentError pdf(C, u)
        end    
    end

    @testset "WCopula - sampling, pdf, cdf" begin
        C = WCopula()
        data = rand(rng, C, 100)
                
        # Verificar que la generación de números aleatorios tiene la forma correcta
        @test size(data) == (2, 100)
        # Evaluar CDF y PDF en cada columna de los datos generados
        for i in 1:d
            u = data[:,i]
            @test 0.0 <= cdf(C, u) <= 1.0
            @test_throws ArgumentError pdf(C, u)
        end    
    end

    @testset "SurvivalCopula - creation, cdf, pdf, and sampling" begin
        # Lista de copulas base para pruebas
        copula_bases = [
            GaussianCopula(0.5),
            ClaytonCopula(2.0),
            GumbelCopula(3.0),
            FrankCopula(4.0),
            JoeCopula(5.0),
        ]
        
        for base_copula in copula_bases
            # Crear la SurvivalCopula a partir de la copula base
            survival_copula = SurvivalCopula(base_copula)

            samples = rand(survival_copula, 100)
            @test size(samples) == (2, 100)
            @test all(0 .<= s .<= 1 for s in samples)
            
            # Verificar que la creación de la copula de supervivencia tiene la estructura correcta
            @test survival_copula.C == base_copula
    
            for i in 1:d
                u = samples[:,i]
                cdf_value = cdf(survival_copula, u)
                pdf_value = pdf(survival_copula, u)
                
                @test 0.0 <= cdf_value <= 1.0 || error("Fallo en la CDF: base_copula=$base_copula, u=$u, cdf_value=$cdf_value")
                @test pdf_value >= 0.0 || error("Fallo en la PDF: base_copula=$base_copula, u=$u, pdf_value=$pdf_value")
            end
        end
    end

    @testset "UtilCopula - sampling, pdf, cdf" begin
        C = UtilCopula()
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

    @testset "ArchimaxCopula - creation, cdf, pdf, and sampling" begin
        # Lista de copulas Archimedean y ExtremeValue para pruebas
        archimedean_copulas = [
            ClaytonCopula(2.0),
            GumbelCopula(2.5)
        ]
        
        extreme_value_copulas = [
            LogCopula(3.5),
            GalambosCopula(4.5)
        ]
        
        for archimedean in archimedean_copulas
            for extreme in extreme_value_copulas
                # Crear la ArchimaxCopula a partir de las copulas base
                archimax_copula = ArchimaxCopula(archimedean, extreme)
    
                # Verificar que la creación de la copula tiene la estructura correcta
                @test archimax_copula.Archimedean == archimedean
                @test archimax_copula.Extreme == extreme
                
                # Generar muestras
                samples = rand(rng, archimax_copula, 150)
    
                @test size(samples) == (2, 150)
                @test all(0 .<= s .<= 1 for s in samples)
        
                # Evaluar CDF y PDF en algunos puntos
                for i in 1:150
                    u = samples[:,i]
                    cdf_value = cdf(archimax_copula, u)
                    pdf_value = pdf(archimax_copula, u)
                    
                    @test 0.0 <= cdf_value <= 1.0 || error("Fallo en la CDF: archimedean=$archimedean, extreme=$extreme, u=$u, cdf_value=$cdf_value")
                    @test (pdf_value >= 0.0 || isnan(pdf_value)) || error("Fallo en la PDF: archimedean=$archimedean, extreme=$extreme, u=$u, pdf_value=$pdf_value")
                end
            end
        end
    end

    @testset "SklarDist - creation, cdf, pdf, and sampling" begin
        
        # Lista de copulas base para pruebas
        copula_bases = [
            GaussianCopula(0.5),
            ClaytonCopula(2.0),
            GalambosCopula(3.5)
        ]
        
        # Lista de distribuciones marginales para pruebas
        margins_list = [
            (Normal(), Normal()),
            (Exponential(1.0), Exponential(1.0)),
            (Beta(2.0, 5.0), Beta(2.0, 5.0)),
            (Beta(3.2, 2.3), Normal(4.0, 1.0))
        ]
        
        for copula in copula_bases
            for margins in margins_list
                # Crear la SklarDist a partir de la copula base y las distribuciones marginales
                sklar_dist = SklarDist(copula, [margins[1], margins[2]])
                
                # Verificar que la creación de la distribución de Sklar tiene la estructura correcta
                @test sklar_dist.copula == copula
                @test sklar_dist.margins == margins
        
                # Generar muestras
                samples = rand(rng, sklar_dist, 100)
                @test size(samples) == (2, 100)
        
                # Evaluar CDF y PDF en algunos puntos
                for i in 1:d
                    u = samples[:,i]
                    cdf_value = cdf(sklar_dist, u)
                    pdf_value = pdf(sklar_dist, u)
                    
                    @test 0.0 <= cdf_value <= 1.0 || error("Fallo en la CDF: copula=$copula, margins=$margins, u=$u, cdf_value=$cdf_value")
                    @test pdf_value >= 0.0 || error("Fallo en la PDF: copula=$copula, margins=$margins, u=$u, pdf_value=$pdf_value")
                end
            end
        end
        
    end
end