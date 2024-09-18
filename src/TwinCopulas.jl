module TwinCopulas
    
    import Base
    import Random
    import SpecialFunctions
    import Roots
    import Distributions
    import StatsBase
    import StatsFuns
    import ForwardDiff
    import Cubature
    import WilliamsonTransforms
    import LinearAlgebra
    import QuadGK
    import PolyLog
    import LogExpFunctions
    import Statistics
    import Interpolations
    
    #principal module
    include("bicopula.jl")
    include("Archimedean.jl")
    include("ExtremeValue.jl")
    include("Elliptical.jl")
    include("SklarDist.jl")
    
    export 
    ArchimedeanCopula,
    ExtremeValueCopula,
    EllipticalCopula,
    SklarDist

    #Utils
    include("Distributions/ExtremeDist.jl")
    include("Distributions/Logarithmic.jl")
    include("Distributions/PStable.jl")
    include("Distributions/RadialDist.jl") #To test API Radial Distributions
    include("Distributions/Sibuya.jl")
    include("Distributions/Stable.jl")

    #Elliptical Copulas
    include("EllipticalCopulas/GaussianCopula.jl")
    include("EllipticalCopulas/tCopula.jl")

    export 

    GaussianCopula,
    tCopula
    

    #Archimedean Copulas
    include("ArchimedeanCopulas/AMHCopula.jl")
    include("ArchimedeanCopulas/ClaytonCopula.jl")
    include("ArchimedeanCopulas/FrankCopula.jl")
    include("ArchimedeanCopulas/GumbelCopula.jl")
    include("ArchimedeanCopulas/GumbelBarnettCopula.jl")
    include("ArchimedeanCopulas/InvGaussianCopula.jl")
    include("ArchimedeanCopulas/JoeCopula.jl")
    include("ArchimedeanCopulas/Nelsen2Copula.jl") #To test API Radial Distributions
    include("ArchimedeanCopulas/UtilCopula.jl")
        
    export
        AMHCopula,
        ClaytonCopula,
        GumbelCopula,
        GumbelBarnettCopula,
        FrankCopula,
        InvGaussianCopula,
        JoeCopula,
        Nelsen2Copula, #To test API Radial Distributions
        UtilCopula #Special case of some copulas

    # Extreme value Copulas
    
    include("ExtremeValueCopulas/AsymGalambosCopula.jl")
    include("ExtremeValueCopulas/AsymLogCopula.jl")
    include("ExtremeValueCopulas/AsymMixedCopula.jl")
    include("ExtremeValueCopulas/BC2Copula.jl")
    include("ExtremeValueCopulas/CuadrasAugeCopula.jl")
    include("ExtremeValueCopulas/GalambosCopula.jl")
    include("ExtremeValueCopulas/HuslerReissCopula.jl")
    include("ExtremeValueCopulas/LogCopula.jl")
    include("ExtremeValueCopulas/MixedCopula.jl")
    include("ExtremeValueCopulas/MOCopula.jl")    
    include("ExtremeValueCopulas/tEVCopula.jl")
    
    export
        AsymGalambosCopula,
        AsymLogCopula,
        AsymMixedCopula,
        BC2Copula,
        CuadrasAugeCopula,
        GalambosCopula,
        HuslerReissCopula,
        LogCopula,
        MixedCopula,
        MOCopula,   
        tEVCopula
    
    #important Copulas
    include("SeveralCopulas/ArchimaxCopula.jl")
    include("SeveralCopulas/EmpiricalCopula.jl")
    include("SeveralCopulas/IndependentCopula.jl")
    include("SeveralCopulas/MCopula.jl")
    include("SeveralCopulas/SurvivalCopula.jl")
    include("SeveralCopulas/WCopula.jl")
    
    export
        ArchimaxCopula,
        EmpiricalCopula,
        IndependentCopula,
        MCopula,
        SurvivalCopula,
        WCopula

    #Other practical Copulas
    include("Others/B11Copula.jl")
    include("Others/FrechetCopula.jl")
    include("Others/MardiaCopula.jl")
    include("Others/MaresiasCopula.jl")
    include("Others/MorgensternCopula.jl")
    include("Others/PlackettCopula.jl")
    include("Others/RafteryCopula.jl")
    
    export 
    B11Copula,
    FrechetCopula,
    MardiaCopula,
    MaresiasCopula,
    MorgensternCopula,
    PlackettCopula,
    RafteryCopula


end
