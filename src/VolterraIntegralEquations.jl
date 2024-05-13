#=
Created on Tue 27 Jun 2023
Updated on Mon 13 May 2024
=#

#=
References:
Braß, H., "On the Principle of Avoiding the Singularity in Quadrature," Zeitschrift für angewandte Mathematik und Mechanik, 75 (1995), S617-S618.

Caputo, M. and Mainardi, F., "Linear models of dissipation in anelastic solids," Rivista del Nuoro Cimento, 1 (1971), 161-198.

Caputo, M. and Mainardi, F., "A new dissipation model based on memory mechanism," Pure and Applied Geophysics, 91 (1971), 134-147.

Cole, K.S. and Cole, R.H., "Dispersion and absorption in dielectrics I. Alternating current characteristics," Journal of Chemical Physics, 9 (1941), 342-351.

Cole, K.S. and Cole, R.H., "Dispersion and absorption in dielectrics II. Direct current characteristics," Journal of Chemical Physics, 10 (1942), 98-105.

Freed, A.D., Soft Solids: A primer to the theoretical mechanics of materials, Modeling and Simulation in Science, Engineering and Technology. Basel: Birkhäuser, 2014.

Freed, A.D. and Rajagopal, K.R. "A viscoelastic model for describing the response of biological fibers," ACTA Mechanica, 227 (2016), 3367-3380.

Fung, Y.-C., "Biorheology of Soft Tissues," Biorheology, 10 (1973), 139-155.

Kohlrausch, R., "Ueber das Dellmann'sche Elektrometer," Annalen der Physik und Chemie, 72 (1847), 353-405.

Maxwell, J.C., "On the dynamical theory of gases," Philosophical Transactions of the Royal Society, London, 157 (1867), 49-88.

Neubert, H.K., "A simple model representing internal damping in solid materials," The Aeronautical Quarterly, 14 (1963), 187-210.

Volterra, V. Theory of functionals and of integral and integro-differential equations. Glasgow: Blackie and Son, 1930.

Williams, G. and Watts, D.C., "Non-symmetrical dielectric relaxation behaviour arising from a simple empirical decay function," Transactions of the Faraday Society, 66 (1970), 80-85.

Williams, M.L., "Structural analysis of viscoelastic materials," AIAA Journal, 2 (1964), 785-808.

Young, A., "Approximate product integration," Proceedings of the Royal Society, London, A-224 (1954), 552-561.

Zener, C., Elasticity and Anelasticity of Metals. Chicago: University of Chicago Press, 1948.
=#

module VolterraIntegralEquations

using
    JSON3,
    MittagLeffler,
    PhysicalFields,
    StructTypes

export
    # Memory functions, which are kernels to Volterra integral equations.
    BOX,    # the 𝑏𝑜𝑥 model of Neuber, a.k.a. Fung's QLV kernel
    CCM,    # 𝐶ole and 𝐶ole's power-law 𝑀odel
    FLS,    # Caputo and Mainardi's 𝐹ractional 𝐿inear 𝑆olid
    KWW,    # 𝐾ohlrausch's and 𝑊illiams & 𝑊atts' stretched exponential
    MCM,    # 𝑀axwell's 𝐶hain 𝑀odel, a.k.a. the Prony series model
    MPL,    # Williams' 𝑀odified 𝑃ower-𝐿aw model
    RFS,    # Freed and Rajagopal's 𝑅egularized 𝐹L𝑆 model
    SLS,    # Zener's 𝑆tandard 𝐿inear 𝑆olid, a.k.a. the kernel of Maxwell-Debye

    # Function used to create weights of quadrature for a given memory function.
    normalizedQuadratureWeights,

    # Solvers for Volterra integral equations of the second kind; specifically,
    #   f′(t) = g′(t) - c(t) ∫₀ᵗ K(t-τ) f′(τ) dτ
    # where
    #   g′ is the time rate-of-change of some control function g(t)
    #   f′ is the time rate-of-change of the response function f(t)
    #   c  is a scalar function, e.g., (E₀ - E∞)/E∞ in linear viscoelasticity
    #   K  is a memory function, i.e., the derivative of a creep function
    # Here f′ and g′ may be scalar, vector or tensor valued.
    # Upon solving f′, the resulting ODE can be integrated to get response f.

    # abstract type
    VolterraIntegralEquation,

    # concrete types with internal constructors
    VolterraIntegralScalarEquation,
    VolterraIntegralVectorEquation,
    VolterraIntegralTensorEquation,

    # methods

    copy,
    deepcopy,
    toFile,
    fromFile,

    advance!,
    update!

#=
-------------------------------------------------------------------------------
=#

#= 
Memory functions are the derivatives of creep functions [Freed, 2014], the latter being more commonly found in the literature. All memory functions are to have an interface of:
    (name, k, τ) = <memoryFunctionName>(systemOfUnits, time, parameters)
which returns a tuple whose first entry is a string specifying the name of this kernel, e.g., "FLS", its second entry contains a value for this memory function `k` evaluated at `time,` and its third entry contains the kernel's controlling characteristic time `τ,` which is the smallest one whenever multiple characteristic times exist. At present, the supplied argument `systemOfUnits` to a kernel call can be either "SI" or "CGS". Its second argument `time` contains the current time, while its final argument `parameters` is a tuple containing this kernel's physical parameters.

Specifically, the following memory functions have been implemented:
    BOX     the box kernel of Neuber, a.k.a. Fung's QLV kernel
    CCM     Cole and Cole's power-law Model
    FLS     Caputo and Mainardi's Fractional Linear Solid
    KWW     Kohlrausch's and Williams & Watts' stretched exponential
    MCM     Maxwell's Chain Model, a.k.a. the Prony series kernel
    MPL     Williams' Modified Power-Law
    RFS     Freed and Rajagopal's Regularized FLS
    SLS     Zener's Standard Linear Solid, a.k.a. the Maxwell-Debye kernel
whose material parameters are supplied via the following tuples:
    BOX     parameters = (τ₁, τ₂)
    CCM     parameters = (α, τ)
    FLS     parameters = (α, τ)
    KWW     parameters = (α, τ)
    MCM     parameters = (c₁, c₂, …, cₙ, τ₁, τ₂ …, τₙ)
    MPL     parameters = (α, τ)
    RFS     parameters = (α, δ, τ)
    SLS     parameters = (τ,)
wherein τ denotes a characteristic time for creep. There are two characteristic times in the BOX model, and n in the MCM, arranged so that 0 < τ₁ < τ₂ < ⋯ < τₙ, with each cᵢ > 0, i = 1, 2, …, n, being a coefficient in the Prony series whose collective sum is 1, i.e., ∑_{i=1}^n cᵢ = 1. Parameter α is the exponent in a power law, while parameter δ shifts time to remove a weak singularity.

The following memory functions are weakly singular at the upper limit of integration in their Volterra integrals:
    CCM, FLS and KWW.
Consequently, the existence of such kernels requires a numerical method that avoids this possible singularity which can occur at the upper limit of integration. Gauss' quadrature for integrating an integral has this property.
=#

"""
Memory function BOX (the box energy function of Neubert)\n
    k = (exp(-t/τ₂) - exp(-t/τ₁)) / (t ln(τ₂/τ₁)) \n
    τ = τ₁\n
Argument `parameters` describes the tuple (τ₁, τ₂) ordered as 0 < τ₁ < τ₂.
"""
function BOX(systemOfUnits::String, time::PhysicalScalar, parameters::Tuple)::Tuple
    # Ensure that a consistent system of physical units is used.
    if (systemOfUnits == "SI") || (systemOfUnits == "si")
        t = toSI(time)
        τ₁ = toSI(parameters[1])
        τ₂ = toSI(parameters[2])
        if (t.units ≠ SECOND) || (τ₁.units ≠ SECOND) || (τ₂.units ≠ SECOND)
            msg = string("Argument time and parameters τ₁ and τ₂ must have units of time.")
            throw(ErrorException(msg))
        end
    elseif (systemOfUnits == "CGS") || (systemOfUnits == "cgs")
        t = toCGS(time)
        τ₁ = toCGS(parameters[1])
        τ₂ = toCGS(parameters[2])
        if ((t.units ≠ CGS_SECOND) || (τ₁.units ≠ CGS_SECOND) ||
            (τ₂.units ≠ CGS_SECOND))
            msg = string("Argument time and parameter τ₁ and τ₂ must have units of time.")
            throw(ErrorException(msg))
        end
    else
        msg = "The assigned physical system of units is unknown."
         throw(ErrorException(msg))
    end

    if (τ₁.value > 0.0) && (τ₂ > τ₁)
        if t.value ≈ 0.0
            k = (1/τ₁ - 1/τ₂) / log(τ₂/τ₁)
        elseif t.value > 0.0
            k = (exp(-t/τ₂) - exp(-t/τ₁)) / (t * log(τ₂/τ₁))
        else
            msg = "Argument time must be non-negative."
            throw(ErrorException(msg))
        end
    else
        msg = "Parameters τ₁ and τ₂ must be ordered as 0 < τ₁ < τ₂."
        throw(ErrorException(msg))
    end
    return ("BOX", k, τ₁)
end # BOX

"""
Memory function CCM (Cole-Cole power-law Model)\n
    k = (t/τ)^α (α / t) / (1 + (t/τ)^α)² \n
    τ = τ\n
Argument `parameters` describes the tuple (α, τ).
"""
function CCM(systemOfUnits::String, time::PhysicalScalar, parameters::Tuple)::Tuple
    # Ensure that a consistent system of physical units is used.
    if (systemOfUnits == "SI") || (systemOfUnits == "si")
        t = toSI(time)
        α = toSI(parameters[1])
        τ = toSI(parameters[2])
        if ((t.units ≠ SECOND) || (α.units ≠ DIMENSIONLESS) ||
            (τ.units ≠ SECOND))
            msg = string("Argument time and parameter τ must have units of time, while parameter α must be dimensionless.")
            throw(ErrorException(msg))
        end
    elseif (systemOfUnits == "CGS") || (systemOfUnits == "cgs")
        t = toCGS(time)
        α = toCGS(parameters[1])
        τ = toCGS(parameters[2])
        if ((t.units ≠ CGS_SECOND) || (α.units ≠ CGS_DIMENSIONLESS) ||
            (τ.units ≠ CGS_SECOND))
            msg = string("Argument time and parameter τ must have units of time, while parameter α must be dimensionless.")
            throw(ErrorException(msg))
        end
    else
        msg = "The assigned physical system of units is unknown."
         throw(ErrorException(msg))
    end

    if (α.value > 0.0) && (τ.value > 0.0)
        if t.value == 0.0
            k = PhysicalScalar(-t.units)
            set!(k, Inf)
        elseif t.value > 0.0
            x = (t/τ)^get(α)
            k = x * (α/t) / ((1 + x) * (1 + x))
        else
            msg = "Argument time must be non-negative."
            throw(ErrorException(msg))
        end
    else
        msg = "Parameters α and τ must be positive."
        throw(ErrorException(msg))
    end
    return ("CCM", k, τ)
end # CCM

"""
Memory function FLS (Fractional Linear Solid)\n
    k = -E_{α,0}(-(t/τ)^α) / t \n
    τ = τ\n
where E_{α, β}(z) denotes the two-parameter Mittag-Leffler function, with β = 0 in this case. Argument `parameters` describes the tuple (α, τ).
"""
function FLS(systemOfUnits::String, time::PhysicalScalar, parameters::Tuple)::Tuple
    # Ensure that a consistent system of physical units is used.
    if (systemOfUnits == "SI") || (systemOfUnits == "si")
        t = toSI(time)
        α = toSI(parameters[1])
        τ = toSI(parameters[2])
        if (t.units ≠ SECOND) || (α.units ≠ DIMENSIONLESS) || (τ.units ≠ SECOND)
            msg = string("Argument time and parameter τ must have units of time, while parameter α must be dimensionless.")
            throw(ErrorException(msg))
        end
    elseif (systemOfUnits == "CGS") || (systemOfUnits == "cgs")
        t = toCGS(time)
        α = toCGS(parameters[1])
        τ = toCGS(parameters[2])
        if ((t.units ≠ CGS_SECOND) || (α.units ≠ CGS_DIMENSIONLESS) ||
            (τ.units ≠ CGS_SECOND))
            msg = string("Argument time and parameter τ must have units of time, while parameter α must be dimensionless.")
            throw(ErrorException(msg))
        end
    else
        msg = "The assigned physical system of units is unknown."
         throw(ErrorException(msg))
    end

    if α.value ≈ 1.0
        (k, τ) = SLS(systemOfUnits, time, (τ,))
    elseif (τ.value > 0.0) && (α.value > 0.0) && (α.value < 1.0)
        if t.value == 0.0
            k = PhysicalScalar(-t.units)
            set!(k, Inf)
        elseif t.value > 0.0
            x = (t / τ)^get(α)
            k = -MittagLeffler.mittleff(get(α), 0.0, -get(x)) / t
        else
            msg = "Argument time must be non-negative."
            throw(ErrorException(msg))
        end
    else
        msg = "Parameter τ must be positive, and parameter α ∈ (0,1]."
        throw(ErrorException(msg))
    end
    return ("FLS", k, τ)
end # FLS

"""
Memory function KWW (stretched exponential of Kohlrausch, Williams and Watts)\n
    k = (t/τ)^α (α/t) exp(-(t/τ)^α) \n
    τ = τ\n
Argument `parameters` describes the tuple (α, τ).
"""
function KWW(systemOfUnits::String, time::PhysicalScalar, parameters::Tuple)::Tuple
    # Ensure that a consistent system of physical units is used.
    if (systemOfUnits == "SI") || (systemOfUnits == "si")
        t = toSI(time)
        α = toSI(parameters[1])
        τ = toSI(parameters[2])
        if (t.units ≠ SECOND) || (α.units ≠ DIMENSIONLESS) || (τ.units ≠ SECOND)
            msg = string("Argument time and parameter τ must have units of time, while parameter α must be dimensionless.")
            throw(ErrorException(msg))
        end
    elseif (systemOfUnits == "CGS") || (systemOfUnits == "cgs")
        t = toCGS(time)
        α = toCGS(parameters[1])
        τ = toCGS(parameters[2])
        if ((t.units ≠ CGS_SECOND) || (α.units ≠ CGS_DIMENSIONLESS) ||
            (τ.units ≠ CGS_SECOND))
            msg = string("Argument time and parameter τ must have units of time, while parameter α must be dimensionless.")
            throw(ErrorException(msg))
        end
    else
        msg = "The assigned physical system of units is unknown."
         throw(ErrorException(msg))
    end

    if α.value ≈ 1.0
        (k, τ) = SLS(systemOfUnits, time, (τ,))
    elseif (τ.value > 0.0) && (α.value > 0.0) && (α.value < 1.0)
        if t.value == 0.0
            k = PhysicalScalar(-t.units)
            set!(k, Inf)
        elseif t.value > 0.0
            k = (t/τ)^α.value * (α/t) * exp(-(t/τ)^α.value)
        else
            msg = "Argument time must be non-negative"
            throw(ErrorException(msg))
        end
    else
        msg = "Parameter τ must be positive and parameter α ∈ (0,1]."
        throw(ErrorException(msg))
    end
    return ("KWW", k, τ)
end # KWW

"""
Memory function MCM (Maxwell Chain Model, which is a Prony series)\n
    k = (c₁/τ₁) exp(-t/τ₁) + ⋯ + (cₙ/τₙ) exp(-t/τₙ) \n
    τ = τ₁\n
Argument `parameters` describes a tuple (c₁, c₂, …, cₙ, τ₁, τ₂, …, τₙ) of length 2n, where ∑_{i=1}^n cᵢ = 1, and where 0 < τ₁ < τ₂ < ⋯ < τₙ.
"""
function MCM(systemOfUnits::String, time::PhysicalScalar, parameters::Tuple)::Tuple
    # Verify the inputs.
    if length(parameters) % 2 == 0
        n = length(parameters) ÷ 2
    else
        msg = "There must be an even number of parameters, viz. paired sets."
        throw(ErrorException(msg))
    end

    # Ensure that a consistent system of physical units is used.
    if (systemOfUnits == "SI") || (systemOfUnits == "si")
        t = toSI(time)
        c = ArrayOfPhysicalScalars(n, DIMENSIONLESS)
        τ = ArrayOfPhysicalScalars(n, SECOND)
        for i in 1:n
            c[i] = toSI(parameters[i])
            τ[i] = toSI(parameters[n+i])
        end
        sum = PhysicalScalar(DIMENSIONLESS)
        for i in 1:n
            sum = sum + c[i]
        end
        if !(sum ≈ 1.0)
            msg = "Prony series coefficients must sum to 1."
            throw(ErrorException(msg))
        end
        for i in 2:n
            if τ[i-1] ≥ τ[i]
                msg = "Prony characteristic times must order 0 < τ₁ < ⋯ < τₙ."
                throw(ErrorException(msg))
            end
        end
    elseif (systemOfUnits == "CGS") || (systemOfUnits == "cgs")
        t = toCGS(time)
        c = ArrayOfPhysicalScalars(n, CGS_DIMENSIONLESS)
        τ = ArrayOfPhysicalScalars(n, CGS_SECOND)
        for i in 1:n
            c[i] = toCGS(parameters[i])
            τ[i] = toCGS(parameters[n+i])
        end
        sum = PhysicalScalar(CGS_DIMENSIONLESS)
        for i in 1:n
            sum = sum + c[i]
        end
        if !(sum ≈ 1.0)
            msg = "Prony series coefficients must sum to 1."
            throw(ErrorException(msg))
        end
        for i in 2:n
            if τ[i-1] ≥ τ[i]
                msg = "Prony characteristic times must order 0 < τ₁ < ⋯ < τₙ."
                throw(ErrorException(msg))
            end
        end
    else
        msg = "The assigned physical system of units is unknown."
         throw(ErrorException(msg))
    end

    if (t.value ≥ 0.0) && (τ[1].value > 0.0)
        k = PhysicalScalar(-t.units)
        if t.value ≈ 0.0
            for i in 1:n
                k = k + c[i] / τ[i]
            end
        else
            for i in 1:n
                k = k + (c[i] / τ[i]) * exp(-t/τ[i])
            end
        end
        τ₁ = τ[1]
    else
        msg = "Argument time must be non-negative, and parameters τᵢ must be positive."
        throw(ErrorException(msg))
    end
    return ("MCM", k, τ₁)
end # MCM

"""
Memory function MPL (Modified Power-Law model of Williams')\n
    k = (α/τ) / (1 + t/τ)^(1+α) \n
    τ = τ\n
Argument `parameters` describes the tuple (α, τ).
"""
function MPL(systemOfUnits::String, time::PhysicalScalar, parameters::Tuple)::Tuple
    # Ensure that a consistent system of physical units is used.
    if (systemOfUnits == "SI") || (systemOfUnits == "si")
        t = toSI(time)
        α = toSI(parameters[1])
        τ = toSI(parameters[2])
        if ((t.units ≠ SECOND) || (α.units ≠ DIMENSIONLESS) ||
            (τ.units ≠ SECOND))
            msg = string("Argument time and parameter τ must have units of time, while parameter α must be dimensionless.")
            throw(ErrorException(msg))
        end
    elseif (systemOfUnits == "CGS") || (systemOfUnits == "cgs")
        t = toCGS(time)
        α = toCGS(parameters[1])
        τ = toCGS(parameters[2])
        if ((t.units ≠ CGS_SECOND) || (α.units ≠ CGS_DIMENSIONLESS) ||
            (τ.units ≠ CGS_SECOND))
            msg = string("Argument time and parameter τ must have units of time, while parameter α must be dimensionless.")
            throw(ErrorException(msg))
        end
    else
        msg = "The assigned physical system of units is unknown."
         throw(ErrorException(msg))
    end

    if (t.value ≥ 0.0) && (α.value > 0.0) && (τ.value > 0.0)
        k = (α/τ) / (1 + t/τ)^(1+α.value)
    else
        msg = "Argument time must be non-negative, and parameters α and τ must be positive."
        throw(ErrorException(msg))
    end
    return ("MPL", k, τ)
end # MPL

"""
Memory function RFS (Regularized Fractional Solid)\n
    k = -E_{α,0}(-((t+δ)/τ)^α) / (E_{α,1}(-(δ/τ)^α)(t+δ)) \n
    τ = τ\n
where E_{α, β}(z) denotes the two-parameter Mittag-Leffler function. Argument `parameters` describes the tuple (α, δ, τ).
"""
function RFS(systemOfUnits::String, time::PhysicalScalar, parameters::Tuple)::Tuple
    # Ensure that a consistent system of physical units is used.
    if (systemOfUnits == "SI") || (systemOfUnits == "si")
        t = toSI(time)
        α = toSI(parameters[1])
        δ = toSI(parameters[2])
        τ = toSI(parameters[3])
        if ((t.units ≠ SECOND) || (α.units ≠ DIMENSIONLESS) ||
            (δ.units ≠ SECOND) || (τ.units ≠ SECOND))
            msg = string("Argument time and parameters τ and δ must have units of time, while parameter α must be dimensionless.")
            throw(ErrorException(msg))
        end
    elseif (systemOfUnits == "CGS") || (systemOfUnits == "cgs")
        t = toCGS(time)
        α = toCGS(parameters[1])
        δ = toCGS(parameters[2])
        τ = toCGS(parameters[3])
        if ((t.units ≠ CGS_SECOND) || (α.units ≠ CGS_DIMENSIONLESS) ||
            (δ.units ≠ CGS_SECOND) || (τ.units ≠ CGS_SECOND))
            msg = string("Argument time and parameters τ and δ must have units of time, while parameter α must be dimensionless.")
            throw(ErrorException(msg))
        end
    else
        msg = "The assigned physical system of units is unknown."
         throw(ErrorException(msg))
    end

    if α.value ≈ 1.0
        (k, τ) = SLS(systemOfUnits, time, (τ,))
    elseif ((t.value ≥ 0.0) && (α.value > 0.0) && (α.value < 1.0) &&
        (δ.value > 0.0) && (τ.value > 0.0))
        x = ((t + δ) / τ)^get(α)
        numerMLF = MittagLeffler.mittleff(get(α), 0.0, -get(x))
        y = (δ / τ)^get(α)
        denomMLF = MittagLeffler.mittleff(get(α), 1.0, -get(y))
        k = -numerMLF / (denomMLF * (δ + t))
    else
        msg = "Argument time must be non-negative, parameters τ and δ must be positive, and α ∈ (0,1]."
        throw(ErrorException(msg))
    end
    return ("RFS", k, τ)
end # RFS

"""
Memory function SLS (Standard Linear Solid)\n
    k = exp(-t/τ) / τ \n
    τ = τ\n
Argument `parameters` describes the tuple (τ,).
"""
function SLS(systemOfUnits::String, time::PhysicalScalar, parameters::Tuple)::Tuple
    # Ensure that a consistent system of physical units is used.
    if (systemOfUnits == "SI") || (systemOfUnits == "si")
        t = toSI(time)
        τ = toSI(parameters[1])
        if (t.units ≠ SECOND) || (τ.units ≠ SECOND)
            msg = string("Argument time and parameter τ must have units of time.")
            throw(ErrorException(msg))
        end
    elseif (systemOfUnits == "CGS") || (systemOfUnits == "cgs")
        t = toCGS(time)
        τ = toCGS(parameters[1])
        if (t.units ≠ CGS_SECOND) || (τ.units ≠ CGS_SECOND)
            msg = string("Argument time and parameter τ must have units of time.")
            throw(ErrorException(msg))
        end
    else
        msg = "The assigned physical system of units is unknown."
         throw(ErrorException(msg))
    end

    if (t.value ≥ 0.0) && (τ.value > 0.0)
        k = exp(-t/τ) / τ
    else
        msg = "time must be non-negative, and τ must be positive."
        throw(ErrorException(msg))
    end
    return ("SLS", k, τ)
end # SLS

#=
-------------------------------------------------------------------------------
=#

# Weights and nodes for Gaussian quadrature, which are used to determine the
# weights of quadrature for our Volterra integral equation solver.

struct GaussQuad
    S::Integer
    w::Vector{Float64}
    ξ::Vector{Float64}
end

w1 = zeros(Float64, 1)
w1[1] = 2.0
ξ1 = zeros(Float64, 1)

w2 = zeros(Float64, 2)
w2[1] = 1.0
w2[2] = 1.0
ξ2 = zeros(Float64, 2)
ξ2[1] = -sqrt(1.0/3.0)
ξ2[2] =  sqrt(1.0/3.0)

w3 = zeros(Float64, 3)
w3[1] = 5.0 / 9.0
w3[2] = 8.0 / 9.0
w3[3] = 5.0 / 9.0
ξ3 = zeros(Float64, 3)
ξ3[1] = -sqrt(3.0/5.0)
ξ3[3] =  sqrt(3.0/5.0)

w4 = zeros(Float64, 4)
w4[1] = (18.0 - sqrt(30.0)) / 36.0
w4[2] = (18.0 + sqrt(30.0)) / 36.0
w4[3] = (18.0 + sqrt(30.0)) / 36.0
w4[4] = (18.0 - sqrt(30.0)) / 36.0
ξ4 = zeros(Float64, 4)
ξ4[1] = -sqrt(3.0/7.0 + (2.0/7.0)*sqrt(6.0/5.0))
ξ4[2] = -sqrt(3.0/7.0 - (2.0/7.0)*sqrt(6.0/5.0))
ξ4[3] =  sqrt(3.0/7.0 - (2.0/7.0)*sqrt(6.0/5.0))
ξ4[4] =  sqrt(3.0/7.0 + (2.0/7.0)*sqrt(6.0/5.0))

w5 = zeros(Float64, 5)
w5[1] = (322.0 - 13.0sqrt(70.0)) / 900.0
w5[2] = (322.0 + 13.0sqrt(70.0)) / 900.0
w5[3] = 128.0 / 225.0
w5[4] = (322.0 + 13.0sqrt(70.0)) / 900.0
w5[5] = (322.0 - 13.0sqrt(70.0)) / 900.0
ξ5 = zeros(Float64, 5)
ξ5[1] = -sqrt(5.0 + 2.0sqrt(10.0/7.0)) / 3.0
ξ5[2] = -sqrt(5.0 - 2.0sqrt(10.0/7.0)) / 3.0
ξ5[4] =  sqrt(5.0 - 2.0sqrt(10.0/7.0)) / 3.0
ξ5[5] =  sqrt(5.0 + 2.0sqrt(10.0/7.0)) / 3.0

GaussQuad1 = GaussQuad(1, w1, ξ1)
GaussQuad2 = GaussQuad(2, w2, ξ2)
GaussQuad3 = GaussQuad(3, w3, ξ3)
GaussQuad4 = GaussQuad(4, w4, ξ4)
GaussQuad5 = GaussQuad(5, w5, ξ5)

# The function used to create normalized weights of quadrature, i.e., they are
# not multiplied by coefficient c. The supplied kernel may be any of those pre-
# programmed above, or a kernel of one's own creation.

"""
Function\n
    W = normalizedQuadratureWeights(systemOfUnits, N, dTime, kernel, parameters)\n
where at present `systemOfUnits` is either "SI" or "CGS". There are to be `N` intervals of size `dTime` that are to span a solution, whose `kernel` has `parameters` described via a tuple of material constants. These weights are written to a file in the user's ./files/ directory for efficient future use.\n
The supplied memory function `kernel` is to have an interface of\n
    (name, k, τ) = kernel(systemOfUnits, time, parameters)\n
where `systemOfUnits` is either "SI" or "CGS", `time` is current time, and `parameters` is a tuple containing this kernel's physical parameters, i.e., its material constants. The returned tuple contains a string specifying the `name` of the kernel being evaluated, the value of kernel `k` being evaluated at `time,` and its characteristic time `τ.`\n
The weights of quadrature returned here are normalized, e.g., the actual weights of quadrature for a linear viscoelastic kernel would be these normalized weights multiplied by a scalar coefficient of (E₀ - E∞)/E∞, which is to be supplied via a function call assigned to field `c` in an object implementing abstract type `VolterraIntegralEquation.`\n
The returned array holds normalized quadrature weights that are to be assigned to field `W` in an object implementing abstract type `VolterraIntegralEquation.`
"""
function normalizedQuadratureWeights(systemOfUnits::String, N::Integer, dTime::PhysicalScalar, kernel::Function, parameters::Tuple)::ArrayOfPhysicalScalars

    # Ensure the system of units is consistent.
    if (systemOfUnits == "SI") || (systemOfUnits == "si")
        dt = toSI(dTime)
        units = "SI"
    elseif (systemOfUnits == "CGS") || (systemOfUnits == "cgs")
        dt = toCGS(dTime)
        units = "CGS"
    else
        msg = "The assigned physical system of units is unknown."
        throw(ErrorException(msg))
    end

    # Verify the inputs.
    if N < 1
        msg = "The number of integration steps N must be positive."
        throw(ErrorException, msg)
    end
    if get(dt) < eps(Float32)
        msg = "The time step size dTime must be positive."
        throw(ErrorException, msg)
    end

    # Check to see if quadrature weights have been previously calculated or not.
    format = 'E'
    precision = 6
    aligned  = false
    (fileName, k, τ) = kernel(systemOfUnits, dt, parameters)
    if !isDimensionless(dt.units+k.units)
        msg = "Units for dTime and for the kernel are not compatible."
        throw(ErrorException, msg)
    end
    fileName = string(fileName, "_", units)
    fileName = string(fileName, "_τ=", PhysicalFields.toString(get(τ); format, precision, aligned))
    fileName = string(fileName, "_dt=", PhysicalFields.toString(get(dt); format, precision, aligned))
    fileName = string(fileName, "_N=", N)
    fileName = string(fileName, ".json")
    dirPath  = string(pwd(), "/files/")
    if !isdir(dirPath)
        mkdir(dirPath)
    end
    my_file = string(dirPath, fileName)
    if isfile(my_file)
        # The quadrature weights exist. Read them in from a file.
        json_stream = PhysicalFields.openJSONReader(dirPath, fileName)
        quadWgts = PhysicalFields.fromFile(ArrayOfPhysicalScalars, json_stream)
        PhysicalFields.closeJSONStream(json_stream)
        return quadWgts
    end

    # Determine weights of quadrature for a Volterra integral equation.
    quadWgts = ArrayOfPhysicalScalars(N, dt.units+k.units)

    tₙ = PhysicalScalar(dt.units)
    for n in 1:N
        tₙ = tₙ + dt
        # Use Gauss' 5ᵗʰ order method. Its accuracy matters.
        quad = GaussQuad5
        sum = PhysicalScalar(k.units)
        for s in 1:quad.S
            tₛ = (N - 0.5(2n - 1 + quad.ξ[s])) * dt
            (name, kₛ, τₛ) = kernel(systemOfUnits, tₛ, parameters)
            sum = sum + quad.w[s] * kₛ
        end
        quadWgts[n] = 0.5dt * sum
    end

    # Write these weights of quadrature to a file.
    json_stream = PhysicalFields.openJSONWriter(dirPath, fileName)
    PhysicalFields.toFile(quadWgts, json_stream)
    PhysicalFields.closeJSONStream(json_stream)
    println("Weights of quadrature have been written to a file.")

    return quadWgts
end # noralizedQuadratureWeights

#=
-------------------------------------------------------------------------------
=#

abstract type VolterraIntegralEquation end

# Scalar-valued Volterra integral equations of the second kind.

struct VolterraIntegralScalarEquation <: VolterraIntegralEquation
    # Dimensioning fields
    dt::PhysicalScalar          # distance separating neighboring solution nodes
    N::Int64                    # number of integration nodes in a solution path
    n::MInteger                 # current node along a solution path
    # Arrays of length N+1 holding control and response fields, and their rates
    f::ArrayOfPhysicalScalars   # array of integrated response function values
    f′::ArrayOfPhysicalScalars  # array of response function rates
    g::ArrayOfPhysicalScalars   # array of integrated control function values
    g′::ArrayOfPhysicalScalars  # array of control function rates
    t::ArrayOfPhysicalScalars   # array of times, the independent variable
    # Array of N normalized weights of quadrature for a product integral
    W::ArrayOfPhysicalScalars   # array holding the quadrature weights

    # constructors

    # For use when first creating this data structure.
    function VolterraIntegralScalarEquation(systemOfUnits::String, N::Int64, dt::PhysicalScalar, f₀::PhysicalScalar, g₀::PhysicalScalar, W::ArrayOfPhysicalScalars)

        # Ensure that a consistent system of physical units is used.
        if (systemOfUnits == "SI") || (systemOfUnits == "si")
            d𝑡 = toSI(dt)
            𝑓₀ = toSI(f₀)
            𝑔₀ = toSI(g₀)
            𝑊 = toSI(W)
        elseif (systemOfUnits == "CGS") || (systemOfUnits == "cgs")
            d𝑡 = toCGS(dt)
            𝑓₀ = toCGS(f₀)
            𝑔₀ = toCGS(g₀)
            𝑊 = toCGS(W)
        else
            msg = "The assigned physical system of units is unknown."
            throw(ErrorException(msg))
        end
        t₀ = PhysicalScalar(d𝑡.units)

        # Verify the remaining inputs.
        if N < 1
            msg = "The number of nodes N must be positive valued."
            throw(ErrorException(msg))
        end
        if 𝑓₀.units ≠ 𝑔₀.units
            msg = "Physical units for initial conditions f₀ and g₀ must be equal."
            throw(ErrorException(msg))
        end
        if !isDimensionless(𝑊)
            msg = "Weights of quadrature W must be dimensionless."
            throw(ErrorException(msg))
        end

        # Create the fields for this data structure,
        f = ArrayOfPhysicalScalars(N+1, 𝑓₀.units)
        f[1] = 𝑓₀
        f′ = ArrayOfPhysicalScalars(N+1, 𝑓₀.units-d𝑡.units)
        g  = ArrayOfPhysicalScalars(N+1, 𝑔₀.units)
        g[1] = 𝑔₀
        g′ = ArrayOfPhysicalScalars(N+1, 𝑔₀.units-d𝑡.units)
        t  = ArrayOfPhysicalScalars(N+1, t₀.units)
        t[1] = t₀
        for n in 1:N
            t[n+1] = n * d𝑡
        end
        n  = MInteger(1)

        new(d𝑡, N, n, f, f′, g, g′, t, 𝑊)
    end

    # Used by JSON3 whenever this data structure is to be created from a file.
    function VolterraIntegralScalarEquation(dt::PhysicalScalar, N::Int64, n::MInteger, f::ArrayOfPhysicalScalars, f′::ArrayOfPhysicalScalars, g::ArrayOfPhysicalScalars, g′::ArrayOfPhysicalScalars, t::ArrayOfPhysicalScalars, W::ArrayOfPhysicalScalars)

        new(dt, N, n, f, f′, g, g′, t, W)
    end
end # VolterraIntegralScalarEquation

# Methods

function Base.:(copy)(vie::VolterraIntegralScalarEquation)::VolterraIntegralScalarEquation
    dt   = copy(vie.dt)
    N    = copy(vie.N)
    n    = copy(vie.n)
    f    = copy(vie.f)
    f′   = copy(vie.f′)
    g    = copy(vie.g)
    g′   = copy(vie.g′)
    t    = copy(vie.t)
    W    = copy(vie.W)
    return VolterraIntegralScalarEquation(dt, N, n, f, f′, g, g′, t, W)
end

function Base.:(deepcopy)(vie::VolterraIntegralScalarEquation)::VolterraIntegralScalarEquation
    dt   = deepcopy(vie.dt)
    N    = deepcopy(vie.N)
    n    = deepcopy(vie.n)
    f    = deepcopy(vie.f)
    f′   = deepcopy(vie.f′)
    g    = deepcopy(vie.g)
    g′   = deepcopy(vie.g′)
    t    = deepcopy(vie.t)
    W    = deepcopy(vie.W)
    return VolterraIntegralScalarEquation(dt, N, n, f, f′, g, g′, t, W)
end

StructTypes.StructType(::Type{VolterraIntegralScalarEquation}) = StructTypes.Struct()

function toFile(vie::VolterraIntegralScalarEquation, json_stream::IOStream)
    if isopen(json_stream)
        JSON3.write(json_stream, vie)
        write(json_stream, '\n')
    else
        msg = "The supplied JSON stream is not open."
        throw(ErrorException(msg))
    end
    flush(json_stream)
    return nothing
end

function fromFile(::Type{VolterraIntegralScalarEquation}, json_stream::IOStream)::VolterraIntegralScalarEquation
    if isopen(json_stream)
        vie = JSON3.read(readline(json_stream), VolterraIntegralScalarEquation)
    else
        msg = "The supplied JSON stream is not open."
        throw(ErrorException(msg))
    end
    return vie
end

# Solver for advancing a solution step-by-step.

function advance!(vie::VolterraIntegralScalarEquation, g′ₙ::PhysicalScalar, cₙ::PhysicalScalar)
    if vie.n > vie.N
        println("The Volterra integral solution has reached its endpoint.")
        return nothing
    end

    # verify inputs
    if g′ₙ.units ≠ vie.f′.units
        msg = "Physical units for g′ₙ must equal those of vie.f′.\n"
        msg = string(msg, "   g′ₙ has units ", PhysicalFields.toString(g′ₙ.units), "\n")
        msg = string(msg, "   f′  has units ", PhysicalFields.toString(vie.f′.units))
        throw(ErrorException(msg))
    end
    if !isDimensionless(cₙ)
        msg = "Coefficient cₙ must be dimensionless scalar."
        throw(ErrorException(msg))
    end

    # update the counter
    set!(vie.n, get(vie.n)+1)

    # Solve a Volterra integral equation to get the response rate.
    n = get(vie.n)
    sum = PhysicalScalar(vie.f′.units)
    for i in 1:n-2
        sum = sum + vie.W[vie.N-i] * vie.f′[n-i]
    end
    vie.f′[n] = (g′ₙ - cₙ*sum) / (1 + cₙ*vie.W[vie.N])
    vie.g′[n] = g′ₙ

    # Integrate the differential equations governing control and response.
    if vie.n == 2
        vie.f[2] = vie.f[1] + 0.5vie.f′[2]*vie.dt
        vie.g[2] = vie.g[1] + 0.5vie.g′[2]*vie.dt
    elseif vie.n == 3
        vie.f[3] = (4/3)*vie.f[2] - (1/3)*vie.f[1] + (2/3)*vie.f′[2]*vie.dt
        vie.g[3] = (4/3)*vie.g[2] - (1/3)*vie.g[1] + (2/3)*vie.g′[2]*vie.dt
    else
        vie.f[n] = ((18/11)*vie.f[n-1] - (9/11)*vie.f[n-2] + (2/11)*vie.f[n-3]
            + (6/11)*vie.f′[n]*vie.dt)
        vie.g[n] = ((18/11)*vie.g[n-1] - (9/11)*vie.g[n-2] + (2/11)*vie.g[n-3]
            + (6/11)*vie.g′[n]*vie.dt)
    end

    return nothing
end # advance!

# Perform an iteration of refinement on a solution at current step n. Call only
# if the control function g′ₙ or coefficient cₙ undergo iterative refinement.

function update!(vie::VolterraIntegralScalarEquation, g′ₙ::PhysicalScalar, cₙ::PhysicalScalar)

    set!(vie.n, get(vie.n)-1)
    advance!(vie, g′ₙ, cₙ)
    return nothing
end # update!

#=
-------------------------------------------------------------------------------
=#

# Vector-valued Volterra integral equations of the second kind.

struct VolterraIntegralVectorEquation <: VolterraIntegralEquation
    # Dimensioning fields
    dt::PhysicalScalar          # distance separating neighboring solution nodes
    N::Int64                    # number of integration nodes in a solution path
    n::MInteger                 # current node along a solution path
    # Arrays of length N+1 holding control and response fields, and their rates
    f::ArrayOfPhysicalVectors   # array of integrated response function values
    f′::ArrayOfPhysicalVectors  # array of response function rates
    g::ArrayOfPhysicalVectors   # array of integrated control function values
    g′::ArrayOfPhysicalVectors  # array of control function rates
    t::ArrayOfPhysicalScalars   # array of times, the independent variable
    # Array of N normalized weights of quadrature for a product integral
    W::ArrayOfPhysicalScalars   # array holding the quadrature weights

    # constructors

    # For use when first creating this data structure.
    function VolterraIntegralVectorEquation(systemOfUnits::String, N::Int64, dt::PhysicalScalar, f₀::PhysicalVector, g₀::PhysicalVector, W::ArrayOfPhysicalScalars)

        # Ensure that a consistent system of physical units is used.
        if (systemOfUnits == "SI") || (systemOfUnits == "si")
            d𝑡 = toSI(dt)
            𝑓₀ = toSI(f₀)
            𝑔₀ = toSI(g₀)
            𝑊 = toSI(W)
        elseif (systemOfUnits == "CGS") || (systemOfUnits == "cgs")
            d𝑡 = toCGS(dt)
            𝑓₀ = toCGS(f₀)
            𝑔₀ = toCGS(g₀)
            𝑊 = toCGS(W)
        else
            msg = "The assigned physical system of units is unknown."
            throw(ErrorException(msg))
        end
        t₀ = PhysicalScalar(d𝑡.units)

        # Verify the remaining inputs.
        if N < 1
            msg = "The number of nodes N must be positive valued."
            throw(ErrorException(msg))
        end
        if 𝑓₀.units ≠ 𝑔₀.units
            msg = "Physical units for initial conditions f₀ and g₀ must be equal."
            throw(ErrorException(msg))
        end
        if 𝑓₀.vector.len ≠ 𝑔₀.vector.len
            msg = "Length of initial conditions, vectors f₀ and g₀, must equal."
            throw(ErrorException(msg))
        end
        if !isDimensionless(𝑊)
            msg = "Weights of quadrature W must be dimensionless."
            throw(ErrorException(msg))
        end

        # Create the fields for this data structure,
        f = ArrayOfPhysicalVectors(N+1, 𝑓₀.vector.len, 𝑓₀.units)
        f[1] = 𝑓₀
        f′ = ArrayOfPhysicalVectors(N+1, 𝑓₀.vector.len, 𝑓₀.units-d𝑡.units)
        g  = ArrayOfPhysicalVectors(N+1, 𝑔₀.vector.len, 𝑔₀.units)
        g[1] = 𝑔₀
        g′ = ArrayOfPhysicalVectors(N+1, 𝑔₀.vector.len, 𝑔₀.units-d𝑡.units)
        t  = ArrayOfPhysicalScalars(N+1, t₀.units)
        t[1] = t₀
        for n in 1:N
            t[n+1] = n * d𝑡
        end
        n  = MInteger(1)

        new(d𝑡, N, n, f, f′, g, g′, t, 𝑊)
    end

    # Used by JSON3 whenever this data structure is to be created from a file.
    function VolterraIntegralVectorEquation(dt::PhysicalScalar, N::Int64, n::MInteger, f::ArrayOfPhysicalVectors, f′::ArrayOfPhysicalVectors, g::ArrayOfPhysicalVectors, g′::ArrayOfPhysicalVectors, t::ArrayOfPhysicalScalars, W::ArrayOfPhysicalScalars)

        new(dt, N, n, f, f′, g, g′, t, W)
    end
end # VolterraIntegralVectorEquation

# Methods

function Base.:(copy)(vie::VolterraIntegralVectorEquation)::VolterraIntegralVectorEquation
    dt   = copy(vie.dt)
    N    = copy(vie.N)
    n    = copy(vie.n)
    f    = copy(vie.f)
    f′   = copy(vie.f′)
    g    = copy(vie.g)
    g′   = copy(vie.g′)
    t    = copy(vie.t)
    W    = copy(vie.W)
    return VolterraIntegralVectorEquation(dt, N, n, f, f′, g, g′, t, W)
end

function Base.:(deepcopy)(vie::VolterraIntegralVectorEquation)::VolterraIntegralVectorEquation
    dt   = deepcopy(vie.dt)
    N    = deepcopy(vie.N)
    n    = deepcopy(vie.n)
    f    = deepcopy(vie.f)
    f′   = deepcopy(vie.f′)
    g    = deepcopy(vie.g)
    g′   = deepcopy(vie.g′)
    t    = deepcopy(vie.t)
    W    = deepcopy(vie.W)
    return VolterraIntegralVectorEquation(dt, N, n, f, f′, g, g′, t, W)
end

StructTypes.StructType(::Type{VolterraIntegralVectorEquation}) = StructTypes.Struct()

function toFile(vie::VolterraIntegralVectorEquation, json_stream::IOStream)
    if isopen(json_stream)
        JSON3.write(json_stream, vie)
        write(json_stream, '\n')
    else
        msg = "The supplied JSON stream is not open."
        throw(ErrorException(msg))
    end
    flush(json_stream)
    return nothing
end

function fromFile(::Type{VolterraIntegralVectorEquation}, json_stream::IOStream)::VolterraIntegralVectorEquation
    if isopen(json_stream)
        vie = JSON3.read(readline(json_stream), VolterraIntegralVectorEquation)
    else
        msg = "The supplied JSON stream is not open."
        throw(ErrorException(msg))
    end
    return vie
end

# Solver for advancing a solution step-by-step.

function advance!(vie::VolterraIntegralVectorEquation, g′ₙ::PhysicalVector, cₙ::PhysicalScalar)
    if vie.n > vie.N
        println("The Volterra integral solution has reached its endpoint.")
        return nothing
    end

    # verify inputs
    if g′ₙ.units ≠ vie.f′.units
        msg = "Physical units for g′ₙ must equal those of vie.f′.\n"
        msg = string(msg, "   g′ₙ has units ", PhysicalFields.toString(g′ₙ.units), "\n")
        msg = string(msg, "   f′  has units ", PhysicalFields.toString(vie.f′.units))
        throw(ErrorException(msg))
    end
    if g′ₙ.vector.len ≠ g′.array.cols
        msg = "Vector g′ₙ has the wrong length."
        throw(ErrorException(msg))
    end
    if !isDimensionless(cₙ)
        msg = "Coefficient cₙ must be dimensionless scalar."
        throw(ErrorException(msg))
    end

    # update the counter
    set!(vie.n, get(vie.n)+1)

    # Solve a Volterra integral equation to get the response rate.
    n = get(vie.n)
    len = vie.f′.vector.len
    sum = PhysicalVector(len, vie.f′.units)
    for i in 1:n-2
        for j in 1:len
            sum[j] = sum[j] + vie.W[vie.N-i] * vie.f′[n-i,j]
        end
    end
    for j in 1:len
        vie.f′[n,j] = (g′ₙ[j] - cₙ*sum[j]) / (1 + cₙ*vie.W[vie.N])
        vie.g′[n,j] = g′ₙ[j]
    end

    # Integrate the differential equations governing control and response.
    for j in 1:len
        if vie.n == 2
            vie.f[2,j] = vie.f[1,j] + 0.5vie.f′[2,j]*vie.dt
            vie.g[2,j] = vie.g[1,j] + 0.5vie.g′[2,j]*vie.dt
        elseif vie.n == 3
            vie.f[3,j] = ((4/3)*vie.f[2,j] - (1/3)*vie.f[1,j]
                + (2/3)*vie.f′[2,j]*vie.dt)
            vie.g[3,j] = ((4/3)*vie.g[2,j] - (1/3)*vie.g[1,j]
                + (2/3)*vie.g′[2,j]*vie.dt)
        else
            vie.f[n,j] = ((18/11)*vie.f[n-1,j] - (9/11)*vie.f[n-2,j]
                + (2/11)*vie.f[n-3,j] + (6/11)*vie.f′[n,j]*vie.dt)
            vie.g[n,j] = ((18/11)*vie.g[n-1,j] - (9/11)*vie.g[n-2,j]
                + (2/11)*vie.g[n-3,j] + (6/11)*vie.g′[n,j]*vie.dt)
        end
    end

    return nothing
end # advance!

# Perform an iteration of refinement on a solution at current step n. Call only
# if the control function g′ₙ or coefficient cₙ undergo iterative refinement.

function update!(vie::VolterraIntegralVectorEquation, g′ₙ::PhysicalVector, cₙ::PhysicalScalar)

    set!(vie.n, get(vie.n)-1)
    advance!(vie, g′ₙ, cₙ)
    return nothing
end # update!

#=
-------------------------------------------------------------------------------
=#

# Tensor-valued Volterra integral equations of the second kind.

struct VolterraIntegralTensorEquation <: VolterraIntegralEquation
    # Dimensioning fields
    dt::PhysicalScalar          # distance separating neighboring solution nodes
    N::Int64                    # number of integration nodes in a solution path
    n::MInteger                 # current node along a solution path
    # Arrays of length N+1 holding control and response fields, and their rates
    f::ArrayOfPhysicalTensors   # array of integrated response function values
    f′::ArrayOfPhysicalTensors  # array of response function rates
    g::ArrayOfPhysicalTensors   # array of integrated control function values
    g′::ArrayOfPhysicalTensors  # array of control function rates
    t::ArrayOfPhysicalScalars   # array of times, the independent variable
    # Array of N normalized weights of quadrature for a product integral
    W::ArrayOfPhysicalScalars   # array holding the quadrature weights

    # constructors

    # For use when first creating this data structure.
    function VolterraIntegralTensorEquation(systemOfUnits::String, N::Int64, dt::PhysicalScalar, f₀::PhysicalTensor, g₀::PhysicalTensor, W::ArrayOfPhysicalScalars)

        # Ensure that a consistent system of physical units is used.
        if (systemOfUnits == "SI") || (systemOfUnits == "si")
            d𝑡 = toSI(dt)
            𝑓₀ = toSI(f₀)
            𝑔₀ = toSI(g₀)
            𝑊 = toSI(W)
        elseif (systemOfUnits == "CGS") || (systemOfUnits == "cgs")
            d𝑡 = toCGS(dt)
            𝑓₀ = toCGS(f₀)
            𝑔₀ = toCGS(g₀)
            𝑊 = toCGS(W)
        else
            msg = "The assigned physical system of units is unknown."
            throw(ErrorException(msg))
        end
        t₀ = PhysicalScalar(d𝑡.units)

        # Verify the remaining inputs.
        if N < 1
            msg = "The number of nodes N must be positive valued."
            throw(ErrorException(msg))
        end
        if 𝑓₀.units ≠ 𝑔₀.units
            msg = "Physical units for initial conditions f₀ and g₀ must be equal."
            throw(ErrorException(msg))
        end
        if (𝑓₀.matrix.rows ≠ 𝑔₀.matrix.rows) || (𝑓₀.matrix.cols ≠ 𝑔₀.matrix.cols)
            msg = "Dimensions of initial conditions, tensors f₀ and g₀, must equal."
            throw(ErrorException(msg))
        end
        if !isDimensionless(𝑊)
            msg = "Weights of quadrature W must be dimensionless."
            throw(ErrorException(msg))
        end

        # Create the fields for this data structure,
        rows = 𝑓₀.matrix.rows
        cols = 𝑓₀.matrix.cols
        f  = ArrayOfPhysicalTensors(N+1, rows, cols, 𝑓₀.units)
        f[1] = 𝑓₀
        f′ = ArrayOfPhysicalTensors(N+1, rows, cols, 𝑓₀.units-d𝑡.units)
        g  = ArrayOfPhysicalTensors(N+1, rows, cols, 𝑔₀.units)
        g[1] = 𝑔₀
        g′ = ArrayOfPhysicalTensors(N+1, rows, cols, 𝑔₀.units-d𝑡.units)
        t  = ArrayOfPhysicalScalars(N+1, t₀.units)
        t[1] = t₀
        for n in 1:N
            t[n+1] = n * d𝑡
        end
        n  = MInteger(1)

        new(d𝑡, N, n, f, f′, g, g′, t, 𝑊)
    end

    # Used by JSON3 whenever this data structure is to be created from a file.
    function VolterraIntegralTensorEquation(dt::PhysicalScalar, N::Int64, n::MInteger, f::ArrayOfPhysicalTensors, f′::ArrayOfPhysicalTensors, g::ArrayOfPhysicalTensors, g′::ArrayOfPhysicalTensors, t::ArrayOfPhysicalScalars, W::ArrayOfPhysicalScalars)

        new(dt, N, n, f, f′, g, g′, t, W)
    end
end # VolterraIntegralTensorEquation

# Methods

function Base.:(copy)(vie::VolterraIntegralTensorEquation)::VolterraIntegralTensorEquation
    dt   = copy(vie.dt)
    N    = copy(vie.N)
    n    = copy(vie.n)
    f    = copy(vie.f)
    f′   = copy(vie.f′)
    g    = copy(vie.g)
    g′   = copy(vie.g′)
    t    = copy(vie.t)
    W    = copy(vie.W)
    return VolterraIntegralTensorEquation(dt, N, n, f, f′, g, g′, t, W)
end

function Base.:(deepcopy)(vie::VolterraIntegralTensorEquation)::VolterraIntegralTensorEquation
    dt   = deepcopy(vie.dt)
    N    = deepcopy(vie.N)
    n    = deepcopy(vie.n)
    f    = deepcopy(vie.f)
    f′   = deepcopy(vie.f′)
    g    = deepcopy(vie.g)
    g′   = deepcopy(vie.g′)
    t    = deepcopy(vie.t)
    W    = deepcopy(vie.W)
    return VolterraIntegralTensorEquation(dt, N, n, f, f′, g, g′, t, W)
end

StructTypes.StructType(::Type{VolterraIntegralTensorEquation}) = StructTypes.Struct()

function toFile(vie::VolterraIntegralTensorEquation, json_stream::IOStream)
    if isopen(json_stream)
        JSON3.write(json_stream, vie)
        write(json_stream, '\n')
    else
        msg = "The supplied JSON stream is not open."
        throw(ErrorException(msg))
    end
    flush(json_stream)
    return nothing
end

function fromFile(::Type{VolterraIntegralTensorEquation}, json_stream::IOStream)::VolterraIntegralTensorEquation
    if isopen(json_stream)
        vie = JSON3.read(readline(json_stream), VolterraIntegralTensorEquation)
    else
        msg = "The supplied JSON stream is not open."
        throw(ErrorException(msg))
    end
    return vie
end

# Solver for advancing a solution step-by-step.

function advance!(vie::VolterraIntegralTensorEquation, g′ₙ::PhysicalTensor, cₙ::PhysicalScalar)
    if vie.n > vie.N
        println("The Volterra integral solution has reached its endpoint.")
        return nothing
    end

    # verify inputs
    if g′ₙ.units ≠ vie.f′.units
        msg = "Physical units for g′ₙ must equal those of vie.f′.\n"
        msg = string(msg, "   g′ₙ has units ", PhysicalFields.toString(g′ₙ.units), "\n")
        msg = string(msg, "   f′  has units ", PhysicalFields.toString(vie.f′.units))
        throw(ErrorException(msg))
    end
    if (g′ₙ.matrix.rows ≠ g′.array.rows) || (g′ₙ.matrix.cols ≠ g′.array.cols)
        msg = "Tensor g′ₙ has the wrong dimensions."
        throw(ErrorException(msg))
    end
    if !isDimensionless(cₙ)
        msg = "Coefficient cₙ must be dimensionless scalar."
        throw(ErrorException(msg))
    end

    # update the counter
    set!(vie.n, get(vie.n)+1)

    # Solve a Volterra integral equation to get the response rate.
    n = get(vie.n)
    rows = vie.f′.matrix.rows
    cols = vie.f′.matrix.cols
    sum  = PhysicalTensor(rows, cols, vie.f′.units)
    for i in 1:n-2
        for j in 1:rows
            for k in 1:cols
                sum[j,k] = sum[j,k] + vie.W[vie.N-i] * vie.f′[n-i,j,k]
            end
        end
    end
    for j in 1:rows
        for k in 1:cols
            vie.f′[n,j,k] = (g′ₙ[j,k] - cₙ*sum[j,k]) / (1 + cₙ*vie.W[vie.N])
            vie.g′[n,j,k] = g′ₙ[j,k]
        end
    end

    # Integrate the differential equations governing the control and response.
    for j in 1:rows
        for k in 1:cols
            if vie.n == 2
                vie.f[2,j,k] = vie.f[1,j,k] + 0.5vie.f′[2,j,k]*vie.dt
                vie.g[2,j,k] = vie.g[1,j,k] + 0.5vie.g′[2,j,k]*vie.dt
            elseif vie.n == 3
                vie.f[3,j,k] = ((4/3)*vie.f[2,j,k] - (1/3)*vie.f[1,j,k]
                    + (2/3)*vie.f′[2,j,k]*vie.dt)
                vie.g[3,j,k] = ((4/3)*vie.g[2,j,k] - (1/3)*vie.g[1,j,k]
                    + (2/3)*vie.g′[2,j,k]*vie.dt)
            else
                vie.f[n,j,k] = ((18/11)*vie.f[n-1,j,k] - (9/11)*vie.f[n-2,j,k]
                    + (2/11)*vie.f[n-3,j,k] + (6/11)*vie.f′[n,j,k]*vie.dt)
                vie.g[n,j,k] = ((18/11)*vie.g[n-1,j,k] - (9/11)*vie.g[n-2,j,k]
                    + (2/11)*vie.g[n-3,j,k] + (6/11)*vie.g′[n,j,k]*vie.dt)
            end
        end
    end

    return nothing
end # advance!

# Perform an iteration of refinement on a solution at current step n. Call only
# if the control function g′ₙ or coefficient cₙ undergo iterative refinement.

function update!(vie::VolterraIntegralTensorEquation, g′ₙ::PhysicalTensor, cₙ::PhysicalScalar)

    set!(vie.n, get(vie.n)-1)
    advance!(vie, g′ₙ, cₙ)
    return nothing
end # update!

end  # module VolterraIntegralEquations
