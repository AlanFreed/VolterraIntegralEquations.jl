#=
Created on Tue 27 Jun 2023
Updated on Wed 31 Jan 2024
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
    #   c  is a scalar function, e.g., (E₀ - E∞)/E∞ in viscoelasticity
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
    (k, τ) = <memoryFunctionName>(systemOfUnits, time, parameters)
which returns a tuple whose first entry is the value of the memory function `k` and whose second entry is its controlling characteristic time `τ`, which is the smallest one whenever multiple characteristic times are present. Here argument `systemOfUnits` is either "SI" or "CGS", argument `time` is current time, and argument `parameters` is a tuple containing this kernel's physical parameters.

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
wherein τ denotes a characteristic time for creep. There are two in the BOX model, and n in the MCM, arranged so that 0 < τ₁ < τ₂ < ⋯ < τₙ, with each cᵢ > 0, i = 1, 2, …, n, being a coefficient in the Prony series whose collective sum is 1, i.e., ∑_{i=1}^n cᵢ = 1. Parameter α is the exponent in a power law, and parameter δ is a shift in time introduced to remove a weak singularity.

The following memory functions are weakly singular at the upper limit of integration in their Volterra integrals: CCM, FLS and KWW.
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
    return (k, τ₁)
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
    return (k, τ)
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
    return (k, τ)
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
    return (k, τ)
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
    return (k, τ₁)
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
    return (k, τ)
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
    return (k, τ)
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
    return (k, τ)
end # SLS

#=
-------------------------------------------------------------------------------
=#

# The function used to create normalized weights of quadrature, i.e., they are
# not multiplied by coefficient c. The supplied kernel may be any of those 
# pre-programmed above, or a kernel of one's own creation.

"""
Function\n
    W = normalizedQuadratureWeights(systemOfUnits, dTime, parameters, kernel, significantFigures)\n
where `systemOfUnits` is either "SI" or "CGS", `dTime` is an uniform increment in time separating nodes from their nearest neighbors, `parameters` is a tuple of material constants to be passed to the memory function `kernel`. The array for weights of quadrature is truncated at a specified number of `significantFigures` in accuracy, which default to 5, but accept values from 2 through 10.\n
The supplied memory function `kernel` is to have an interface of\n
    (k, τ) = kernel(systemOfUnits, time, parameters)\n
where `systemOfUnits` is either "SI" or "CGS", `time` is current time, and `parameters` is a tuple containing this kernel's physical parameters, i.e., its material constants. The returned tuple contains values for the kernel `k` and its characteristic time `τ.`\n
The weights of quadrature returned here are normalized, e.g., the actual weights of quadrature for a viscoelastic kernel would be these normalized weights multiplied by a scalar coefficient of (E₀ - E∞)/E∞, which is to be returned from a function assigned to field `c` in an object implementing abstract type `VolterraIntegralEquation.`\n
The returned array holds Nₘₐₓ normalized quadrature weights that is to be assigned to field `W` in an object implementing abstract type `VolterraIntegralEquation.`
"""
function normalizedQuadratureWeights(systemOfUnits::String, dTime::PhysicalScalar, parameters::Tuple, kernel::Function, Nₘₐₓ::Integer, significantFigures::Integer=5)::ArrayOfPhysicalTensors

    # Ensure the system of units is consistent.
    if (systemOfUnits == "SI") || (systemOfUnits == "si")
        dt = toSI(dTime)
    elseif (systemOfUnits == "CGS") || (systemOfUnits == "cgs")
        dt = toCGS(dTime)
    else
        msg = "The assigned physical system of units is unknown."
        throw(ErrorException(msg))
    end

    # Determine the truncation length for the array of weights.
    (k, τ) = kernel(systemOfUnits, dt, parameters)
    L = ceil(τ/dt)
    if L < 10
        msg = string("WARNING: There are ", Int64(get(L)), " integration steps per unit\n")
        msg = string(msg, "characteristic time. There should be at least 10.")
        println(msg)
    end
    if significantFigures < 2
        SF = 2
    elseif significantFigures > 10
        SF = 10
    else
        SF = significantFigures
    end
    N = 1
    (k, τ) = kernel(systemOfUnits, N*L*dt, parameters)
    while get(k) > 10.0^(-SF)
        N = N + 1
        (k, τ) = kernel(systemOfUnits, N*L*dt, parameters)
    end
    if N * L < Nₘₐₓ
        Nₘₐₓ = N * L
    end

    # Basic arrays needed to create the weights of quadrature.

    # Inverse of the 3x3 midpoint Vandermonde matrix X, i.e., Xinv.
    Xinv = zeros(Float64, 3, 3)
    Xinv[1,1] =  0.0
    Xinv[1,2] = -0.5
    Xinv[1,3] =  0.5
    Xinv[2,1] =  1.0
    Xinv[2,2] =  0.0
    Xinv[2,3] = -1.0
    Xinv[3,1] =  0.0
    Xinv[3,2] =  0.5
    Xinv[3,3] =  0.5

    # The three Gauss-quadrature matrices used to create moment μ₁.
    m = zeros(Float64, 3, 3, 3)  # for indices i, j, s in what follows
    m[1,1,1] = 1.0
    m[1,2,1] = 1.0
    m[1,3,1] = 1.0
    m[2,1,1] = -sqrt(3.0/80.0)
    m[2,2,1] = -sqrt(27.0/80.0)
    m[2,3,1] = -sqrt(15.0/16.0)
    m[3,1,1] = 3.0 / 80.0
    m[3,2,1] = 27.0 / 80.0
    m[3,3,1] = 15.0 / 16.0

    m[1,1,2] = 1.0
    m[1,2,2] = 1.0
    m[1,3,2] = 1.0
    m[2,1,2] = 0.0
    m[2,2,2] = 0.0
    m[2,3,2] = 0.0
    m[3,1,2] = 0.0
    m[3,2,2] = 0.0
    m[3,3,2] = 0.0

    m[1,1,3] = 1.0
    m[1,2,3] = 1.0
    m[1,3,3] = 1.0
    m[2,1,3] = sqrt(3.0/80.0)
    m[2,2,3] = sqrt(27.0/80.0)
    m[2,3,3] = sqrt(15.0/16.0)
    m[3,1,3] = 3.0 / 80.0
    m[3,2,3] = 27.0 / 80.0
    m[3,3,3] = 15.0 / 16.0

    # The three Gauss-quadrature vectors used to create moments μᵢ, i > 1.
    v = zeros(Float64, 3,3)  # for indices i, s in what follows
    v[1,1] = 1.0
    v[2,1] = -sqrt(27.0/20.0)
    v[3,1] = 27.0 / 20.0

    v[1,2] = 1.0
    v[2,2] = 0.0
    v[3,2] = 0.0

    v[1,3] = 1.0
    v[2,3] = sqrt(27.0/20.0)
    v[3,3] = 27.0 / 20.0

    # The weights and nodes of Gaussian quadrature.
    w = zeros(Float64, 3)
    w[1] = 5.0 / 9.0
    w[2] = 8.0 / 9.0
    w[3] = 5.0 / 9.0

    x = zeros(Float64, 3)
    x[1] = -sqrt(3.0/5.0)
    x[2] = 0.0
    x[3] = sqrt(3.0/5.0)

    #=
    ----------------------------------------------------------------------------
    =#

    # Create the first moment matrix.

    (k, τ) = kernel(systemOfUnits, dt, parameters)  # call made to establish units.
    μ₁ = PhysicalTensor(3, 3, dt.units+k.units)
    for j in 1:3
        coef = (j - 0.5)*dt / 6
        for i in 1:3
            sum = PhysicalScalar(k.units)
            for s in 1:3
                t = (j - 0.5) * (1 - x[s]) * dt / 6
                (k, τ) = kernel(systemOfUnits, t, parameters)
                sum = sum + w[s] * m[i,j,s] * k
            end
            μ₁[i,j] = coef * sum
        end
    end
    quadMatrix = PhysicalTensor(3, 3, μ₁.units)
    for i in 1:3
        for j in 1:3
            sum = PhysicalScalar(μ₁.units)
            for s in 1:3
                sum = sum + Xinv[i,s] * μ₁[s,j]
            end
            quadMatrix[i,j] = sum
        end
    end
    quadWgts = ArrayOfPhysicalTensors(Nₘₐₓ, 3, 3, μ₁.units)
    quadWgts[1] = quadMatrix

    # Create the remaining moment matrices.
    coef = dt / 2
    for n in 2:Nₘₐₓ
        μₙ = PhysicalTensor(3, 3, μ₁.units)
        for i in 1:3
            for j in 1:3
                sum = PhysicalScalar(k.units)
                for s in 1:3
                    t = (n - (5 - j)/3 - x[s]/2) * dt
                    (k, τ) = kernel(systemOfUnits, t, parameters)
                    sum = sum + w[s] * v[i,s] * k
                end
                μₙ[i,j] = coef * sum
            end
        end
        quadMatrix = PhysicalTensor(3, 3, μₙ.units)
        for i in 1:3
            for j in 1:3
                sum = PhysicalScalar(μₙ.units)
                for s in 1:3
                    sum = sum + Xinv[i,s] * μₙ[s,j]
                end
                quadMatrix[i,j] = sum
            end
        end
        quadWgts[n] = quadMatrix
    end

    return quadWgts
end # noralizedQuadratureWeights

#=
-------------------------------------------------------------------------------
=#

abstract type VolterraIntegralEquation end

# Scalar-valued Volterra integral equations of the second kind.

struct VolterraIntegralScalarEquation <: VolterraIntegralEquation
    # Dimensioning fields
    dt::PhysicalScalar          # distance separating global integration nodes
    N::Integer                  # number of integration nodes in a solution path
    Nₘₐₓ::Integer               # maximum number of nodes whose history is kept
    n::MInteger                 # current node along a solution path
    # Arrays of length N+1 holding integrated variable rates at the global nodes
    f::ArrayOfPhysicalScalars   # array of integrated response function values
    g::ArrayOfPhysicalScalars   # array of integrated control function values
    t::ArrayOfPhysicalScalars   # array of times, the independent variable
    # Array of length 3N holding response function rates at over local intervals
    f′::ArrayOfPhysicalScalars  # history array of response function rates
    # Array of Nₘₐₓ normalized weights of quadrature for a product integral
    W::ArrayOfPhysicalTensors   # array of matrices holding quadrature weights

    # constructors

    # For use when first creating this data structure.
    function VolterraIntegralScalarEquation(systemOfUnits::String, N::Integer, dt::PhysicalScalar, f₀::PhysicalScalar, g₀::PhysicalScalar, W::ArrayOfPhysicalTensors)

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
        if !isDimensionless(𝑊) || (𝑊.array.rows ≠ 3) || (𝑊.array.cols ≠ 3)
            msg = "Weights of quadrature W must be dimensionless 3x3 matrices."
            throw(ErrorException(msg))
        end

        # Create the fields for this data structure,
        Nₘₐₓ = 𝑊.array.pgs
        f = ArrayOfPhysicalScalars(N+1, 𝑓₀.units)
        f[1] = 𝑓₀
        g = ArrayOfPhysicalScalars(N+1, 𝑔₀.units)
        g[1] = 𝑔₀
        t = ArrayOfPhysicalScalars(N+1, t₀.units)
        t[1] = t₀
        for n in 1:N
            t[n+1] = n * d𝑡
        end
        n  = MInteger(1)
        f′ = ArrayOfPhysicalScalars(3N, 𝑓₀.units-d𝑡.units)

        new(d𝑡, N, Nₘₐₓ, n, f, g, t, f′, 𝑊)
    end

    # Used by JSON3 whenever this data structure is to be created from a file.
    function VolterraIntegralScalarEquation(dt::PhysicalScalar, N::Integer, Nₘₐₓ::Integer, n::MInteger, f::ArrayOfPhysicalScalars, g::ArrayOfPhysicalScalars, t::ArrayOfPhysicalScalars, f′::ArrayOfPhysicalScalars, W::ArrayOfPhysicalTensors)

        new(dt, N, Nₘₐₓ, n, f, g, t, f′, W)
    end
end # VolterraIntegralScalarEquation

# Methods

function Base.:(copy)(vie::VolterraIntegralScalarEquation)::VolterraIntegralScalarEquation
    dt   = copy(vie.dt)
    N    = copy(vie.N)
    Nₘₐₓ = copy(vie.Nₘₐₓ)
    n    = copy(vie.n)
    t    = copy(vie.t)
    f′   = copy(vie.f′)
    W    = copy(vie.W)
    return VolterraIntegralScalarEquation(dt, N, Nₘₐₓ, n, f, g, t, f′, W)
end

function Base.:(deepcopy)(vie::VolterraIntegralScalarEquation)::VolterraIntegralScalarEquation
    dt   = deepcopy(vie.dt)
    N    = deepcopy(vie.N)
    Nₘₐₓ = deepcopy(vie.Nₘₐₓ)
    n    = deepcopy(vie.n)
    t    = deepcopy(vie.t)
    f′   = deepcopy(vie.f′)
    W    = deepcopy(vie.W)
    return VolterraIntegralScalarEquation(dt, N, Nₘₐₓ, n, f, g, t, f′, W)
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

function advance!(vie::VolterraIntegralScalarEquation, g′ₙ::ArrayOfPhysicalScalars, cₙ::ArrayOfPhysicalScalars)
    if vie.n > vie.N
        println("The Volterra integral solution has reached its endpoint.")
        return nothing
    end

    # verify inputs
    if (g′ₙ.array.len ≠ 3) || (cₙ.array.len ≠ 3)
        msg = "The control function g′ₙ and coefficient cₙ must be of lengths 3."
        throw(ErrorException(msg))
    end
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

    # Create the matrix coefficient for the right-hand side (rhs) product.
    I = PhysicalTensor(3, 3, vie.W.units)
    for i in 1:3
        I[i,i] = PhysicalScalar(1.0, vie.W.units)
    end
    W₁ = vie.W[1]
    W₁inv = inv(I + cₙ*transpose(W₁))

    # Create the temporary working arrays, which are of length 3.
    b′ = ArrayOfPhysicalScalars(3, g′ₙ.units)
    f′ = ArrayOfPhysicalScalars(3, g′ₙ.units)
    x′ = ArrayOfPhysicalScalars(3, g′ₙ.units)

    # Create the vector for the rhs product.

    # First, add in the control contribution to the rhs vector.
    for i in 1:3
        b′[i] = g′ₙ[i]
    end

    # Second, incorporate the history effects acting on this rhs vector.
    if vie.n ≤ vie.Nₘₐₓ
        # Advance the solution along a path with full history.
        for m in 1:vie.n-1
            W = vie.W[vie.n-m+1]
            f′[1] = vie.f′[3(m-1)+1]
            f′[2] = vie.f′[3(m-1)+2]
            f′[3] = vie.f′[3(m-1)+3]
            for i in 1:3
                for j in 1:3
                    b′[i] = b′[i] - cₙ[i] * W[j,i] * f′[j]
                end
            end
        end
    else  # vie.n > vie.Nₘₐₓ
        # Advance the solution along a path with truncated history.
        for m in 1:vie.Nₘₐₓ-1
            W = vie.W[vie.Nₘₐₓ-m+1]
            f′[1] = vie.f′[3(m+vie.n-vie.Nₘₐₓ-1)+1]
            f′[2] = vie.f′[3(m+vie.n-vie.Nₘₐₓ-1)+2]
            f′[3] = vie.f′[3(m+vie.n-vie.Nₘₐₓ-1)+3]
            for i in 1:3
                for j in 1:3
                    b′[i] = b′[i] - cₙ[i] * W[j,i] * f′[j]
                end
            end
        end
    end

    # Finally, solve A x′ = b′ for x′, i.e., solve the linear system.
    for i in 1:3
        for j in 1:3
            x′[i] = x′[i] + W₁inv[i,j] * b′[j]
        end
    end

    # Assign this solution to its location for nᵗʰ step in the history vector.
    for i in 1:3
        vie.f′[3(vie.n-1)+i] = x′[i]
    end

    n = get(vie.n)
    # Integrate rate expressions describing the control and response functions.
    vie.f[n+1] = (vie.f[n] + (vie.dt/8) *
        (3vie.f′[3(n-1)+1] + 2vie.f′[3(n-1)+2] + 3vie.f′[3(n-1)+3]))
    vie.g[n+1] = vie.g[n] + (vie.dt/8) * (3g′ₙ[1] + 2g′ₙ[2] + 3g′ₙ[3])

    # Update the counter.
    set!(vie.n, n+1)

    return nothing
end # advance!

# Perform an iteration of refinement on a solution at current step n. Call only
# if the control function g′ₙ or coefficient cₙ undergo iterative refinement.

function update!(vie::VolterraIntegralScalarEquation, g′ₙ::ArrayOfPhysicalScalars, cₙ::ArrayOfPhysicalScalars)

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
    dt::PhysicalScalar          # distance separating global integration nodes
    N::Integer                  # number of integration nodes in a solution path
    Nₘₐₓ::Integer               # maximum number of nodes whose history is kept
    n::MInteger                 # current node along a solution path
    # Arrays of length N+1 holding integrated variable rates at the global nodes
    f::ArrayOfPhysicalVectors   # array of integrated response function values
    g::ArrayOfPhysicalVectors   # array of integrated control function values
    t::ArrayOfPhysicalScalars   # array of times, the independent variable
    # Array of length 3N holding response function rates at over local intervalS
    f′::ArrayOfPhysicalVectors  # history array of response function rates
    # Array of Nₘₐₓ normalized weights of quadrature for a product integral
    W::ArrayOfPhysicalTensors   # array of matrices holding quadrature weights

    # constructors

    # For use when first creating this data structure.
    function VolterraIntegralVectorEquation(systemOfUnits::String, N::Integer, dt::PhysicalScalar, f₀::PhysicalVector, g₀::PhysicalVector, W::ArrayOfPhysicalTensors)

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
        if (𝑓₀.units ≠ 𝑔₀.units) || (𝑓₀.vector.len ≠ 𝑔₀.vector.len)
            msg = "Units and dimensions for initial conditions f₀ and g₀ must be equal."
            throw(ErrorException(msg))
        end
        if !isDimensionless(𝑊) || (𝑊.array.rows ≠ 3) || (𝑊.array.cols ≠ 3)
            msg = "Weights of quadrature W must be dimensionless 3x3 matrices."
            throw(ErrorException(msg))
        end

        # Create the fields for this data structure,
        Nₘₐₓ = 𝑊.array.pgs
        f = ArrayOfPhysicalVectors(N+1, 𝑓₀.vector.len, 𝑓₀.units)
        f[1] = 𝑓₀
        g = ArrayOfPhysicalVectors(N+1, 𝑔₀.vector.len, 𝑔₀.units)
        g[1] = 𝑔₀
        t = ArrayOfPhysicalScalars(N+1, t₀.units)
        t[1] = t₀
        for n in 1:N
            t[n+1] = n * d𝑡
        end
        f′ = ArrayOfPhysicalVectors(3N, 𝑓₀.vector.len, 𝑓₀.units-d𝑡.units)
        n  = MInteger(1)
        new(d𝑡, N, Nₘₐₓ, n, f, g, t, f′, 𝑊)
    end

    # Used by JSON3 whenever this data structure is to be created from a file.
    function VolterraIntegralVectorEquation(dt::PhysicalScalar, N::Integer, Nₘₐₓ::Integer, n::MInteger, f::ArrayOfPhysicalVectors, g::ArrayOfPhysicalVectors, t::ArrayOfPhysicalScalars, f′::ArrayOfPhysicalVectors, W::ArrayOfPhysicalTensors)
        new(dt, N, Nₘₐₓ, n, f, g, t, f′, W)
    end
end # VolterraIntegralVectorEquation

# Methods

function Base.:(copy)(vie::VolterraIntegralVectorEquation)::VolterraIntegralVectorEquation
    dt   = copy(vie.dt)
    N    = copy(vie.N)
    Nₘₐₓ = copy(vie.Nₘₐₓ)
    n    = copy(vie.n)
    t    = copy(vie.t)
    f′   = copy(vie.f′)
    W    = copy(vie.W)
    return VolterraIntegralVectorEquation(dt, N, Nₘₐₓ, n, f, g, t, f′, W)
end

function Base.:(deepcopy)(vie::VolterraIntegralVectorEquation)::VolterraIntegralVectorEquation
    dt   = deepcopy(vie.dt)
    N    = deepcopy(vie.N)
    Nₘₐₓ = deepcopy(vie.Nₘₐₓ)
    n    = deepcopy(vie.n)
    t    = deepcopy(vie.t)
    f′   = deepcopy(vie.f′)
    W    = deepcopy(vie.W)
    return VolterraIntegralVectorEquation(dt, N, Nₘₐₓ, n, f, g, t, f′, W)
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

function advance!(vie::VolterraIntegralVectorEquation, g′ₙ::ArrayOfPhysicalVectors, cₙ::ArrayOfPhysicalScalars)
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
    if (g′ₙ.array.rows ≠ 3) || (cₙ.array.len ≠ 3)
        msg = "The control function g′ₙ and coefficient cₙ must be of lengths 3."
        throw(ErrorException(msg))
    end
    if g′ₙ.array.cols ≠ vie.f′.array.cols
        msg = "The length of each vector in g′ₙ must be of length vector vie.f′."
        throw(ErrorException(msg))
    end
    if !isDimensionless(cₙ)
        msg = "Coefficient cₙ must be dimensionless array of scalars."
        throw(ErrorException(msg))
    end

    # Create the matrix coefficient for the right-hand side (rhs) product.
    I = PhysicalTensor(3, 3, vie.W.units)
    for i in 1:3
        I[i,i] = PhysicalScalar(1.0, vie.W.units)
    end
    W₁ = vie.W[1]
    W₁inv = inv(I + cₙ*transpose(W₁))

    # Create the temporary working arrays, which are of length 3.
    b′ = ArrayOfPhysicalVectors(3, g′ₙ.array.cols, g′ₙ.units)
    f′ = ArrayOfPhysicalVectors(3, g′ₙ.array.cols, g′ₙ.units)
    x′ = ArrayOfPhysicalVectors(3, g′ₙ.array.cols, g′ₙ.units)

    # Create the vector for the rhs product.

    # First, add in the control contribution to the rhs vector.
    for i in 1:3
        b′[i] = g′ₙ[i]
    end

    # Second, incorporate the history effects acting on this rhs vector.
    if vie.n ≤ vie.Nₘₐₓ
        # Advance the solution along a path with full history.
        for m in 1:vie.n-1
            W = vie.W[vie.n-m+1]
            f′[1] = vie.f′[3(m-1)+1]
            f′[2] = vie.f′[3(m-1)+2]
            f′[3] = vie.f′[3(m-1)+3]
            for i in 1:3
                for j in 1:3
                    b′[i] = b′[i] - cₙ[i] * W[j,i] * f′[j]
                end
            end
        end
    else  # vie.n > vie.Nₘₐₓ
        # Advance the solution along a path with truncated history.
        for m in 1:vie.Nₘₐₓ-1
            W = vie.W[vie.Nₘₐₓ-m+1]
            f′[1] = vie.f′[3(m+vie.n-vie.Nₘₐₓ-1)+1]
            f′[2] = vie.f′[3(m+vie.n-vie.Nₘₐₓ-1)+2]
            f′[3] = vie.f′[3(m+vie.n-vie.Nₘₐₓ-1)+3]
            for i in 1:3
                for j in 1:3
                    b′[i] = b′[i] - cₙ[i] * W[j,i] * f′[j]
                end
            end
        end
    end

    # Finally, solve A x′ = b′ for x′, i.e., solve the linear system.
    for i in 1:3
        for j in 1:3
            x′[i] = x′[i] + W₁inv[i,j] * b′[j]
        end
    end

    # Assign this solution to its location for nᵗʰ step in the history vector.
    for i in 1:3
        vie.f′[3(vie.n-1)+i] = x′[i]
    end

    n = get(vie.n)
    # Integrate rate expressions describing the control and response functions.
    vie.f[n+1] = (vie.f[n] + (vie.dt/8) *
        (3vie.f′[3(n-1)+1] + 2vie.f′[3(n-1)+2] + 3vie.f′[3(n-1)+3]))
    vie.g[n+1] = vie.g[n] + (vie.dt/8) * (3g′ₙ[1] + 2g′ₙ[2] + 3g′ₙ[3])

    # Update the counter.
    set!(vie.n, n+1)

    return nothing
end # advance!

# Perform an iteration of refinement on a solution at current step n. Call only
# if the control function g′ₙ or coefficient cₙ undergo iterative refinement.

function update!(vie::VolterraIntegralVectorEquation, g′ₙ::ArrayOfPhysicalVectors, cₙ::ArrayOfPhysicalScalars)

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
    dt::PhysicalScalar          # distance separating global integration nodes
    N::Integer                  # number of integration nodes in a solution path
    Nₘₐₓ::Integer               # maximum number of nodes whose history is kept
    n::MInteger                 # current node along a solution path
    # Arrays of length N+1 holding integrated variable rates at the global nodes
    f::ArrayOfPhysicalTensors   # array of integrated response function values
    g::ArrayOfPhysicalTensors   # array of integrated control function values
    t::ArrayOfPhysicalScalars   # array of times, the independent variable
    # Array of length 3N holding response function rates at over local intervalS
    f′::ArrayOfPhysicalTensors  # history array of response function rates
    # Array of Nₘₐₓ normalized weights of quadrature for a product integral
    W::ArrayOfPhysicalTensors   # array of matrices holding quadrature weights

    # constructors

    # For use when first creating this data structure.
    function VolterraIntegralTensorEquation(systemOfUnits::String, N::Integer, dt::PhysicalScalar, f₀::PhysicalTensor, g₀::PhysicalTensor, W::ArrayOfPhysicalTensors)

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
        if (𝑓₀.units ≠ 𝑔₀.units) || (𝑓₀.matrix.rows ≠ 𝑔₀.matrix.rows) || (𝑓₀.matrix.cols ≠ 𝑔₀.matrix.cols)
            msg = "Units and dimensions for initial conditions f₀ and g₀ must be equal."
            throw(ErrorException(msg))
        end
        if !isDimensionless(𝑊) || (𝑊.array.rows ≠ 3) || (𝑊.array.cols ≠ 3)
            msg = "Weights of quadrature W must be dimensionless 3x3 matrices."
            throw(ErrorException(msg))
        end

        # Create the fields for this data structure,
        Nₘₐₓ = 𝑊.array.pgs
        f = ArrayOfPhysicalTensors(N+1, 𝑓₀.matrix.rows, 𝑓₀.matrix.cols, 𝑓₀.units)
        f[1] = 𝑓₀
        g = ArrayOfPhysicalTensors(N+1, 𝑔₀.matrix.rows, 𝑔₀.matrix.cols, 𝑔₀.units)
        g[1] = 𝑔₀
        t = ArrayOfPhysicalScalars(N+1, t₀.units)
        t[1] = t₀
        for n in 1:N
            t[n+1] = n * d𝑡
        end
        f′ = ArrayOfPhysicalTensors(3N, 𝑓₀.matrix.rows, 𝑓₀.matrix.cols, 𝑓₀.units-d𝑡.units)
        n  = MInteger(1)
        new(d𝑡, N, Nₘₐₓ, n, f, g, t, f′, 𝑊)
    end

    # Used by JSON3 whenever this data structure is to be created from a file.
    function VolterraIntegralTensorEquation(dt::PhysicalScalar, N::Integer, Nₘₐₓ::Integer, n::MInteger, f::ArrayOfPhysicalTensors, g::ArrayOfPhysicalTensors, t::ArrayOfPhysicalScalars, f′::ArrayOfPhysicalTensors, W::ArrayOfPhysicalTensors)
        new(dt, N, Nₘₐₓ, n, f, g, t, f′, W)
    end
end # VolterraIntegralTensorEquation

# Methods

function Base.:(copy)(vie::VolterraIntegralTensorEquation)::VolterraIntegralTensorEquation
    dt   = copy(vie.dt)
    N    = copy(vie.N)
    Nₘₐₓ = copy(vie.Nₘₐₓ)
    n    = copy(vie.n)
    t    = copy(vie.t)
    f′   = copy(vie.f′)
    W    = copy(vie.W)
    return VolterraIntegralTensorEquation(dt, N, Nₘₐₓ, n, f, g, t, f′, W)
end

function Base.:(deepcopy)(vie::VolterraIntegralTensorEquation)::VolterraIntegralTensorEquation
    dt   = deepcopy(vie.dt)
    N    = deepcopy(vie.N)
    Nₘₐₓ = deepcopy(vie.Nₘₐₓ)
    n    = deepcopy(vie.n)
    t    = deepcopy(vie.t)
    f′   = deepcopy(vie.f′)
    W    = deepcopy(vie.W)
    return VolterraIntegralTensorEquation(dt, N, Nₘₐₓ, n, f, g, t, f′, W)
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

function advance!(vie::VolterraIntegralTensorEquation, g′ₙ::ArrayOfPhysicalTensors, cₙ::ArrayOfPhysicalScalars)
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
    if (g′ₙ.array.pgs ≠ 3) || (cₙ.array.len ≠ 3)
        msg = "The control function g′ₙ and coefficient cₙ must be of lengths 3."
        throw(ErrorException(msg))
    end
    if (g′ₙ.array.rows ≠ vie.f′.array.rows) || (g′ₙ.array.cols ≠ vie.f′.array.cols)
        msg = "The dimension of tensor g′ₙ must be that of tensor vie.f′."
        throw(ErrorException(msg))
    end
    if !isDimensionless(cₙ)
        msg = "Coefficient cₙ must be dimensionless scalar."
        throw(ErrorException(msg))
    end

    # Create the matrix coefficient for the right-hand side (rhs) product.
    I = PhysicalTensor(3, 3, vie.W.units)
    for i in 1:3
        I[i,i] = PhysicalScalar(1.0, vie.W.units)
    end
    W₁ = vie.W[1]
    W₁inv = inv(I + cₙ*transpose(W₁))

    # Create the temporary working arrays, which are of length 3.
    b′ = ArrayOfPhysicalTensors(3, g′ₙ.array.rows, g′ₙ.array.cols, g′ₙ.units)
    f′ = ArrayOfPhysicalTensors(3, g′ₙ.array.rows, g′ₙ.array.cols, g′ₙ.units)
    x′ = ArrayOfPhysicalTensors(3, g′ₙ.array.rows, g′ₙ.array.cols, g′ₙ.units)

    # Create the vector for the rhs product.

    # First, add in the control contribution to the rhs vector.
    for i in 1:3
        b′[i] = g′ₙ[i]
    end

    # Second, incorporate the history effects acting on this rhs vector.
    if vie.n ≤ vie.Nₘₐₓ
        # Advance the solution along a path with full history.
        for m in 1:vie.n-1
            W = vie.W[vie.n-m+1]
            f′[1] = vie.f′[3(m-1)+1]
            f′[2] = vie.f′[3(m-1)+2]
            f′[3] = vie.f′[3(m-1)+3]
            for i in 1:3
                for j in 1:3
                    b′[i] = b′[i] - cₙ[i] * W[j,i] * f′[j]
                end
            end
        end
    else  # vie.n > vie.Nₘₐₓ
        # Advance the solution along a path with truncated history.
        for m in 1:vie.Nₘₐₓ-1
            W = vie.W[vie.Nₘₐₓ-m+1]
            f′[1] = vie.f′[3(m+vie.n-vie.Nₘₐₓ-1)+1]
            f′[2] = vie.f′[3(m+vie.n-vie.Nₘₐₓ-1)+2]
            f′[3] = vie.f′[3(m+vie.n-vie.Nₘₐₓ-1)+3]
            for i in 1:3
                for j in 1:3
                    b′[i] = b′[i] - cₙ[i] * W[j,i] * f′[j]
                end
            end
        end
    end

    # Finally, solve A x′ = b′ for x′, i.e., solve the linear system.
    for i in 1:3
        for j in 1:3
            x′[i] = x′[i] + W₁inv[i,j] * b′[j]
        end
    end

    # Assign this solution to its location for nᵗʰ step in the history vector.
    for i in 1:3
        vie.f′[3(vie.n-1)+i] = x′[i]
    end

    n = get(vie.n)
    # Integrate rate expressions describing the control and response functions.
    vie.f[n+1] = (vie.f[n] + (vie.dt/8) *
        (3vie.f′[3(n-1)+1] + 2vie.f′[3(n-1)+2] + 3vie.f′[3(n-1)+3]))
    vie.g[n+1] = vie.g[n] + (vie.dt/8) * (3g′ₙ[1] + 2g′ₙ[2] + 3g′ₙ[3])

    # Update the counter.
    set!(vie.n, n+1)

    return nothing
end # advance!

# Perform an iteration of refinement on a solution at current step n. Call only
# if the control function g′ₙ or coefficient cₙ undergo iterative refinement.

function update!(vie::VolterraIntegralTensorEquation, g′ₙ::ArrayOfPhysicalTensors, cₙ::ArrayOfPhysicalScalars)

    set!(vie.n, get(vie.n)-1)
    advance!(vie, g′ₙ, cₙ)
    return nothing
end # update!

end # module VolterraIntegralEquations
