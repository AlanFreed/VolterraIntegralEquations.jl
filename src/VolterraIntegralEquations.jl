#=
Created on Tue 27 Jun 2023
Updated on Fri 20 Oct 2023
-------------------------------------------------------------------------------
This software, like the language it is written in, is published under the MIT
License, https://opensource.org/licenses/MIT.

Copyright (c) 2023:
Alan Freed and John Clayton

Permission is hereby granted, free of charge, to any person obtaining a copy of
this software and associated documentation files (the "Software"), to deal in
the Software without restriction, including without limitation the rights to
use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies
of the Software, and to permit persons to whom the Software is furnished to do
so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS
FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR
COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER
IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
-------------------------------------------------------------------------------
=#

#=
References:
Caputo, M. and Mainardi, F., "Linear models of dissipation in anelastic solids," Rivista del Nuoro Cimento, 1 (1971), 161-198.

Caputo, M. and Mainardi, F., "A new dissipation model based on memory mechanism," Pure and Applied Geophysics, 91 (1971), 134-147.

Cole, K.S. and Cole, R.H., "Dispersion and absorption in dielectrics I. Alternating current characteristics," Journal of Chemical Physics, 9 (1941), 342-351.

Cole, K.S. and Cole, R.H., "Dispersion and absorption in dielectrics II. Direct current characteristics," Journal of Chemical Physics, 10 (1942), 98-105.

Freed, A.D., Soft Solids: A primer to the theoretical mechanics of materials, Modeling and Simulation in Science, Engineering and Technology. Basel: Birkhäuser, 2014.

Freed, A.D. and Rajagopal, K.R. "A viscoelastic model for describing the response of biological fibers," ACTA Mechanica, 227 (2016), 3367-3380.

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
    SLS,    # Zener's 𝑆tandard 𝐿inear 𝑆olid, a.k.a. the kernel of Maxwell-Debye
    FLS,    # Caputo and Mainardi's 𝐹ractional 𝐿inear 𝑆olid
    RFS,    # Freed and Rajagopal's 𝑅egularized 𝐹L𝑆
    KWW,    # 𝐾ohlrausch's and 𝑊illiams & 𝑊atts' stretched exponential
    CCM,    # 𝐶ole and 𝐶ole's power-law 𝑀odel
    MPL,    # Williams' 𝑀odified 𝑃ower-𝐿aw model
    BOX,    # the 𝑏𝑜𝑥 model of Neuber, a.k.a. Fung's QLV kernel
    MCM,    # 𝑀axwell's 𝐶hain 𝑀odel, a.k.a. the Prony series model

    # Function used to create weights of quadrature for a given memory function.
    normalizedQuadratureWeights,

    # Solvers for Volterra integral equations of the second kind; specifically,
    #   f′(t) = g′(t) - c ∫₀ᵗ K(t-τ) f′(τ) dτ
    # where
    #   g′ is the time rate-of-change of some control function g(t)
    #   f′ is the time rate-of-change of the response function f(t)
    #   c  is a scalar coefficent, e.g., (E₀ - E∞)/E∞ in viscoelasticity
    #   K  is a memory function, i.e., the derivative of a relaxation function
    # Here f′ and g′ may be scalar, vector or tensor valued.
    # Upon solving f′, the resulting ODE can be integrated to get response f.

    # abstract type
    VolterraIntegralEquation,

    # concrete types with internal constructors
    VolterraIntegralScalarEquation,
    VolterraIntegralVectorEquation,
    VolterraIntegralTensorEquation,

    # methods

    toFile,
    fromFile,

    advance!,
    update!

#=
-------------------------------------------------------------------------------
=#

#= Memory functions are the derivatives of relaxation functions [Freed, 2014], the latter being more commonly found in the literature. All memory functions have an interface of:
    k = <memoryFunctionName>(systemOfUnits, time, parameters)\n
where `systemOfUnits` is either "SI" or "CGS", `time` is current time, and `parameters` is a tuple containing this kernel's physical parameters. Specifically, the memory functions implemented here include:
    BOX     the box model of Neuber, a.k.a. Fung's QLV kernel
    CCM     Cole and Cole's power-law Model
    FLS     Caputo and Mainardi's Fractional Linear Solid
    KWW     Kohlrausch's and Williams & Watts' stretched exponential
    MCM     Maxwell's Chain Model, a.k.a. the Prony series model
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
wherein τ denotes a characteristic time. There are two in the BOX model, and n in the MCM, arranged so that 0 < τ₁ < τ₂ < ⋯ < τₙ, with the cᵢ > 0, i = 1, 2, …, n being the coefficients of a Prony series whose sum is 1, i.e., ∑_{i=1}^n cᵢ = 1. Parameter α is the exponent in a power law, and parameter δ is a shift in time introduced to remove a weak singularity.\n
The following memory functions are weakly singular at the upper limit of integration in their Volterra integrals: BOX, CCM, FLS and KWW.
=#

"""
Memory function BOX (the box energy function of Neubert)\n
    k = (exp(-t/τ₂) - exp(-t/τ₁)) / (t ln(τ₂/τ₁)) \n
Argument `parameters` describes the tuple (τ₁, τ₂) ordered as 0 < τ₁ < τ₂.
"""
function BOX(systemOfUnits::String, time::PhysicalScalar, parameters::Tuple)::PhysicalScalar
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

    if (t.value ≥ 0.0) && (τ₁.value > 0.0) && (τ₂ > τ₁)
        return (exp(-t/τ₂) - exp(-t/τ₁)) / (t * log(τ₂/τ₁))
    else
        msg = "Argument time must be non-negative, and parameters τ₁ and τ₂ must be ordered as 0 < τ₁ < τ₂."
        throw(ErrorException(msg))
    end
end # BOX

"""
Memory function CCM (Cole-Cole power-law Model)\n
    k = (t/τ)^α (α / t) / (1 + (t/τ)^α)² \n
Argument `parameters` describes the tuple (α, τ).
"""
function CCM(systemOfUnits::String, time::PhysicalScalar, parameters::Tuple)::PhysicalScalar
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
        x = (t/τ)^get(α)
        return x * (α/t) / ((1 + x) * (1 + x))
    else
        msg = "Argument time must be non-negative, while parameters α and τ must be positive."
        throw(ErrorException(msg))
    end
end # CCM

"""
Memory function FLS (Fractional Linear Solid)\n
    k = -E_{α,0}(-(t/τ)^α) / t \n
where E_{α, β}(z) denotes the two-parameter Mittag-Leffler function, with β = 0 in this case. Argument `parameters` describes the tuple (α, τ).
"""
function FLS(systemOfUnits::String, time::PhysicalScalar, parameters::Tuple)::PhysicalScalar
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

    if ((t.value ≥ 0.0) && (τ.value > 0.0) &&
        (α.value > 0.0) && (α.value < 1.0))
        x = (t / τ)^get(α)
        return -MittagLeffler.mittleff(get(α), 0.0, -get(x)) / t
    else
        msg = "Argument time must be non-negative, τ must be positive, and α ∈ (0,1)."
        throw(ErrorException(msg))
    end
end # FLS

"""
Memory function KWW (stretched exponential of Kohlrausch, Williams and Watts)\n
    k = (t/τ)^α (α/t) exp(-(t/τ)^α) \n
Argument `parameters` describes the tuple (α, τ).
"""
function KWW(systemOfUnits::String, time::PhysicalScalar, parameters::Tuple)::PhysicalScalar
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

    if ((t.value ≥ 0.0) && (τ.value > 0.0) &&
        (α.value > 0.0) && (α.value < 1.0))
        return (t/τ)^α.value * (α/t) * exp(-(t/τ)^α.value)
    else
        msg = "Argument time must be non-negative, while parameters τ must be positive and α ∈ (0,1)."
        throw(ErrorException(msg))
    end
end # KWW

"""
Memory function MCM (Maxwell Chain Model, which is a Prony series)\n
    k = (c₁/τ₁) exp(-t/τ₁) + ⋯ + (cₙ/τₙ) exp(-t/τₙ) \n
Argument `parameters` describes a tuple (c₁, c₂, …, cₙ, τ₁, τ₂, …, τₙ) of length 2n, where ∑_{i=1}^n cᵢ = 1, and where 0 < τ₁ < τ₂ < ⋯ < τₙ.
"""
function MCM(systemOfUnits::String, time::PhysicalScalar, parameters::Tuple)::PhysicalScalar
    # Verify the inputs.
    if length(parameters) % 2 == 0
        n = length(parameters) ÷ 2
    else
        msg = "There must be an even number of parameters."
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
            if τ[i-1] < τ[i]
                msg = "Prony characteristic times must order 0 < τ₁ < ⋯ < τₙ."
                throw(ErrorException(msg))
            end
        end
        pronySeries = PhysicalScalar(PhysicalUnits("SI", 0, 0, 0, -1, 0, 0, 0))
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
        pronySeries = PhysicalScalar(PhysicalUnits("CGS", 0, 0, 0, -1, 0, 0, 0))
    else
        msg = "The assigned physical system of units is unknown."
         throw(ErrorException(msg))
    end

    if (t.value ≥ 0.0) && (τ[1].value > 0.0)
        for i in 1:n
            pronySeries = pronySeries + (c[i]/τ[i]) * exp(-t/τ[i])
        end
        return pronySeries
    else
        msg = "Argument time must be non-negative, while parameter τ₁ must be positive."
        throw(ErrorException(msg))
    end
end # MCM

"""
Memory function MPL (Modified Power-Law model of Williams')\n
    k = (α/τ) / (1 + t/τ)^(1+α) \n
Argument `parameters` describes the tuple (α, τ).
"""
function MPL(systemOfUnits::String, time::PhysicalScalar, parameters::Tuple)::PhysicalScalar
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
        return (α/τ) / (1 + t/τ)^(1+α.value)
    else
        msg = "Argument time must be non-negative, while parameters α and τ must be positive."
        throw(ErrorException(msg))
    end
end # MPL

"""
Memory function RFS (Regularized Fractional Solid)\n
    k = -E_{α,0}(-((t+δ)/τ)^α) / (E_{α,1}(-(δ/τ)^α)(t+δ)) \n
where E_{α, β}(z) denotes the two-parameter Mittag-Leffler function. Argument `parameters` describes the tuple (α, δ, τ).
"""
function RFS(systemOfUnits::String, time::PhysicalScalar, parameters::Tuple)::PhysicalScalar
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

    if ((t.value ≥ 0.0) && (α.value > 0.0) && (α.value < 1.0) &&
        (δ.value > 0.0) && (τ.value > 0.0))
        x = ((t + δ) / τ)^get(α)
        numerMLF = MittagLeffler.mittleff(get(α), 0.0, -get(x))
        y = (δ / τ)^get(α)
        denomMLF = MittagLeffler.mittleff(get(α), 1.0, -get(y))
        return -numerMLF / (denomMLF * (δ + t))
    else
        msg = "Argument time must be non-negative, parameters τ and δ must be positive, and α ∈ (0,1)."
        throw(ErrorException(msg))
    end
end # RFS

"""
Memory funtion SLS (Standard Linear Solid)\n
    k = exp(-t/τ) / τ \n
Argument `parameters` describes the tuple (τ,).
"""
function SLS(systemOfUnits::String, time::PhysicalScalar, parameters::Tuple)::PhysicalScalar
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
        return exp(-t/τ) / τ
    else
        msg = "time must be non-negative, and τ must be positive."
        throw(ErrorException(msg))
    end
end # SLS

#=
-------------------------------------------------------------------------------
=#

# The function used to create weights of quadrature. The supplied kernel `K`
# may be any of those preprogrammed above, or one of your own creation.

"""
Function\n
    W = normalizedQuadratureWeights(K, systemOfUnits, N, dTime, parameters)\n
where `K` is a memory function, `systemOfUnits` is either "SI" or "CGS", `N` is the number of nodes of integration, i.e., length of the returned array of weights, `dTime` is an uniform increment in time separating nodes from their nearest neighobors, and `parameters` is a tuple of material constants to be passed to kernel `K`.\n
The supplied memory function is to have an interface of\n
    k = K(systemOfUnits, time, parameters)\n
where `systemOfUnits` is either "SI" or "CGS", `time` is current time, and `parameters` is a tuple containing this kernel's physical parameters, i.e., its material constants.\n
The weights of quadrature returned here are normalized, e.g., the actual weights of quadrature for a viscoelastic kernel would be these normalized weights multiplied by a scalar coefficient of (E₀ - E∞)/E∞, which is to be assigned to field `c` in an object of type VolterraIntegralEquation.\n
The returned array of normalized quadrature weights is to be assigned to field `W` in an object of type VolterraIntegralEquation.
"""
function normalizedQuadratureWeights(K::Function, systemOfUnits::String, N::Integer, dTime::PhysicalScalar, parameters::Tuple)::ArrayOfPhysicalTensors

    # Ensure the system of units is consistent.
    if (systemOfUnits == "SI") || (systemOfUnits == "si")
        dt = toSI(dTime)
    elseif (systemOfUnits == "CGS") || (systemOfUnits == "cgs")
        dt = toCGS(dTime)
    else
        msg = "The assigned physical system of units is unknown."
        throw(ErrorException(msg))
    end

    # Basic arrays needed to create the weights of quadrature.

    # Inverse of the 3x3 midpoint alternant matrix X, i.e., Xinv.
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

    # The three Gauss quadrature matrices.
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

    # The three Gauss quadrature vectors.
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

    k = K(systemOfUnits, dt, parameters)  # call made to establish units.
    μ₁ = PhysicalTensor(3, 3, dt.units+k.units)
    for j in 1:3
        coef = (j - 0.5)*dt / 6
        for i in 1:3
            sum = PhysicalScalar(k.units)
            for s in 1:3
                t = (j - 0.5) * (1 - x[s]) * dt / 6
                sum = sum + w[s] * m[i,j,s] * K(systemOfUnits, t, parameters)
            end
            μ₁[i,j] = coef * sum
        end
    end
    quadMatrix = PhysicalTensor(3, 3, μ₁.units)
    for i in 1:3
        for j in 1:3
            sum = PhysicalScalar(μ₁.units)
            for k in 1:3
                sum = sum + Xinv[i,k] * μ₁[k,j]
            end
            quadMatrix[i,j] = sum
        end
    end
    quadWgts = ArrayOfPhysicalTensors(N, 3, 3, μ₁.units)
    quadWgts[1] = quadMatrix

    # Create the remaining moment matrices.
    coef = dt / 2
    for n in 2:N
        μₙ = PhysicalTensor(3, 3, μ₁.units)
        for i in 1:3
            for j in 1:3
                sum = PhysicalScalar(k.units)
                for s in 1:3
                    t = (n - (5 - j)/3 - x[s]/2) * dt
                    sum = sum + w[s] * v[i,s] * K(systemOfUnits, t, parameters)
                end
                μₙ[i,j] = coef * sum
            end
        end
        quadMatrix = PhysicalTensor(3, 3, μₙ.units)
        for i in 1:3
            for j in 1:3
                sum = PhysicalScalar(μₙ.units)
                for k in 1:3
                    sum = sum + Xinv[i,k] * μₙ[k,j]
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
    n::MInteger                 # current node along a solution path
    N::Integer                  # number of integration nodes in solution path
    Nₘₐₓ::Integer               # maximum number of nodes whose history is kept
    dt::PhysicalScalar          # distance separating global integration nodes
    # Arrays of length N+1 holding integrated variable rates at the global nodes
    f::ArrayOfPhysicalScalars   # array of integrated response function values
    g::ArrayOfPhysicalScalars   # array of integrated control function values
    t::ArrayOfPhysicalScalars   # array of times, the independent variable
    # Array of length 3N holding response function rates at the local nodes
    f′::ArrayOfPhysicalScalars  # array of response function rates
    # Coefficient that scales the normalized weights of quadrature
    c::PhysicalScalar           # e.g., c = (E₀ - E∞)/E∞ in viscoelaticity
    # Array of Nₘₐₓ normalized weights of quadrature for a product integral
    W::ArrayOfPhysicalTensors   # array of matrices holding quadrature weights

    # constructors

    # For use when first creating this data structure.
    function VolterraIntegralScalarEquation(systemOfUnits::String, N::Integer, dt::PhysicalScalar, f₀::PhysicalScalar, g₀::PhysicalScalar, c::PhysicalScalar, W::ArrayOfPhysicalTensors)

        # Ensure that a consistent system of physical units is used.
        if (systemOfUnits == "SI") || (systemOfUnits == "si")
            d𝑡 = toSI(dt)
            𝑓₀ = toSI(f₀)
            𝑔₀ = toSI(g₀)
            𝑐 = toSI(c)
            𝑊 = toSI(W)
            t₀ = PhysicalScalar(d𝑡.units)
        elseif (systemOfUnits == "CGS") || (systemOfUnits == "cgs")
            d𝑡 = toCGS(dt)
            𝑓₀ = toCGS(f₀)
            𝑔₀ = toCGS(g₀)
            𝑐 = toCGS(c)
            𝑊 = toCGS(W)
            t₀ = PhysicalScalar(d𝑡.units)
        else
            msg = "The assigned physical system of units is unknown."
            throw(ErrorException(msg))
        end

        # Verify the remaining inputs.
        if N < 1
            msg = "The number of nodes N must be positive valued."
            throw(ErrorException(msg))
        end
        if 𝑓₀.units ≠ 𝑔₀.units
            msg = "Physical units for initial conditions f₀ and g₀ must be equal."
            throw(ErrorException(msg))
        end
        if !isDimensionless(c)
            msg = "Coefficient c scaling the weights of quadrature must be dimensionless."
            throw(ErrorException(msg))
        end
        if !isDimensionless(𝑊) || (𝑊.array.rows ≠ 3) || (𝑊.array.cols ≠ 3)
            msg = "Weights of quadrature W must be dimensionless 3x3 matrices."
            throw(ErrorException(msg))
        end

        # Create the fields for this data structure,
        n = MInteger(1)
        Nₘₐₓ = 𝑊.array.pgs
        f = ArrayOfPhysicalScalars(N+1, 𝑓₀.units)
        f[1] = 𝑓₀
        g = ArrayOfPhysicalScalars(N+1, 𝑔₀.units)
        g[1] = 𝑔₀
        t = ArrayOfPhysicalScalars(N+1, t₀.units)
        t[1] = t₀
        for n in 2:N+1
            t[n] = t[n-1] + d𝑡
        end
        f′ = ArrayOfPhysicalScalars(3N, 𝑓₀.units-d𝑡.units)
        new(n, N, Nₘₐₓ, d𝑡, f, g, t, f′, 𝑐, 𝑊)
    end

    # Used by JSON3 whenever this data structure is to be created from a file.
    function VolterraIntegralScalarEquation(n::MInteger, N::Integer, Nₘₐₓ::Integer, dt::PhysicalScalar, f::ArrayOfPhysicalScalars, g::ArrayOfPhysicalScalars, t::ArrayOfPhysicalScalars, f′::ArrayOfPhysicalScalars, c::PhysicalScalar, W::ArrayOfPhysicalTensors)
        new(n, N, Nₘₐₓ, dt, f, g, t, f′, c, W)
    end
end # VolterraIntegralScalarEquation

# Vector-valued Volterra integral equations of the second kind.

struct VolterraIntegralVectorEquation <: VolterraIntegralEquation
    # Dimensioning fields
    n::MInteger                 # current node along a solution path
    N::Integer                  # number of integration nodes in solution path
    Nₘₐₓ::Integer               # maximum number of nodes whose history is kept
    dt::PhysicalScalar          # distance separating global integration nodes
    # Arrays of length N+1 holding integrated variable rates at the global nodes
    f::ArrayOfPhysicalVectors   # array of integrated response function values
    g::ArrayOfPhysicalVectors   # array of integrated control function values
    t::ArrayOfPhysicalScalars   # array of times, the independent variable
    # Array of length 3N holding response function rates at the local nodes
    f′::ArrayOfPhysicalVectors  # array of response function rates
    # Coefficient that scales the normalized weights of quadrature
    c::PhysicalScalar           # e.g., c = (E₀ - E∞)/E∞ in viscoelaticity
    # Array of Nₘₐₓ normalized weights of quadrature for a product integral
    w::ArrayOfPhysicalTensors   # array of matrices holding quadrature weights

    # constructors

    # For use when first creating this data structure.
    function VolterraIntegralVectorEquation(systemOfUnits::String, N::Integer, dt::PhysicalScalar, f₀::PhysicalVector, g₀::PhysicalVector, c::PhysicalScalar, W::ArrayOfPhysicalTensors)

        # Ensure that a consistent system of physical units is used.
        if (systemOfUnits == "SI") || (systemOfUnits == "si")
            d𝑡 = toSI(dt)
            𝑓₀ = toSI(f₀)
            𝑔₀ = toSI(g₀)
            𝑐 = toSI(c)
            𝑊 = toSI(W)
            t₀ = PhysicalScalar(d𝑡.units)
        elseif (systemOfUnits == "CGS") || (systemOfUnits == "cgs")
            d𝑡 = toCGS(dt)
            𝑓₀ = toCGS(f₀)
            𝑔₀ = toCGS(g₀)
            𝑐 = toCGS(c)
            𝑊 = toCGS(W)
            t₀ = PhysicalScalar(d𝑡.units)
        else
            msg = "The assigned physical system of units is unknown."
            throw(ErrorException(msg))
        end

        # Verify the remaining inputs.
        if N < 1
            msg = "The number of nodes N must be positive valued."
            throw(ErrorException(msg))
        end
        if (𝑓₀.units ≠ 𝑔₀.units) || (𝑓₀.vector.len ≠ 𝑔₀.vector.len)
            msg = "Units and dimensions for initial conditions f₀ and g₀ must be equal."
            throw(ErrorException(msg))
        end
        if !isDimensionless(c)
            msg = "Coefficient c scaling the weights of quadrature must be dimensionless."
            throw(ErrorException(msg))
        end
        if !isDimensionless(𝑊) || (𝑊.array.rows ≠ 3) || (𝑊.array.cols ≠ 3)
            msg = "Weights of quadrature W must be dimensionless 3x3 matrices."
            throw(ErrorException(msg))
        end

        # Create the fields for this data structure,
        n = MInteger(1)
        Nₘₐₓ = 𝑊.array.pgs
        f = ArrayOfPhysicalVectors(N+1, 𝑓₀.vector.len, 𝑓₀.units)
        f[1] = 𝑓₀
        g = ArrayOfPhysicalVectors(N+1, 𝑔₀.vector.len, 𝑔₀.units)
        g[1] = 𝑔₀
        t = ArrayOfPhysicalScalars(N+1, t₀.units)
        t[1] = t₀
        for n in 2:N+1
            t[n] = t[n-1] + d𝑡
        end
        f′ = ArrayOfPhysicalVectors(3N, 𝑓₀.vector.len, 𝑓₀.units-d𝑡.units)
        new(n, N, Nₘₐₓ, d𝑡, f, g, t, f′, 𝑐, 𝑊)
    end

    # Used by JSON3 whenever this data structure is to be created from a file.
    function VolterraIntegralVectorEquation(n::MInteger, N::Integer, Nₘₐₓ::Integer, dt::PhysicalScalar, f::ArrayOfPhysicalVectors, g::ArrayOfPhysicalVectors, t::ArrayOfPhysicalScalars, f′::ArrayOfPhysicalVectors, c::PhysicalScalar, W::ArrayOfPhysicalTensors)
        new(n, N, Nₘₐₓ, dt, f₀, g₀, t₀, f, g, t, f′, c, W)
    end
end # VolterraIntegralVectorEquation

# Tensor-valued Volterra integral equations of the second kind.

struct VolterraIntegralTensorEquation <: VolterraIntegralEquation
    # Dimensioning fields
    n::MInteger                 # current node along a solution path
    N::Integer                  # number of integration nodes in solution path
    Nₘₐₓ::Integer               # maximum number of nodes whose history is kept
    dt::PhysicalScalar          # distance separating global integration nodes
    # Arrays of length N+1 holding integrated variable rates at the global nodes
    f::ArrayOfPhysicalTensors   # array of integrated response function values
    g::ArrayOfPhysicalTensors   # array of integrated control function values
    t::ArrayOfPhysicalScalars   # array of times, the independent variable
    # Array of length 3N holding response function rates at the local nodes
    f′::ArrayOfPhysicalTensors  # array of response function rates
    # Coefficient that scales the normalized weights of quadrature
    c::PhysicalScalar           # e.g., c = (E₀ - E∞)/E∞ in viscoelaticity
    # Array of Nₘₐₓ normalized weights of quadrature for a product integral
    w::ArrayOfPhysicalTensors   # array of matrices holding quadrature weights

    # constructors

    # For use when first creating this data structure.
    function VolterraIntegralTensorEquation(systemOfUnits::String, N::Integer, dt::PhysicalScalar, f₀::PhysicalTensor, g₀::PhysicalTensor, c::PhysicalScalar, W::ArrayOfPhysicalTensors)

        # Ensure that a consistent system of physical units is used.
        if (systemOfUnits == "SI") || (systemOfUnits == "si")
            d𝑡 = toSI(dt)
            𝑓₀ = toSI(f₀)
            𝑔₀ = toSI(g₀)
            𝑐 = toSI(c)
            𝑊 = toSI(W)
            t₀ = PhysicalScalar(d𝑡.units)
        elseif (systemOfUnits == "CGS") || (systemOfUnits == "cgs")
            d𝑡 = toCGS(dt)
            𝑓₀ = toCGS(f₀)
            𝑔₀ = toCGS(g₀)
            𝑐 = toCGS(c)
            𝑊 = toCGS(W)
            t₀ = PhysicalScalar(d𝑡.units)
        else
            msg = "The assigned physical system of units is unknown."
            throw(ErrorException(msg))
        end

        # Verify the remaining inputs.
        if N < 1
            msg = "The number of nodes N must be positive valued."
            throw(ErrorException(msg))
        end
        if ((𝑓₀.units ≠ 𝑔₀.units) ||
            (𝑓₀.matrix.rows ≠ 𝑔₀.matrix.rows) ||
            (𝑓₀.matrix.cols ≠ 𝑔₀.matrix.cols))
            msg = "Units and dimensions for initial conditions f₀ and g₀ must be equal."
            throw(ErrorException(msg))
        end
        if !isDimensionless(c)
            msg = "Coefficient c scaling the weights of quadrature must be dimensionless."
            throw(ErrorException(msg))
        end
        if !isDimensionless(𝑊) || (𝑊.array.rows ≠ 3) || (𝑊.array.cols ≠ 3)
            msg = "Weights of quadrature W must be dimensionless 3x3 matrices."
            throw(ErrorException(msg))
        end

        # Create the fields for this data structure,
        n = MInteger(1)
        Nₘₐₓ = 𝑊.array.pgs
        f = ArrayOfPhysicalTensors(N+1, 𝑓₀.vector.len, 𝑓₀.units)
        f[1] = 𝑓₀
        g = ArrayOfPhysicalTensors(N+1, 𝑔₀.vector.len, 𝑔₀.units)
        g[1] = 𝑔₀
        t = ArrayOfPhysicalScalars(N+1, t₀.units)
        t[1] = t₀
        for n in 2:N+1
            t[n] = t[n-1] + d𝑡
        end
        f′ = ArrayOfPhysicalTensors(3N, 𝑓₀.matrix.rows, 𝑓₀.matrix.cols, 𝑓₀.units-d𝑡.units)
        new(n, N, Nₘₐₓ, d𝑡, f, g, t, f′, 𝑐, 𝑊)
    end

    # Used by JSON3 whenever this data structure is to be created from a file.
    function VolterraIntegralTensorEquation(n::MInteger, N::Integer, Nₘₐₓ::Integer, dt::PhysicalScalar, f::ArrayOfPhysicalTensors, g::ArrayOfPhysicalTensors, t::ArrayOfPhysicalScalars, f′::ArrayOfPhysicalTensors, c::PhysicalScalar, W::ArrayOfPhysicalTensors)
        new(n, N, Nₘₐₓ, dt, f, g, t, f′, c, W)
    end
end # VolterraIntegralTensorEquation

#=
-------------------------------------------------------------------------------
=#

# Reading from and writing to a JSON file.

StructTypes.StructType(::Type{VolterraIntegralScalarEquation}) = StructTypes.Struct()
StructTypes.StructType(::Type{VolterraIntegralVectorEquation}) = StructTypes.Struct()
StructTypes.StructType(::Type{VolterraIntegralTensorEquation}) = StructTypes.Struct()

function toFile(y::VolterraIntegralScalarEquation, json_stream::IOStream)
    if isopen(json_stream)
        JSON3.write(json_stream, y)
        write(json_stream, '\n')
    else
        msg = "The supplied JSON stream is not open."
        throw(ErrorException(msg))
    end
    flush(json_stream)
    return nothing
end

function toFile(y::VolterraIntegralVectorEquation, json_stream::IOStream)
    if isopen(json_stream)
        JSON3.write(json_stream, y)
        write(json_stream, '\n')
    else
        msg = "The supplied JSON stream is not open."
        throw(ErrorException(msg))
    end
    flush(json_stream)
    return nothing
end

function toFile(y::VolterraIntegralTensorEquation, json_stream::IOStream)
    if isopen(json_stream)
        JSON3.write(json_stream, y)
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
        ps = JSON3.read(readline(json_stream), VolterraIntegralScalarEquation)
    else
        msg = "The supplied JSON stream is not open."
        throw(ErrorException(msg))
    end
    return ps
end

function fromFile(::Type{VolterraIntegralVectorEquation}, json_stream::IOStream)::VolterraIntegralVectorEquation
    if isopen(json_stream)
        pv = JSON3.read(readline(json_stream), VolterraIntegralVectorEquation)
    else
        msg = "The supplied JSON stream is not open."
        throw(ErrorException(msg))
    end
    return pv
end

function fromFile(::Type{VolterraIntegralTensorEquation}, json_stream::IOStream)::VolterraIntegralTensorEquation
    if isopen(json_stream)
        pt = JSON3.read(readline(json_stream), VolterraIntegralTensorEquation)
    else
        msg = "The supplied JSON stream is not open."
        throw(ErrorException(msg))
    end
    return pt
end

#=
-------------------------------------------------------------------------------
=#

# methods

function advance!(VIE::VolterraIntegralEquation, g′ₙ::Tuple)
    # verify input
    (g′ₙ₁, g′ₙ₂, g′ₙ₃) = g′ₙ
    f′₁ = VIE.f′[1]
    if ((g′ₙ₁.units ≠ f′₁.units) ||
        (g′ₙ₂.units ≠ f′₁.units) ||
        (g′ₙ₃.units ≠ f′₁.units))
        msg = string("The units of g′ₙ are ", toString(g′ₙ₃.units), " and should be ", toString(f′₁.units), ".")
        throw(ErrorException(msg))
    end
    if isa(VIE, VolterraIntegralScalarEquation)
        if (!isa(g′ₙ₁, PhysicalScalar) ||
            !isa(g′ₙ₂, PhysicalScalar) ||
            !isa(g′ₙ₃, PhysicalScalar))
            msg = "Control rates g′ₙ must be a tuple of 3 PhysicalScalars."
            throw(ErrorException(msg))
        end
    elseif isa(VIE, VolterraIntegralVectorEquation)
        if (!isa(g′ₙ₁, PhysicalVector) ||
            !isa(g′ₙ₂, PhysicalVector) || (g′ₙ₂.vector.len ≠ g′ₙ₁.vector.len) ||
            !isa(g′ₙ₃, PhysicalVector) || (g′ₙ₃.vector.len ≠ g′ₙ₁.vector.len))
            msg = "Control rates g′ₙ must be a tuple of 3 PhysicalVectors."
            throw(ErrorException(msg))
        end
        if g′ₙ₁.vector.len ≠ f′₁.array.cols
            msg = "The dimensions of vectors g′ₙ are inadmissible."
            throw(ErrorException(msg))
        end
    elseif isa(VIE, VolterraIntegralTensorEquation)
        if (!isa(g′ₙ₁, PhysicalTensor) ||
            !isa(g′ₙ₂, PhysicalTensor) ||
                (g′ₙ₂.matrix.rows ≠ g′ₙ₁.matrix.rows) ||
                (g′ₙ₂.matrix.cols ≠ g′ₙ₁.matrix.cols) ||
            !isa(g′ₙ₃, PhysicalTensor) ||
                (g′ₙ₃.matrix.rows ≠ g′ₙ₁.matrix.rows) ||
                (g′ₙ₃.matrix.cols ≠ g′ₙ₁.matrix.cols))
            msg = "Control rates g′ₙ must be a tuple of 3 PhysicalTensors."
            throw(ErrorException(msg))
        end
        if ((g′ₙ₁.matrix.rows ≠ f′₁.array.rows) ||
            (g′ₙ₁.matrix.cols ≠ f′₁.array.cols))
            msg = "The dimensions of tensors g′ₙ are inadmissible."
            throw(ErrorException(msg))
        end
    else
        msg = "The supplied Volterra integral equation VIE is of unknown type."
        throw(ErrorException(msg))
    end

    # Create the temporary working arrays, which are of length 3.
    if isa(VIE, VolterraIntegralScalarEquation)
        zero = PhysicalScalar(g′ₙ₁.units)
        f′ = ArrayOfPhysicalScalars(3, g′ₙ₁.units)
        x′ = ArrayOfPhysicalScalars(3, g′ₙ₁.units)
        y′ = ArrayOfPhysicalScalars(3, g′ₙ₁.units)
    elseif isa(VIE, VolterraIntegralVectorEquation)
        zero = PhysicalVector(g′ₙ₁.vector.len, g′ₙ₁.units)
        f′ = ArrayOfPhysicalVectors(3, g′ₙ₁.vector.len, g′ₙ₁.units)
        x′ = ArrayOfPhysicalVectors(3, g′ₙ₁.vector.len, g′ₙ₁.units)
        y′ = ArrayOfPhysicalVectors(3, g′ₙ₁.vector.len, g′ₙ₁.units)
    elseif isa(VIE, VolterraIntegralTensorEquation)
        zero = PhysicalTensor(g′ₙ₁.matrix.rows, g′ₙ₁.matrix.cols, g′ₙ₁.units)
        f′ = ArrayOfPhysicalTensors(3, g′ₙ₁.matrix.rows, g′ₙ₁.matrix.cols, g′ₙ₁.units)
        x′ = ArrayOfPhysicalTensors(3, g′ₙ₁.matrix.rows, g′ₙ₁.matrix.cols, g′ₙ₁.units)
        y′ = ArrayOfPhysicalTensors(3, g′ₙ₁.matrix.rows, g′ₙ₁.matrix.cols, g′ₙ₁.units)
    else
        msg = "The supplied Volterra integral equation VIE does not exist."
        throw(ErrorException(msg))
    end

    # Create the matrix coefficient for the right-hand side (rhs) product.
    I = PhysicalTensor(3, 3, VIE.W.units)
    one = PhysicalScalar(1.0, VIE.W.units)
    for i in 1:3
        I[i,i] = one
    end
    W₁ = VIE.W[1]
    W₁inv = inv(I + transpose(W₁))

    # Create the vector for the rhs product.
    # First, add in the control contribution to the rhs vector.
    for i in 1:3
        x′[i] = g′ₙ[i]
    end

    # Second, incorporate the history effects acting on this rhs vector.
    n = get(VIE.n)
    if n ≤ VIE.Nₘₐₓ
        # Advance the solution along a path with full history.
        for m in 1:n-1
            W = VIE.W[n-m+1]
            f′[1] = VIE.f′[3(m-1)+1]
            f′[2] = VIE.f′[3(m-1)+2]
            f′[3] = VIE.f′[3(m-1)+3]
            for i in 1:3
                for j in 1:3
                    x′[i] = x′[i] - W[j,i] * f′[j]
                end
            end
        end
    else # VIE.n > VIE.Nₘₐₓ
        # Advance the solution along a path with a truncated history.
        for m in 1:VIE.Nₘₐₓ-1
            W = VIE.W[VIE.Nₘₐₓ-m+1]
            f′[1] = VIE.f′[3(m+n-VIE.Nₘₐₓ-1)+1]
            f′[2] = VIE.f′[3(m+n-VIE.Nₘₐₓ-1)+2]
            f′[3] = VIE.f′[3(m+n-VIE.Nₘₐₓ-1)+3]
            for i in 1:3
                for j in 1:3
                    x′[i] = x′[i] - W[j,i] * f′[j]
                end
            end
        end
    end

    # Finally, compute matrix-vector product, i.e., solve the linear equations.
    for i in 1:3
        y′[i] = zero
        for j in 1:3
            y′[i] = y′[i] + W₁inv[i,j] * x′[j]
        end
    end

    # Assign this solution to its location for nᵗʰ step in the history vector.
    for i in 1:3
        VIE.f′[3(n-1)+i] = y′[i]
    end

    # Integrate rate expressions describing the control and response functions.
    VIE.f[n+1] = (VIE.f[n] + (VIE.dt/8) *
        (3VIE.f′[3(n-1)+1] + 2VIE.f′[3(n-1)+2] + 3VIE.f′[3(n-1)+3]))
    VIE.g[n+1] = VIE.g[n] + (VIE.dt/8) * (3g′ₙ[1] + 2g′ₙ[2] + 3g′ₙ[3])

    # Update the counter.
    if VIE.n < VIE.N
        n = n + 1
        set!(VIE.n, n)
    else
        println("The Volterra integral solution has reached its endpoint.")
    end

    return nothing
end # advance!

function update!(VIE::VolterraIntegralEquation, g′ₙ::Tuple)
    # Call only if control functions g′ₙ require iterative refinement.
    n = get(VIE.n)
    n = n - 1
    set!(VIE.n, n)
    advance!(VIE, g′ₙ)
    return nothing
end # update!

end # module VolterraIntegralEquations
