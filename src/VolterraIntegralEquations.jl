#=
Created on Tue 27 Jun 2023
Updated on Tue 22 Oct 2024
=#

"""
# VolterraIntegralEquations

Package `VolterraIntegralEquations.jl` provides a numerical method [1,15] for solving Volterra integral equations of the second kind, including those with singular kernels at the upper limit of integration, viz.,
```julia
f′(t) = g′(t) - c(t) ∫₀ᵗ K(t-τ) f′(τ) dτ
```
where

    g′ is the time rate-of-change of some control function g(t).
    f′ is the time rate-of-change of its response function f(t).
    c  is a scalar function, e.g., (E₀ - E∞)/E∞ in linear viscoelasticity.
    K  is a memory function, i.e., the derivative of a creep function.
    
Here f′ and g′ may be scalar, vector or tensor valued time derivatives. Upon solving f′, the resulting ODE is then integrated using a stable, backward, difference formula (BDF) to get the response function f.

To be able to use this package one must first import the following package into your project
```
using Pkg
Pkg.add(url = "https://github.com/AlanFreed/PhysicalFields.jl")
```
In the code fragments that follow, `PF` is used as an alias for `PhysicalFields`.

See the "README.md" file in the root directory of `VolterraIntegralEquations.jl` for a more detailed description of what this module does, and how it does it.

## Memory Functions

Memory function `K` is the derivative of a creep function [6], the latter being more commonly found in the literature. All memory functions are to be implemented via the interface:
```julia
(name, k, τ) = <memoryFunctionName>(system_of_units, time, parameters)
```
which returns a tuple with entries:

    name is what the kernel is called, typically an acronym, e.g., "FLS".
    k    is a value for memory function K, evaluated at the specified time.
    τ    is the controlling characteristic time--the smallest characteristic
         time whenever multiple ones occur in a kernel.

The first argument to a kernel call, i.e., `system_of_units`, is a `String` that can be either "SI" or "CGS", at present. (Other systems of units, viz., instances of type `PF.PhysicalUnits`, may be added at a later time.) Its second argument, viz., `time`, is an instance of type `PF.PhysicalScalar` that contains the current time, while its final argument, `parameters`, is a tuple containing this kernel's physical parameters.

Specifically, the following memory functions have been implemented:

    BOX  The box kernel of Neubert [11], a.k.a. Fung's [8] QLV kernel.
    CCM  Cole and Cole's [4,5] power-law Model.
    FLS  Caputo and Mainardi's [2,3] Fractional Linear Solid.
    KWW  Kohlrausch's [9] and Williams & Watts' [13] stretched exponential.
    MCM  Maxwell's [10] Chain Model, a.k.a. a Prony series.
    MPL  Williams' [14] Modified Power-Law.
    RFS  Freed and Rajagopal's [7] Regularized FLS.
    SLS  Zener's [16] Standard Linear Solid, a.k.a. the Maxwell-Debye kernel.
    
whose material parameters are supplied via the following tuples:

    BOX  parameters = (τ₁, τ₂)
    CCM  parameters = (α, τ)
    FLS  parameters = (α, τ)
    KWW  parameters = (α, τ)
    MCM  parameters = (c₁, c₂, …, cₙ, τ₁, τ₂ …, τₙ)
    MPL  parameters = (α, τ)
    RFS  parameters = (α, δ, τ)
    SLS  parameters = (τ,)
    
wherein `τ` denotes a characteristic time for creep, which has physical units of time. There are two characteristic times in the BOX model, and n characteristic times in the MCM, arranged so that
```julia
0 < τ₁ < τ₂ < ⋯ < τₙ with cᵢ > 0, i = 1, 2, …, n, and ∑ᵢ₌₁ⁿ cᵢ = 1
```
with the `cᵢ` being coefficients in the Prony series. Parameter `α` is an exponent in a power law, and parameter `δ` shifts time to forego a weak singularity.

### Weights of Quadrature

In the solution methodology of Young [15], integrals of these kernels become weights of quadrature in a numerical method that approximates Volterra integral equations of the second kind. (They become Volterra integral equations of the first kind whenever `g'(t) = 0`.) These weights are created via the function
```julia
W = normalizedQuadratureWeights(system_of_units, N, dtime, kernel, parameters)
```
where

1. `system_of_units` is a String that specifies the physical system of units being adopted, which at present may be either "SI" or "CGS".
2. `N` is an integer representing the number of time steps required to span the path to be traversed by a solution.
3. `dtime` is a `PF.PhysicalScalar` (i.e., a number with physical units, viz., `PF.PhysicalUnits`) representing a step size in time, with time being uniformly segmented over the span of its solution by this increment. The independent variable is denoted here as time, because it is typically time, but this need not be the case. 
4. `kernel` is that function which is to be called internally to quantify this memory function.
5. `parameters` is a Tuple of the kernel's parameters, as specified above.
6. `W` is an `PF.ArrayOfPhysicalScalars` that contains the associated weights of quadrature.

The weights of quadrature `W` returned here are normalized, e.g., the actual weights of quadrature for a linear viscoelastic kernel would be these normalized weights multiplied by a scalar coefficient of `(E₀ - E∞)/E∞`. These normalized quadrature weights are to be assigned to field `W` in an object implementing abstract type `VolterraIntegralEquation`, while coefficients that scale these weights, like `(E₀ - E∞)/E∞`, are to be assigned to the field `c` in this object.

## Types

The abstract type for solvers of Volterra integral equations is
```julia
abstract type VolterraIntegralEquation end
```
There are three solvers exported here that solve Volterra integral equations; they are:

### VolterraIntegralScalarEquation

```julia
struct VolterraIntegralScalarEquation <: VolterraIntegralEquation
    # Dimensioning fields
    dt::PF.PhysicalScalar          # distance separating neighboring solution nodes
    N::Int                         # number of integration nodes along a solution path
    n::PF.MInteger                 # current node along a solution path
    # Arrays of length N+1 holding control and response fields, and their rates
    f::PF.ArrayOfPhysicalScalars   # array of integrated response function values
    f′::PF.ArrayOfPhysicalScalars  # array of response function rates
    g::PF.ArrayOfPhysicalScalars   # array of integrated control function values
    g′::PF.ArrayOfPhysicalScalars  # array of control function rates
    t::PF.ArrayOfPhysicalScalars   # array of times, the independent variable
    # Array of N normalized weights of quadrature for a product integral
    W::PF.ArrayOfPhysicalScalars   # array holding the weights of quadrature
end
```
with constructor
```julia
vie = VolterraIntegralScalarEquation(system_of_units, N, dt, f₀, g₀, W)
```
where

1. `vie` is a new instance of type `VolterraIntegralScalarEquation`.
2. `system_of_units` is a string, viz., either "SI" or "CGS".
3. `N`  is the number of nodes of quadrature to be imposed.
4. `dt` is a step size separating all neighboring nodes--a `PF.PhyscalScalar`.
5. `f₀` is an initial condition for the response function--a `PF.PhysicalScalar`.
6. `g₀` is an initial condition for the control function--a `PF.PhysicalScalar`.
7. `W`  holds the weights of quadrature--an `PF.ArrayOfPhysicalScalars`.

### VolterraIntegralVectorEquation

```julia
struct VolterraIntegralVectorEquation <: VolterraIntegralEquation
    # Dimensioning fields
    dt::PF.PhysicalScalar          # distance separating neighboring solution nodes
    N::Int                         # number of integration nodes along a solution path
    n::PF.MInteger                 # current node along a solution path
    # Arrays of length N+1 holding control and response fields, and their rates
    f::PF.ArrayOfPhysicalVectors   # array of integrated response function values
    f′::PF.ArrayOfPhysicalVectors  # array of response function rates
    g::PF.ArrayOfPhysicalVectors   # array of integrated control function values
    g′::PF.ArrayOfPhysicalVectors  # array of control function rates
    t::PF.ArrayOfPhysicalScalars   # array of times, the independent variable
    # Array of N normalized weights of quadrature for a product integral
    W::PF.ArrayOfPhysicalScalars   # array holding the weights of quadrature
end
```
with constructor
```julia
vie = VolterraIntegralVectorEquation(system_of_units, N, dt, f₀, g₀, W)
```
where

1. `vie` is a new instance of type `VolterraIntegralVectorEquation`.
2. `system_of_units` is a string, viz., either "SI" or "CGS".
3. `N`  is the number of nodes of quadrature to be imposed.
4. `dt` is a step size separating all neighboring nodes--a `PF.PhyscalScalar`.
5. `f₀` is an initial condition for the response function--a `PF.PhysicalVector`.
6. `g₀` is an initial condition for the control function--a `PF.PhysicalVector`.
7. `W`  holds the weights of quadrature--an `PF.ArrayOfPhysicalScalars`.

### VolterraIntegralTensorEquation

```julia
struct VolterraIntegralTensorEquation <: VolterraIntegralEquation
    # Dimensioning fields
    dt::PF.PhysicalScalar          # distance separating neighboring solution nodes
    N::Int                         # number of integration nodes along a solution path
    n::PF.MInteger                 # current node along a solution path
    # Arrays of length N+1 holding control and response fields, and their rates
    f::PF.ArrayOfPhysicalTensors   # array of integrated response function values
    f′::PF.ArrayOfPhysicalTensors  # array of response function rates
    g::PF.ArrayOfPhysicalTensors   # array of integrated control function values
    g′::PF.ArrayOfPhysicalTensors  # array of control function rates
    t::PF.ArrayOfPhysicalScalars   # array of times, the independent variable
    # Array of N normalized weights of quadrature for a product integral
    W::PF.ArrayOfPhysicalScalars   # array holding the weights of quadrature
end
```
with constructor
```julia
vie = VolterraIntegralTensorEquation(system_of_units, N, dt, f₀, g₀, W)
```
where

1. `vie` is a new instance of type `VolterraIntegralTensorEquation`.
2. `system_of_units` is a string, viz., either "SI" or "CGS".
3. `N`  is the number of nodes of quadrature to be imposed.
4. `dt` is a step size separating all neighboring nodes--a `PF.PhyscalScalar`.
5. `f₀` is an initial condition for the response function--a `PF.PhysicalTensor`.
6. `g₀` is an initial condition for the control function--a `PF.PhysicalTensor`.
7. `W`  holds the weights of quadrature--an `PF.ArrayOfPhysicalScalars`.

## Methods

In the following methods, variable `vie` denotes an instance of either `VolterraIntegralScalarEquation` or `VolterraIntegralVectorEquation` or `VolterraIntegralTensorEquation`, as determined by its context.

For creating a copy
```julia
cc = copy(vie)
```
returns a copy of the Volterra integral equation (vie).

For storing the object to file
```julia
toFile(vie, json_stream)
```
and then reading it back in from file
```julia
vie = fromFile(VolterraIntegralScalarEquation, json_stream)  # for scalars
vie = fromFile(VolterraIntegralVectorEquation, json_stream)  # for vectors
vie = fromFile(VolterraIntegralTensorEquation, json_stream)  # for tensors
```
wherein `json_stream` is an instance of `IOStream`. (See package `PhysicalFields` for functions to create such streams.)

To solve such a Volterra integral equation, call
```julia
advance!(vie, g′ₙ, cₙ)
```
which advances a solution step-by-step, i.e., from step n-1 to step n, while
```julia
update!(vie, g′ₙ, cₙ)
```
reevaluates a solution at step n, which only needs to be called whenever either the rate of a control function `g′ₙ` or the coefficient `cₙ` to the integral vary, i.e., they have an implicit dependence upon the solution. Here `g′ₙ` will be either a `PF.PhysicalScalar`, a `PF.PhysicalVector` or a `PF.PhysicalTensor`, depending upon the type of `vie`, while `cₙ` will always be a `PF.PhysicalScalar`.

## References

1. Braß, H., "On the Principle of Avoiding the Singularity in Quadrature," *Zeitschrift für angewandte Mathematik und Mechanik*, **75** (1995), S617-S618.
2. Caputo, M. and Mainardi, F., "Linear models of dissipation in anelastic solids," *Rivista del Nuoro Cimento*, **1** (1971), 161-198.
3. Caputo, M. and Mainardi, F., "A new dissipation model based on memory mechanism," *Pure and Applied Geophysics*, **91** (1971), 134-147.
4. Cole, K.S. and Cole, R.H., "Dispersion and absorption in dielectrics I. Alternating current characteristics," *Journal of Chemical Physics*, **9** (1941), 342-351.
5. Cole, K.S. and Cole, R.H., "Dispersion and absorption in dielectrics II. Direct current characteristics," *Journal of Chemical Physics*, **10** (1942), 98-105.
6. Freed, A.D., Soft Solids: A primer to the theoretical mechanics of materials, Modeling and Simulation. In Science, Engineering and Technology. Basel: Birkhäuser, 2014.
7. Freed, A.D. and Rajagopal, K.R. "A viscoelastic model for describing the response of biological fibers," *ACTA Mechanica*, **227** (2016), 3367-3380.
8. Fung, Y.-C., "Biorheology of Soft Tissues," *Biorheology*, **10** (1973), 139-155.
9. Kohlrausch, R., "Ueber das Dellmann'sche Elektrometer," *Annalen der Physik und Chemie*, **72** (1847), 353-405.
10. Maxwell, J.C., "On the dynamical theory of gases," *Philosophical Transactions of the Royal Society*, London, **157** (1867), 49-88.
11. Neubert, H.K., "A simple model representing internal damping in solid materials," *The Aeronautical Quarterly*, **14** (1963), 187-210.
12. Volterra, V. Theory of functionals and of integral and integro-differentialequations. Glasgow: Blackie and Son, 1930.
13. Williams, G. and Watts, D.C., "Non-symmetrical dielectric relaxation behaviour arising from a simple empirical decay function," *Transactions of the Faraday Society*, **66** (1970), 80-85.
14. Williams, M.L., "Structural analysis of viscoelastic materials," *AIAA Journal*, **2** (1964), 785-808.
15. Young, A., "Approximate product integration," *Proceedings of the Royal Society*, London, **A-224** (1954), 552-561.
16. Zener, C., Elasticity and Anelasticity of Metals. Chicago: University of Chicago Press, 1948.
"""
module VolterraIntegralEquations

using
    JSON3,
    MittagLeffler,
    PhysicalFields,
    StructTypes
    
import
    PhysicalFields as PF

export
    # Memory functions, which are kernels to Volterra integral equations.
    BOX,    # the 𝑏𝑜𝑥 model of Neuber, a.k.a. Fung's QLV kernel
    CCM,    # 𝐶ole and 𝐶ole's power-law 𝑀odel
    FLS,    # Caputo and Mainardi's 𝐹ractional 𝐿inear 𝑆olid
    KWW,    # 𝐾ohlrausch's and 𝑊illiams & 𝑊atts' stretched exponential
    MCM,    # 𝑀axwell's 𝐶hain 𝑀odel, a.k.a. the Prony series model
    MPL,    # Williams' 𝑀odified 𝑃ower-𝐿aw model
    RFS,    # Freed and Rajagopal's 𝑅egularized 𝐹L𝑆 model
    SLS,    # Zener's 𝑆tandard 𝐿inear 𝑆olid, a.k.a. the Maxwell-Debye kernel

    # Function to create weights of quadrature for a given memory function.
    normalizedQuadratureWeights,

    # Solver for Volterra integral equations of the second kind; specifically,
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
    toFile,
    fromFile,

    advance!,
    update!

#=
-------------------------------------------------------------------------------
=#

"""
Memory function BOX (the box energy function of Neubert)

    k = (exp(-t/τ₂) - exp(-t/τ₁)) / (t ln(τ₂/τ₁))
    τ = τ₁

is described by
```julia
(name, k, τ) = BOX(system_of_units, time, parameters)
```
with arguments

    `system_of_units` is, at present, either "SI" or "CGS".
    `time` is an instance of `PhysicalScalar`.
    `parameters` describes the tuple (τ₁, τ₂) ordered as 0 < τ₁ < τ₂..
    
where

    `name` is returned as "BOX".
"""
function BOX(system_of_units::String,
             time::PhysicalScalar,
             parameters::Tuple)::Tuple
    # Ensure that a consistent system of physical units is used.
    if system_of_units == "SI" || system_of_units == "si"
        t  = toSI(time)
        τ₁ = toSI(parameters[1])
        τ₂ = toSI(parameters[2])
        if t.units ≠ SECOND || τ₁.units ≠ SECOND || τ₂.units ≠ SECOND
            msg = "Argument time and parameters τ₁ and τ₂ "
            msg = string(msg, "must have units of time.")
            error(msg)
        end
    elseif system_of_units == "CGS" || system_of_units == "cgs"
        t  = toCGS(time)
        τ₁ = toCGS(parameters[1])
        τ₂ = toCGS(parameters[2])
        if (t.units ≠ CGS_SECOND || τ₁.units ≠ CGS_SECOND ||
            τ₂.units ≠ CGS_SECOND)
            msg = "Argument time and parameter τ₁ and τ₂ "
            msg = string(msg, "must have units of time.")
            error(msg)
        end
    else
        error("The assigned physical system of units is unknown.")
    end

    if τ₁.value > 0.0 && τ₂ > τ₁
        if t.value ≈ 0.0
            k = (1/τ₁ - 1/τ₂) / log(τ₂/τ₁)
        elseif t.value > 0.0
            k = (exp(-t/τ₂) - exp(-t/τ₁)) / (t * log(τ₂/τ₁))
        else
            error("Argument time must be non-negative.")
        end
    else
        error("Parameters τ₁ and τ₂ must be ordered as 0 < τ₁ < τ₂.")
    end
    return ("BOX", k, τ₁)
end # BOX

"""
Memory function CCM (Cole-Cole power-law Model)

    k = (t/τ)^α (α / t) / (1 + (t/τ)^α)²
    τ = τ

is described by
```julia
(name, k, τ) = CCM(system_of_units, time, parameters)
```
with arguments

    `system_of_units` is, at present, either "SI" or "CGS".
    `time` is an instance of `PhysicalScalar`.
    `parameters` describes the tuple (α, τ).
    
where

    `name` is returned as "CCM".
"""
function CCM(system_of_units::String, 
             time::PhysicalScalar, 
             parameters::Tuple)::Tuple
    # Ensure that a consistent system of physical units is used.
    if system_of_units == "SI" || system_of_units == "si"
        t = toSI(time)
        α = toSI(parameters[1])
        τ = toSI(parameters[2])
        if t.units ≠ SECOND || α.units ≠ DIMENSIONLESS || τ.units ≠ SECOND
            msg = "Argument time and parameter τ must have units of time, "
            msg = string(msg, "while parameter α must be dimensionless.")
            error(msg)
        end
    elseif system_of_units == "CGS" || system_of_units == "cgs"
        t = toCGS(time)
        α = toCGS(parameters[1])
        τ = toCGS(parameters[2])
        if (t.units ≠ CGS_SECOND || α.units ≠ CGS_DIMENSIONLESS ||
            τ.units ≠ CGS_SECOND)
            msg = "Argument time and parameter τ must have units of time, "
            msg = string(msg, "while parameter α must be dimensionless.")
            error(msg)
        end
    else
        error("The assigned physical system of units is unknown.")
    end

    if α.value > 0.0 && τ.value > 0.0
        if t.value == 0.0
            k = PhysicalScalar(-t.units)
            set!(k, Inf)
        elseif t.value > 0.0
            x = (t/τ)^get(α)
            k = x * (α/t) / ((1 + x) * (1 + x))
        else
            error("Argument time must be non-negative.")
        end
    else
        error("Parameters α and τ must be positive.")
    end
    return ("CCM", k, τ)
end # CCM

"""
Memory function FLS (Fractional Linear Solid)

    k = -E_{α,0}(-(t/τ)^α) / t
    τ = τ

where E_{α,β}(z) denotes the two-parameter Mittag-Leffler function, with β = 0 in this case. This kernel follows from 
```julia
(name, k, τ) = FLS(system_of_units, time, parameters)
```
with arguments

    `system_of_units` is, at present, either "SI" or "CGS".
    `time` is an instance of `PhysicalScalar`.
    `parameters` describes the tuple (α, τ).
    
where

    `name` is returned as "FLS".
"""
function FLS(system_of_units::String, time::PhysicalScalar, parameters::Tuple)::Tuple
    # Ensure that a consistent system of physical units is used.
    if system_of_units == "SI" || system_of_units == "si"
        t = toSI(time)
        α = toSI(parameters[1])
        τ = toSI(parameters[2])
        if t.units ≠ SECOND || α.units ≠ DIMENSIONLESS || τ.units ≠ SECOND
            msg = "Argument time and parameter τ must have units of time, "
            msg = string(msg, "while parameter α must be dimensionless.")
            error(msg)
        end
    elseif system_of_units == "CGS" || system_of_units == "cgs"
        t = toCGS(time)
        α = toCGS(parameters[1])
        τ = toCGS(parameters[2])
        if (t.units ≠ CGS_SECOND || α.units ≠ CGS_DIMENSIONLESS ||
            τ.units ≠ CGS_SECOND)
            msg = "Argument time and parameter τ must have units of time, "
            msg = string(msg, "while parameter α must be dimensionless.")
            error(msg)
        end
    else
        error("The assigned physical system of units is unknown.")
    end

    if α.value ≈ 1.0
        (k, τ) = SLS(system_of_units, time, (τ,))
    elseif τ.value > 0.0 && α.value > 0.0 && α.value < 1.0
        if t.value == 0.0
            k = PhysicalScalar(-t.units)
            set!(k, Inf)
        elseif t.value > 0.0
            x = (t / τ)^get(α)
            k = -MittagLeffler.mittleff(get(α), 0.0, -get(x)) / t
        else
            error("Argument time must be non-negative.")
        end
    else
        error("Parameter τ must be positive, and parameter α ∈ (0,1].")
    end
    return ("FLS", k, τ)
end # FLS

"""
Memory function KWW (stretched exponential of Kohlrausch, Williams and Watts)

    k = (t/τ)^α (α/t) exp(-(t/τ)^α)
    τ = τ
    
is described by
```julia
(name, k, τ) = KWW(system_of_units, time, parameters)
```
with arguments

    `system_of_units` is, at present, either "SI" or "CGS".
    `time` is an instance of `PhysicalScalar`.
    `parameters` describes the tuple (α, τ).
    
where

    `name` is returned as "KWW".
"""
function KWW(system_of_units::String, 
             time::PhysicalScalar, 
             parameters::Tuple)::Tuple
    # Ensure that a consistent system of physical units is used.
    if system_of_units == "SI" || system_of_units == "si"
        t = toSI(time)
        α = toSI(parameters[1])
        τ = toSI(parameters[2])
        if t.units ≠ SECOND || α.units ≠ DIMENSIONLESS || τ.units ≠ SECOND
            msg = "Argument time and parameter τ must have units of time, "
            msg = string(msg, "while parameter α must be dimensionless.")
            error(msg)
        end
    elseif system_of_units == "CGS" || system_of_units == "cgs"
        t = toCGS(time)
        α = toCGS(parameters[1])
        τ = toCGS(parameters[2])
        if (t.units ≠ CGS_SECOND || α.units ≠ CGS_DIMENSIONLESS ||
            τ.units ≠ CGS_SECOND)
            msg = "Argument time and parameter τ must have units of time, "
            msg = string(msg, "while parameter α must be dimensionless.")
            error(msg)
        end
    else
        error("The assigned physical system of units is unknown.")
    end

    if α.value ≈ 1.0
        (k, τ) = SLS(system_of_units, time, (τ,))
    elseif τ.value > 0.0 && α.value > 0.0 && α.value < 1.0
        if t.value == 0.0
            k = PhysicalScalar(-t.units)
            set!(k, Inf)
        elseif t.value > 0.0
            k = (t/τ)^α.value * (α/t) * exp(-(t/τ)^α.value)
        else
            error("Argument time must be non-negative.")
        end
    else
        error("Parameter τ must be positive and parameter α ∈ (0,1].")
    end
    return ("KWW", k, τ)
end # KWW

"""
Memory function MCM (Maxwell Chain Model, which is a Prony series)

    k = (c₁/τ₁) exp(-t/τ₁) + ⋯ + (cₙ/τₙ) exp(-t/τₙ)
    τ = τ₁

is described by
```julia
(name, k, τ) = MCM(system_of_units, time, parameters)
```
with arguments

    `system_of_units` is, at present, either "SI" or "CGS".
    `time` is an instance of `PhysicalScalar`.
    `parameters` describes a tuple (c₁, c₂, …, cₙ, τ₁, τ₂, …, τₙ) of length 2n, 
         where ∑ᵢ₌₁ⁿ cᵢ = 1, and where 0 < τ₁ < τ₂ < ⋯ < τₙ.
    
where

    `name` is returned as "MCM".
"""
function MCM(system_of_units::String, 
             time::PhysicalScalar, 
             parameters::Tuple)::Tuple
    # Verify the inputs.
    if length(parameters) % 2 == 0
        n = length(parameters) ÷ 2
    else
        error("There must be an even number of parameters, viz. paired sets.")
    end

    # Ensure that a consistent system of physical units is used.
    if system_of_units == "SI" || system_of_units == "si"
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
            error("Prony series coefficients must sum to 1.")
        end
        for i in 2:n
            if τ[i-1] ≥ τ[i]
                error("Prony characteristic times must order 0 < τ₁ < ⋯ < τₙ.")
            end
        end
    elseif system_of_units == "CGS" || system_of_units == "cgs"
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
            error("Prony series coefficients must sum to 1.")
        end
        for i in 2:n
            if τ[i-1] ≥ τ[i]
                error("Prony characteristic times must order 0 < τ₁ < ⋯ < τₙ.")
            end
        end
    else
        error("The assigned physical system of units is unknown.")
    end

    if t.value ≥ 0.0 && τ[1].value > 0.0
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
        error("Arguments time and τᵢ must be positive.")
    end
    return ("MCM", k, τ₁)
end # MCM

"""
Memory function MPL (Modified Power-Law model of Williams')

    k = (α/τ) / (1 + t/τ)^(1+α)
    τ = τ

is described by
```julia
(name, k, τ) = MPL(system_of_units, time, parameters)
```
with arguments

    `system_of_units` is, at present, either "SI" or "CGS".
    `time` is an instance of `PhysicalScalar`.
    `parameters` describes the tuple (α, τ).
    
where

    `name` is returned as "MPL".
"""
function MPL(system_of_units::String, 
             time::PhysicalScalar, 
             parameters::Tuple)::Tuple
    # Ensure that a consistent system of physical units is used.
    if system_of_units == "SI" || system_of_units == "si"
        t = toSI(time)
        α = toSI(parameters[1])
        τ = toSI(parameters[2])
        if t.units ≠ SECOND || α.units ≠ DIMENSIONLESS || τ.units ≠ SECOND
            msg = "Argument time and parameter τ must have units of time, "
            msg = string(msg, "while parameter α must be dimensionless.")
            error(msg)
        end
    elseif system_of_units == "CGS" || system_of_units == "cgs"
        t = toCGS(time)
        α = toCGS(parameters[1])
        τ = toCGS(parameters[2])
        if (t.units ≠ CGS_SECOND || α.units ≠ CGS_DIMENSIONLESS ||
            τ.units ≠ CGS_SECOND)
            msg = "Argument time and parameter τ must have units of time, "
            msg = string(msg, "while parameter α must be dimensionless.")
            error(msg)
        end
    else
        error("The assigned physical system of units is unknown.")
    end

    if (t.value ≥ 0.0) && (α.value > 0.0) && (τ.value > 0.0)
        k = (α/τ) / (1 + t/τ)^(1+α.value)
    else
        error("Arguments time, α and τ must be positive.")
    end
    return ("MPL", k, τ)
end # MPL

"""
Memory function RFS (Regularized Fractional Solid)

    k = -E_{α,0}(-((t+δ)/τ)^α) / (E_{α,1}(-(δ/τ)^α)(t+δ))
    τ = τ

where E_{α, β}(z) denotes the two-parameter Mittag-Leffler function. This kernel follows from
```julia
(name, k, τ) = RFS(system_of_units, time, parameters)
```
with arguments

    `system_of_units` is, at present, either "SI" or "CGS".
    `time` is an instance of `PhysicalScalar`.
    `parameters` describes the tuple (α, δ, τ).
    
where

    `name` is returned as "RFS".
"""
function RFS(system_of_units::String, 
             time::PhysicalScalar, 
             parameters::Tuple)::Tuple
    # Ensure that a consistent system of physical units is used.
    if system_of_units == "SI" || system_of_units == "si"
        t = toSI(time)
        α = toSI(parameters[1])
        δ = toSI(parameters[2])
        τ = toSI(parameters[3])
        if (t.units ≠ SECOND || α.units ≠ DIMENSIONLESS ||
            δ.units ≠ SECOND || τ.units ≠ SECOND)
            msg = "Argument time and parameters τ and δ must have units of "
            msg = string(msg, "time, while parameter α must be dimensionless.")
            error(msg)
        end
    elseif system_of_units == "CGS" || system_of_units == "cgs"
        t = toCGS(time)
        α = toCGS(parameters[1])
        δ = toCGS(parameters[2])
        τ = toCGS(parameters[3])
        if (t.units ≠ CGS_SECOND || α.units ≠ CGS_DIMENSIONLESS ||
            δ.units ≠ CGS_SECOND || τ.units ≠ CGS_SECOND)
            msg = "Argument time and parameters τ and δ must have units of "
            msg = string(msg, "time, while parameter α must be dimensionless.")
            throw(ErrorException(msg))
        end
    else
        error("The assigned physical system of units is unknown.")
    end

    if α.value ≈ 1.0
        (k, τ) = SLS(system_of_units, time, (τ,))
    elseif ((t.value ≥ 0.0) && (α.value > 0.0) && (α.value < 1.0) &&
        (δ.value > 0.0) && (τ.value > 0.0))
        x = ((t + δ) / τ)^get(α)
        numerMLF = MittagLeffler.mittleff(get(α), 0.0, -get(x))
        y = (δ / τ)^get(α)
        denomMLF = MittagLeffler.mittleff(get(α), 1.0, -get(y))
        k = -numerMLF / (denomMLF * (δ + t))
    else
        msg = "Argument time must be non-negative, while parameters "
        msg = string(msg, "τ and δ must be positive, and α ∈ (0,1].")
        throw(ErrorException(msg))
    end
    return ("RFS", k, τ)
end # RFS

"""
Memory function SLS (Standard Linear Solid)

    k = exp(-t/τ) / τ
    τ = τ

is described by
```julia
(name, k, τ) = SLS(system_of_units, time, parameters)
```
with arguments

    `system_of_units` is, at present, either "SI" or "CGS".
    `time` is an instance of `PhysicalScalar`.
    `parameters` is described by the tuple (τ,).
    
where

    `name` is returned as "SLS".
"""
function SLS(system_of_units::String, 
             time::PhysicalScalar, 
             parameters::Tuple)::Tuple
    # Ensure that a consistent system of physical units is used.
    if system_of_units == "SI" || system_of_units == "si"
        t = toSI(time)
        τ = toSI(parameters[1])
        if t.units ≠ SECOND || τ.units ≠ SECOND
            error("Argument time and parameter τ must have units of time.")
        end
    elseif system_of_units == "CGS" || system_of_units == "cgs"
        t = toCGS(time)
        τ = toCGS(parameters[1])
        if t.units ≠ CGS_SECOND || τ.units ≠ CGS_SECOND
            error("Argument time and parameter τ must have units of time.")
        end
    else
        error("The assigned physical system of units is unknown.")
    end

    if (t.value ≥ 0.0) && (τ.value > 0.0)
        k = exp(-t/τ) / τ
    else
        error("time must be non-negative, and τ must be positive.")
    end
    return ("SLS", k, τ)
end # SLS

#=
-------------------------------------------------------------------------------
=#

# Weights *w* and nodes *ξ* for a Gaussian quadrature of order *S*. Gaussian
# quadrature is used to approximate integrals when deriving the weights of
# quadrature used to approximate Volterra integral equations. Quadrature rules
# for `S ∈ (1, 2, …, 5)` are provided.

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
Function
```julia
W = normalizedQuadratureWeights(system_of_units, N, dtime, kernel, parameters)
```
where, at present, `system_of_units` is either "SI" or "CGS". There are to be `N` intervals of size `dtime` that are to span a solution, whose `kernel` has `parameters` described via a tuple of material constants. These weights are written to a file in the user's **./files/** directory for efficient future use.

The supplied memory function `kernel` is to have an interface of
```julia
(name, k, τ) = kernel(system_of_units, time, parameters)
```
where `system_of_units` is either "SI" or "CGS", `time` is current time, and `parameters` is a tuple containing this kernel's physical parameters, i.e., its material constants. The returned tuple contains a string specifying the `name` or acronym of the kernel being evaluated, the value `k` of kernel `K` being evaluated at `time`, and its characteristic time `τ`.

The weights of quadrature returned here are normalized, e.g., the actual weights of quadrature for a linear viscoelastic kernel would be these normalized weights multiplied by a scalar coefficient of `(E₀ - E∞)/E∞`, which is to be supplied via a function call assigned to field `c` in an object implementing abstract type `VolterraIntegralEquation.`

The returned array holds normalized quadrature weights that are to be assigned to field `W` in an object implementing abstract type `VolterraIntegralEquation`.
"""
function normalizedQuadratureWeights(system_of_units::String, 
                                     N::Int, 
                                     dtime::PhysicalScalar, 
                                     kernel::Function,
                                     parameters::Tuple)::ArrayOfPhysicalScalars

    # Ensure the system of units is consistent.
    if system_of_units == "SI" || system_of_units == "si"
        dt = toSI(dtime)
        units = "SI"
    elseif system_of_units == "CGS" || system_of_units == "cgs"
        dt = toCGS(dtime)
        units = "CGS"
    else
        error("The assigned physical system of units is unknown.")
    end

    # Verify the inputs.
    if N < 1
        error("The number of integration steps N must be positive.")
    end
    if get(dt) < eps(Float32)
        error("The time step size dTime must be positive.")
    end

    # Check if quadrature weights have been previously calculated or not.
    (fileName, k, τ) = kernel(system_of_units, dt, parameters)
    if !isDimensionless(dt.units+k.units)
        error("Units for dtime and for the kernel are not compatible.")
    end
    fileName = string(fileName, "_", units)
    fileName = string(fileName, "_τ=", PF.toString(get(τ)))
    fileName = string(fileName, "_dt=", PF.toString(get(dt)))
    fileName = string(fileName, "_N=", N)
    fileName = string(fileName, ".json")
    dirPath  = string(pwd(), "/files/")
    if !isdir(dirPath)
        mkdir(dirPath)
    end
    my_file = string(dirPath, fileName)
    if isfile(my_file)
        # The quadrature weights exist. Read them in from a file.
        json_stream = PF.openJSONReader(dirPath, fileName)
        quadWgts = PF.fromFile(ArrayOfPhysicalScalars, json_stream)
        PF.closeJSONStream(json_stream)
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
            (name, kₛ, τₛ) = kernel(system_of_units, tₛ, parameters)
            sum = sum + quad.w[s] * kₛ
        end
        quadWgts[n] = 0.5dt * sum
    end

    # Write these weights of quadrature to a file.
    json_stream = PF.openJSONWriter(dirPath, fileName)
    PF.toFile(quadWgts, json_stream)
    PF.closeJSONStream(json_stream)
    println("Weights of quadrature have been written to a file.")

    return quadWgts
end # noralizedQuadratureWeights

#=
-------------------------------------------------------------------------------
=#

"""
The super type of all types that solve Volterra integral equations.
```julia
abstract type VolterraIntegralEquation end
```
"""
abstract type VolterraIntegralEquation end

# Scalar-valued Volterra integral equations of the second kind.

"""
```julia
struct VolterraIntegralScalarEquation <: VolterraIntegralEquation
    # Dimensioning fields
    dt::PhysicalScalar          # distance separating solution neighbor nodes
    N::Int                      # number of integration nodes in solution path
    n::MInteger                 # current node along a solution path
    # Arrays of length N+1 holding control and response fields, and their rates
    f::ArrayOfPhysicalScalars   # array of integrated response function values
    f′::ArrayOfPhysicalScalars  # array of response function rates
    g::ArrayOfPhysicalScalars   # array of integrated control function values
    g′::ArrayOfPhysicalScalars  # array of control function rates
    t::ArrayOfPhysicalScalars   # array of times, the independent variable
    # Array of N normalized weights of quadrature for a product integral
    W::ArrayOfPhysicalScalars   # array holding the quadrature weights
end
```
with constructor
```julia
vie = VolterraIntegralScalarEquation(system_of_units, N, dt, f₀, g₀, W)
```
where
    *vie* is a new instance of type `VolterraIntegralScalarEquation`.
    *system_of_units* is a string, viz., either "SI" or "CGS".
    *N*  is the number of nodes of quadrature to be imposed.
    *dt* is a step size separating all neighboring nodes--a `PhyscalScalar`.
    *f₀* is an initial condition for the response function--a `PhysicalScalar`.
    *g₀* is an initial condition for the control function--a `PhysicalScalar`.
    *W*  holds the weights of quadrature--an `ArrayOfPhysicalScalars`.
"""
struct VolterraIntegralScalarEquation <: VolterraIntegralEquation
    # Dimensioning fields
    dt::PhysicalScalar          # distance separating solution neighbor nodes
    N::Int                      # number of integration nodes in solution path
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
    function VolterraIntegralScalarEquation(system_of_units::String, 
                                            N::Int, 
                                            dt::PhysicalScalar,
                                            f₀::PhysicalScalar,
                                            g₀::PhysicalScalar,
                                            W::ArrayOfPhysicalScalars)

        # Ensure that a consistent system of physical units is used.
        if system_of_units == "SI" || system_of_units == "si"
            dt = toSI(dt)
            f₀ = toSI(f₀)
            g₀ = toSI(g₀)
            W  = toSI(W)
        elseif system_of_units == "CGS" || system_of_units == "cgs"
            dt = toCGS(dt)
            f₀ = toCGS(f₀)
            g₀ = toCGS(g₀)
            W  = toCGS(W)
        else
            error("The assigned physical system of units is unknown.")
        end
        t₀ = PhysicalScalar(dt.units)

        # Verify the remaining inputs.
        if N < 1
            error("The number of nodes N must be positive valued.")
        end
        if f₀.units ≠ g₀.units
            error("Units for initial conditions f₀ and g₀ must be equal.")
        end
        if !isDimensionless(W)
            error("Weights of quadrature W must be dimensionless.")
        end

        # Create the fields for this data structure,
        f = ArrayOfPhysicalScalars(N+1, f₀.units)
        f[1] = f₀
        f′ = ArrayOfPhysicalScalars(N+1, f₀.units-dt.units)
        g  = ArrayOfPhysicalScalars(N+1, g₀.units)
        g[1] = g₀
        g′ = ArrayOfPhysicalScalars(N+1, g₀.units-dt.units)
        t  = ArrayOfPhysicalScalars(N+1, t₀.units)
        t[1] = t₀
        for n in 1:N
            t[n+1] = n * dt
        end
        n = MInteger(1)

        new(dt, N, n, f, f′, g, g′, t, W)::VolterraIntegralScalarEquation
    end

    # Used by JSON3 whenever this data structure is to be created from a file.
    function VolterraIntegralScalarEquation(dt::PhysicalScalar, 
                                            N::Int, 
                                            n::MInteger,
                                            f::ArrayOfPhysicalScalars,
                                            f′::ArrayOfPhysicalScalars,
                                            g::ArrayOfPhysicalScalars,
                                            g′::ArrayOfPhysicalScalars,
                                            t::ArrayOfPhysicalScalars,
                                            W::ArrayOfPhysicalScalars)

        new(dt, N, n, f, f′, g, g′, t, W)::VolterraIntegralScalarEquation
    end
end # VolterraIntegralScalarEquation

# Methods

function Base.:(copy)(vie::VolterraIntegralScalarEquation)::VolterraIntegralScalarEquation
    dt = copy(vie.dt)
    N  = copy(vie.N)
    n  = copy(vie.n)
    f  = copy(vie.f)
    f′ = copy(vie.f′)
    g  = copy(vie.g)
    g′ = copy(vie.g′)
    t  = copy(vie.t)
    W  = copy(vie.W)
    return VolterraIntegralScalarEquation(dt, N, n, f, f′, g, g′, t, W)
end

StructTypes.StructType(::Type{VolterraIntegralScalarEquation}) = StructTypes.Struct()

function toFile(vie::VolterraIntegralScalarEquation, json_stream::IOStream)
    if isopen(json_stream)
        JSON3.write(json_stream, vie)
        write(json_stream, '\n')
    else
        error("The supplied JSON stream is not open.")
    end
    flush(json_stream)
    return nothing
end

function fromFile(::Type{VolterraIntegralScalarEquation}, json_stream::IOStream)::VolterraIntegralScalarEquation
    if isopen(json_stream)
        vie = JSON3.read(readline(json_stream), VolterraIntegralScalarEquation)
    else
        error("The supplied JSON stream is not open.")
    end
    return vie
end

"""
Function
```julia
function advance!(vie::VolterraIntegralScalarEquation, 
                  g′ₙ::PhysicalScalar, 
                  cₙ::PhysicalScalar)
```
advances a solution step-by-step, i.e., from step n-1 to step n.
"""
function advance!(vie::VolterraIntegralScalarEquation, 
                  g′ₙ::PhysicalScalar, 
                  cₙ::PhysicalScalar)
    if vie.n > vie.N
        println("The Volterra integral solution has reached its endpoint.")
        return nothing
    end

    # verify inputs
    if g′ₙ.units ≠ vie.f′.units
        msg = "Physical units for g′ₙ must equal those of vie.f′.\n"
        msg = string(msg, "   g′ₙ has units ", PF.toString(g′ₙ.units), "\n")
        msg = string(msg, "   f′  has units ", PF.toString(vie.f′.units))
        error(msg)
    end
    if !isDimensionless(cₙ)
        error("Coefficient cₙ must be dimensionless scalar.")
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
        vie.f[3] = (4vie.f[2] - vie.f[1] + 2vie.f′[3]*vie.dt) / 3
        vie.g[3] = (4vie.g[2] - vie.g[1] + 2vie.g′[3]*vie.dt) / 3
    else
        vie.f[n] = ((18vie.f[n-1] - 9vie.f[n-2] + 2vie.f[n-3]
            + 6vie.f′[n]*vie.dt) / 11)
        vie.g[n] = ((18vie.g[n-1] - 9vie.g[n-2] + 2vie.g[n-3]
            + 6vie.g′[n]*vie.dt) / 11)
    end

    return nothing
end # advance!

"""
Function
```julia
function update!(vie::VolterraIntegralScalarEquation, 
                 g′ₙ::PhysicalScalar,
                 cₙ::PhysicalScalar)
```
performs an iteration of refinement on a solution at current step n. Call only if the control function *g′ₙ* or coefficient *cₙ* undergo iterative refinement.
"""
function update!(vie::VolterraIntegralScalarEquation, 
                 g′ₙ::PhysicalScalar,
                 cₙ::PhysicalScalar)

    set!(vie.n, get(vie.n)-1)
    advance!(vie, g′ₙ, cₙ)
    return nothing
end # update!

#=
-------------------------------------------------------------------------------
=#

# Vector-valued Volterra integral equations of the second kind.

"""
```julia
struct VolterraIntegralVectorEquation <: VolterraIntegralEquation
    # Dimensioning fields
    dt::PhysicalScalar          # distance separating solution neighbor nodes
    N::Int                      # number of integration nodes in solution path
    n::MInteger                 # current node along a solution path
    # Arrays of length N+1 holding control and response fields, and their rates
    f::ArrayOfPhysicalVectors   # array of integrated response function values
    f′::ArrayOfPhysicalVectors  # array of response function rates
    g::ArrayOfPhysicalVectors   # array of integrated control function values
    g′::ArrayOfPhysicalVectors  # array of control function rates
    t::ArrayOfPhysicalScalars   # array of times, the independent variable
    # Array of N normalized weights of quadrature for a product integral
    W::ArrayOfPhysicalScalars   # array holding the quadrature weights
end
```
with constructor
```julia
vie = VolterraIntegralVectorEquation(system_of_units, N, dt, f₀, g₀, W)
```
where
    *vie* is a new instance of type `VolterraIntegralVectorEquation`.
    *system_of_units* is a string, viz., either "SI" or "CGS".
    *N*  is the number of nodes of quadrature to be imposed.
    *dt* is a step size separating all neighboring nodes--a `PhyscalScalar`.
    *f₀* is an initial condition for the response function--a `PhysicalVector`.
    *g₀* is an initial condition for the control function--a `PhysicalVector`.
    *W*  holds the weights of quadrature--an `ArrayOfPhysicalScalars`.
"""
struct VolterraIntegralVectorEquation <: VolterraIntegralEquation
    # Dimensioning fields
    dt::PhysicalScalar          # distance separating solution neighbor nodes
    N::Int                      # number of integration nodes in solution path
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
    function VolterraIntegralVectorEquation(system_of_units::String, 
                                            N::Int, 
                                            dt::PhysicalScalar,
                                            f₀::PhysicalVector,
                                            g₀::PhysicalVector,
                                            W::ArrayOfPhysicalScalars)

        # Ensure that a consistent system of physical units is used.
        if system_of_units == "SI" || system_of_units == "si"
            dt = toSI(dt)
            f₀ = toSI(f₀)
            g₀ = toSI(g₀)
            W  = toSI(W)
        elseif system_of_units == "CGS" || system_of_units == "cgs"
            dt = toCGS(dt)
            f₀ = toCGS(f₀)
            g₀ = toCGS(g₀)
            W  = toCGS(W)
        else
            error("The assigned physical system of units is unknown.")
        end
        t₀ = PhysicalScalar(dt.units)

        # Verify the remaining inputs.
        if N < 1
            error("The number of nodes N must be positive valued.")
        end
        if f₀.units ≠ g₀.units
            error("Units for initial conditions f₀ and g₀ must be equal.")
        end
        if f₀.vector.len ≠ g₀.vector.len
            error("Initial condition vectors f₀ and g₀ must equal in length.")
        end
        if !isDimensionless(W)
            error("Weights of quadrature W must be dimensionless.")
        end

        # Create the fields for this data structure,
        f = ArrayOfPhysicalVectors(N+1, f₀.vector.len, f₀.units)
        f[1] = f₀
        f′ = ArrayOfPhysicalVectors(N+1, f₀.vector.len, f₀.units-dt.units)
        g  = ArrayOfPhysicalVectors(N+1, g₀.vector.len, g₀.units)
        g[1] = g₀
        g′ = ArrayOfPhysicalVectors(N+1, g₀.vector.len, g₀.units-dt.units)
        t  = ArrayOfPhysicalScalars(N+1, t₀.units)
        t[1] = t₀
        for n in 1:N
            t[n+1] = n * d𝑡
        end
        n = MInteger(1)

        new(dt, N, n, f, f′, g, g′, t, W)::VolterraIntegralVectorEquation
    end

    # Used by JSON3 whenever this data structure is to be created from a file.
    function VolterraIntegralVectorEquation(dt::PhysicalScalar, 
                                            N::Int, 
                                            n::MInteger,
                                            f::ArrayOfPhysicalVectors,
                                            f′::ArrayOfPhysicalVectors,
                                            g::ArrayOfPhysicalVectors,
                                            g′::ArrayOfPhysicalVectors,
                                            t::ArrayOfPhysicalScalars,
                                            W::ArrayOfPhysicalScalars)

        new(dt, N, n, f, f′, g, g′, t, W)::VolterraIntegralVectorEquation
    end
end # VolterraIntegralVectorEquation

# Methods

function Base.:(copy)(vie::VolterraIntegralVectorEquation)::VolterraIntegralVectorEquation
    dt = copy(vie.dt)
    N  = copy(vie.N)
    n  = copy(vie.n)
    f  = copy(vie.f)
    f′ = copy(vie.f′)
    g  = copy(vie.g)
    g′ = copy(vie.g′)
    t  = copy(vie.t)
    W  = copy(vie.W)
    return VolterraIntegralVectorEquation(dt, N, n, f, f′, g, g′, t, W)
end

StructTypes.StructType(::Type{VolterraIntegralVectorEquation}) = StructTypes.Struct()

function toFile(vie::VolterraIntegralVectorEquation, json_stream::IOStream)
    if isopen(json_stream)
        JSON3.write(json_stream, vie)
        write(json_stream, '\n')
    else
        error("The supplied JSON stream is not open.")
    end
    flush(json_stream)
    return nothing
end

function fromFile(::Type{VolterraIntegralVectorEquation}, json_stream::IOStream)::VolterraIntegralVectorEquation
    if isopen(json_stream)
        vie = JSON3.read(readline(json_stream), VolterraIntegralVectorEquation)
    else
        error("The supplied JSON stream is not open.")
    end
    return vie
end

"""
Function
```julia
function advance!(vie::VolterraIntegralVectorEquation, 
                  g′ₙ::PhysicalVector, 
                  cₙ::PhysicalScalar)
```
advances a solution step-by-step, i.e., from step n-1 to step n.
"""
function advance!(vie::VolterraIntegralVectorEquation, 
                  g′ₙ::PhysicalVector, 
                  cₙ::PhysicalScalar)
    if vie.n > vie.N
        println("The Volterra integral solution has reached its endpoint.")
        return nothing
    end

    # verify inputs
    if g′ₙ.units ≠ vie.f′.units
        msg = "Physical units for g′ₙ must equal those of vie.f′.\n"
        msg = string(msg, "   g′ₙ has units ", PF.toString(g′ₙ.units), "\n")
        msg = string(msg, "   f′  has units ", PF.toString(vie.f′.units))
        error(msg)
    end
    if g′ₙ.vector.len ≠ g′.array.cols
        error("Vector g′ₙ has the wrong length.")
    end
    if !isDimensionless(cₙ)
        error("Coefficient cₙ must be dimensionless scalar.")
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
            vie.f[3,j] = ((4vie.f[2,j] - vie.f[1,j]
                + 2vie.f′[3,j]*vie.dt) / 3)
            vie.g[3,j] = ((4vie.g[2,j] - vie.g[1,j]
                + 2vie.g′[3,j]*vie.dt) / 3)
        else
            vie.f[n,j] = ((18vie.f[n-1,j] - 9vie.f[n-2,j]
                + 2vie.f[n-3,j] + 6vie.f′[n,j]*vie.dt) / 11)
            vie.g[n,j] = ((18vie.g[n-1,j] - 9vie.g[n-2,j]
                + 2vie.g[n-3,j] + 6vie.g′[n,j]*vie.dt) / 11)
        end
    end

    return nothing
end # advance!

"""
Function
```julia
function update!(vie::VolterraIntegralVectorEquation, 
                 g′ₙ::PhysicalVector,
                 cₙ::PhysicalScalar)
```
performs an iteration of refinement on a solution at current step n. Call only if the control function *g′ₙ* or coefficient *cₙ* undergo iterative refinement.
"""
function update!(vie::VolterraIntegralVectorEquation, 
                 g′ₙ::PhysicalVector, 
                 cₙ::PhysicalScalar)

    set!(vie.n, get(vie.n)-1)
    advance!(vie, g′ₙ, cₙ)
    return nothing
end # update!

#=
-------------------------------------------------------------------------------
=#

# Tensor-valued Volterra integral equations of the second kind.

"""
```julia
struct VolterraIntegralTensorEquation <: VolterraIntegralEquation
    # Dimensioning fields
    dt::PhysicalScalar          # distance separating solution neighbor nodes
    N::Int                      # number of integration nodes in solution path
    n::MInteger                 # current node along a solution path
    # Arrays of length N+1 holding control and response fields, and their rates
    f::ArrayOfPhysicalTensors   # array of integrated response function values
    f′::ArrayOfPhysicalTensors  # array of response function rates
    g::ArrayOfPhysicalTensors   # array of integrated control function values
    g′::ArrayOfPhysicalTensors  # array of control function rates
    t::ArrayOfPhysicalScalars   # array of times, the independent variable
    # Array of N normalized weights of quadrature for a product integral
    W::ArrayOfPhysicalScalars   # array holding the quadrature weights
end
```
with constructor
```julia
vie = VolterraIntegralTensorEquation(system_of_units, N, dt, f₀, g₀, W)
```
where
    *vie* is a new instance of type `VolterraIntegralTensorEquation`.
    *system_of_units* is a string, viz., either "SI" or "CGS".
    *N*  is the number of nodes of quadrature to be imposed.
    *dt* is a step size separating all neighboring nodes--a `PhyscalScalar`.
    *f₀* is an initial condition for the response function--a `PhysicalTensor`.
    *g₀* is an initial condition for the control function--a `PhysicalTensor`.
    *W*  holds the weights of quadrature--an `ArrayOfPhysicalScalars`.
"""
struct VolterraIntegralTensorEquation <: VolterraIntegralEquation
    # Dimensioning fields
    dt::PhysicalScalar          # distance separating solution neighbor nodes
    N::Int                      # number of integration nodes in solution path
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
    function VolterraIntegralTensorEquation(system_of_units::String, 
                                            N::Int, 
                                            dt::PhysicalScalar,
                                            f₀::PhysicalTensor,
                                            g₀::PhysicalTensor,
                                            W::ArrayOfPhysicalScalars)

        # Ensure that a consistent system of physical units is used.
        if system_of_units == "SI" || system_of_units == "si"
            dt = toSI(dt)
            f₀ = toSI(f₀)
            g₀ = toSI(g₀)
            W  = toSI(W)
        elseif system_of_units == "CGS" || system_of_units == "cgs"
            dt= toCGS(dt)
            f₀ = toCGS(f₀)
            g₀ = toCGS(g₀)
            W  = toCGS(W)
        else
            error("The assigned physical system of units is unknown.")
        end
        t₀ = PhysicalScalar(dt.units)

        # Verify the remaining inputs.
        if N < 1
            error("The number of nodes N must be positive valued.")
        end
        if f₀.units ≠ g₀.units
            error("Units for initial conditions f₀ and g₀ must be equal.")
        end
        if (f₀.matrix.rows ≠ g₀.matrix.rows || 
            f₀.matrix.cols ≠ g₀.matrix.cols)
            msg = "Dimensions for the initial conditions, tensors f₀ and g₀, "
            msg = string(msg, "must equal.")
            error(msg)
        end
        if !isDimensionless(𝑊)
            error("Weights of quadrature W must be dimensionless.")
        end

        # Create the fields for this data structure,
        rows = f₀.matrix.rows
        cols = f₀.matrix.cols
        f  = ArrayOfPhysicalTensors(N+1, rows, cols, f₀.units)
        f[1] = f₀
        f′ = ArrayOfPhysicalTensors(N+1, rows, cols, f₀.units-dt.units)
        g  = ArrayOfPhysicalTensors(N+1, rows, cols, g₀.units)
        g[1] = g₀
        g′ = ArrayOfPhysicalTensors(N+1, rows, cols, g₀.units-dt.units)
        t  = ArrayOfPhysicalScalars(N+1, t₀.units)
        t[1] = t₀
        for n in 1:N
            t[n+1] = n * dt
        end
        n = MInteger(1)

        new(dt, N, n, f, f′, g, g′, t, W)::VolterraIntegralTensorEquation
    end

    # Used by JSON3 whenever this data structure is to be created from a file.
    function VolterraIntegralTensorEquation(dt::PhysicalScalar, 
                                            N::Int, 
                                            n::MInteger,
                                            f::ArrayOfPhysicalTensors,
                                            f′::ArrayOfPhysicalTensors,
                                            g::ArrayOfPhysicalTensors,
                                            g′::ArrayOfPhysicalTensors,
                                            t::ArrayOfPhysicalScalars,
                                            W::ArrayOfPhysicalScalars)

        new(dt, N, n, f, f′, g, g′, t, W)::VolterraIntegralTensorEquation
    end
end # VolterraIntegralTensorEquation

# Methods

function Base.:(copy)(vie::VolterraIntegralTensorEquation)::VolterraIntegralTensorEquation
    dt = copy(vie.dt)
    N  = copy(vie.N)
    n  = copy(vie.n)
    f  = copy(vie.f)
    f′ = copy(vie.f′)
    g  = copy(vie.g)
    g′ = copy(vie.g′)
    t  = copy(vie.t)
    W  = copy(vie.W)
    return VolterraIntegralTensorEquation(dt, N, n, f, f′, g, g′, t, W)
end

StructTypes.StructType(::Type{VolterraIntegralTensorEquation}) = StructTypes.Struct()

function toFile(vie::VolterraIntegralTensorEquation, json_stream::IOStream)
    if isopen(json_stream)
        JSON3.write(json_stream, vie)
        write(json_stream, '\n')
    else
        error("The supplied JSON stream is not open.")
    end
    flush(json_stream)
    return nothing
end

function fromFile(::Type{VolterraIntegralTensorEquation}, json_stream::IOStream)::VolterraIntegralTensorEquation
    if isopen(json_stream)
        vie = JSON3.read(readline(json_stream), VolterraIntegralTensorEquation)
    else
        error("The supplied JSON stream is not open.")
    end
    return vie
end

"""
Function
```julia
function advance!(vie::VolterraIntegralTensorEquation, 
                  g′ₙ::PhysicalTensor, 
                  cₙ::PhysicalScalar)
```
advances a solution step-by-step, i.e., from step n-1 to step n.
"""
function advance!(vie::VolterraIntegralTensorEquation, 
                  g′ₙ::PhysicalTensor, 
                  cₙ::PhysicalScalar)
    if vie.n > vie.N
        println("The Volterra integral solution has reached its endpoint.")
        return nothing
    end

    # verify inputs
    if g′ₙ.units ≠ vie.f′.units
        msg = "Physical units for g′ₙ must equal those of vie.f′.\n"
        msg = string(msg, "   g′ₙ has units ", PF.toString(g′ₙ.units), "\n")
        msg = string(msg, "   f′  has units ", PF.toString(vie.f′.units))
        error(msg)
    end
    if (g′ₙ.matrix.rows ≠ g′.array.rows) || (g′ₙ.matrix.cols ≠ g′.array.cols)
        error("Tensor g′ₙ has the wrong dimensions.")
    end
    if !isDimensionless(cₙ)
        error("Coefficient cₙ must be dimensionless scalar.")
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
                vie.f[3,j,k] = ((4vie.f[2,j,k] - vie.f[1,j,k]
                    + 2vie.f′[3,j,k]*vie.dt) / 3)
                vie.g[3,j,k] = ((4vie.g[2,j,k] - vie.g[1,j,k]
                    + 2vie.g′[3,j,k]*vie.dt) / 3)
            else
                vie.f[n,j,k] = ((18vie.f[n-1,j,k] - 9vie.f[n-2,j,k]
                    + 2vie.f[n-3,j,k] + 6vie.f′[n,j,k]*vie.dt) / 11)
                vie.g[n,j,k] = ((18vie.g[n-1,j,k] - 9vie.g[n-2,j,k]
                    + 2vie.g[n-3,j,k] + 6vie.g′[n,j,k]*vie.dt) / 11)
            end
        end
    end

    return nothing
end # advance!

"""
Function
```julia
function update!(vie::VolterraIntegralTensorEquation, 
                 g′ₙ::PhysicalTensor,
                 cₙ::PhysicalScalar)
```
performs an iteration of refinement on a solution at current step n. Call only if the control function *g′ₙ* or coefficient *cₙ* undergo iterative refinement.
"""
function update!(vie::VolterraIntegralTensorEquation, 
                 g′ₙ::PhysicalTensor, 
                 cₙ::PhysicalScalar)

    set!(vie.n, get(vie.n)-1)
    advance!(vie, g′ₙ, cₙ)
    return nothing
end # update!

end  # module VolterraIntegralEquations

