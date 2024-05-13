module testVolterraIntEqs

using
    CairoMakie,       # Pixel based figure construction.
    PhysicalFields,
    VolterraIntegralEquations

export
    Abel,
    memoryFns

function memoryFns(myDirPath::String)
    N = 150

    # viscoelastic material constants for elastin
    ν   = PhysicalScalar(0.59, CGS_DIMENSIONLESS)
    τ_ϵ = PhysicalScalar(9.0E-4, CGS_SECOND)
    τ_σ = PhysicalScalar(7.0E-8, CGS_SECOND)
    δ   = 0.1219 * τ_ϵ
    E_∞ = PhysicalScalar(2.3E6, BARYE)
    E_0 = (τ_ϵ/τ_σ)^get(ν) * E_∞

    # material constants for our model
    α = ν
    τ = τ_ϵ
    E = E_0
    c = (E_0 - E_∞) / E_∞

    # assumed material constants for the Maxwell Chain Model.
    c₁ = PhysicalScalar(0.25, CGS_DIMENSIONLESS)
    c₂ = PhysicalScalar(0.25, CGS_DIMENSIONLESS)
    c₃ = PhysicalScalar(0.25, CGS_DIMENSIONLESS)
    c₄ = PhysicalScalar(0.25, CGS_DIMENSIONLESS)
    τ₁ = τ_σ
    τ  = exp10(log10(get(τ_σ)) + (log10(get(τ_ϵ)) - log10(get(τ_σ)))/3)
    τ₂ = PhysicalScalar(τ, CGS_SECOND)
    τ  = exp10(log10(get(τ_σ)) + 2(log10(get(τ_ϵ)) - log10(get(τ_σ)))/3)
    τ₃ = PhysicalScalar(τ, CGS_SECOND)
    τ₄ = τ_ϵ
    println("For the Prony series, let c_i = 0.25, i = 1,...,4, and")
    println("   τ₁ = ", toString(τ₁))
    println("   τ₂ = ", toString(τ₂))
    println("   τ₃ = ", toString(τ₃))
    println("   τ₄ = ", toString(τ₄))
    println()

    # Create the arrays to hold values for the memory function.
    reciprocalTime = PhysicalUnits("CGS", 0, 0, 0, -1, 0, 0, 0)
    # weakly singular kernels
    arrayCCM = ArrayOfPhysicalScalars(N, reciprocalTime)
    arrayFLS = ArrayOfPhysicalScalars(N, reciprocalTime)
    arrayKWW = ArrayOfPhysicalScalars(N, reciprocalTime)
    # non-singular kernels
    arrayBOX = ArrayOfPhysicalScalars(N+1, reciprocalTime)
    arrayMCM = ArrayOfPhysicalScalars(N, reciprocalTime)
    arrayMPL = ArrayOfPhysicalScalars(N+1, reciprocalTime)
    arrayRFS = ArrayOfPhysicalScalars(N+1, reciprocalTime)
    arraySLS = ArrayOfPhysicalScalars(N+1, reciprocalTime)

    # Populate the memory function arrays
    time  = PhysicalScalar(CGS_SECOND)
    dTime = 3 * τ_ϵ / N
    (k, tau) = VolterraIntegralEquations.BOX("CGS", time, (τ_σ, τ_ϵ))
    arrayBOX[1] = k
    (k, tau) = VolterraIntegralEquations.MPL("CGS", time, (α, τ_ϵ))
    arrayMPL[1] = k
    (k, tau) = VolterraIntegralEquations.RFS("CGS", time, (α, δ, τ_ϵ))
    arrayRFS[1] = k
    (k, tau) = VolterraIntegralEquations.SLS("CGS", time, (τ_ϵ,))
    arraySLS[1] = k
    for n in 1:N
        time = time + dTime
        # weakly singular kernels
        (k, tau) = VolterraIntegralEquations.CCM("CGS", time, (α, τ_ϵ))
        arrayCCM[n] = k
        (k, tau) = VolterraIntegralEquations.FLS("CGS", time, (α, τ_ϵ))
        arrayFLS[n] = k
        (k, tau) = VolterraIntegralEquations.KWW("CGS", time, (α, τ_ϵ))
        arrayKWW[n] = k
        # non-singular kernels
        (k, tau) = VolterraIntegralEquations.BOX("CGS", time, (τ_σ, τ_ϵ))
        arrayBOX[n] = k
        (k, tau) = VolterraIntegralEquations.MCM("CGS", time, (c₁, c₂, c₃, c₄, τ₁, τ₂, τ₃, τ₄))
        arrayMCM[n] = k
        (k, tau) = VolterraIntegralEquations.MPL("CGS", time, (α, τ_ϵ))
        arrayMPL[n+1] = k
        (k, tau) = VolterraIntegralEquations.RFS("CGS", time, (α, δ, τ_ϵ))
        arrayRFS[n+1] = k
        (k, tau) = VolterraIntegralEquations.SLS("CGS", time, (τ_ϵ,))
        arraySLS[n+1] = k
    end

    # Create the arrays for plotting.
    t₁  = zeros(Float64, N)
    ccm = zeros(Float64, N)
    fls = zeros(Float64, N)
    kww = zeros(Float64, N)
    t₂  = zeros(Float64, N+1)
    box = zeros(Float64, N+1)
    mcm = zeros(Float64, N)
    mpl = zeros(Float64, N+1)
    rfs = zeros(Float64, N+1)
    sls = zeros(Float64, N+1)
    set!(time, 0.0)
    t₂[1]  = get(time) / get(τ_ϵ)
    box[1] = get(τ_ϵ) * get(arrayBOX[1])
    mcm[1] = get(τ_ϵ) * get(arrayMCM[1])
    mpl[1] = get(τ_ϵ) * get(arrayMPL[1])
    rfs[1] = get(τ_ϵ) * get(arrayRFS[1])
    sls[1] = get(τ_ϵ) * get(arraySLS[1])
    for n in 1:N
        # weakly singular kernels
        time   = time + dTime
        t₁[n]  = get(time) / get(τ_ϵ)
        ccm[n] = get(τ_ϵ) * get(arrayCCM[n])
        fls[n] = get(τ_ϵ) * get(arrayFLS[n])
        kww[n] = get(τ_ϵ) * get(arrayKWW[n])
        # non-singular kernels
        t₂[n+1]  = t₁[n]
        box[n]   = get(τ_ϵ) * get(arrayBOX[n+1])
        mcm[n]   = get(τ_ϵ) * get(arrayMCM[n])
        mpl[n+1] = get(τ_ϵ) * get(arrayMPL[n+1])
        rfs[n+1] = get(τ_ϵ) * get(arrayRFS[n+1])
        sls[n+1] = get(τ_ϵ) * get(arraySLS[n+1])
    end

    # Plot the non-singular kernels.

    fig1 = Figure(; size = (809, 500)) # (500ϕ, 500), ϕ is golden ratio
    ax1 = Axis(fig1[1, 1];
        xlabel = "time ÷ characteristic time  (t / τ_ϵ)",
        ylabel = "characteristic time × memory function  (τ_ϵ × k)",
        title = "Non-singular Memory Functions for Elastin",
        titlesize = 24,
        xlabelsize = 20,
        ylabelsize = 20)
    lines!(ax1, t₂, box;
        linewidth = 3,
        linestyle = :solid,
        color = :green,
        label = "BOX")
    lines!(ax1, t₁, mcm;
        linewidth = 3,
        linestyle = :solid,
        color = :cyan,
        label = "MCM")
    lines!(ax1, t₂, mpl;
        linewidth = 3,
        linestyle = :solid,
        color = :blue,
        label = "MPL")
    lines!(ax1, t₂, rfs;
        linewidth = 3,
        linestyle = :solid,
        color = :red,
        label = "RFS")
    lines!(ax1, t₂, sls;
        linewidth = 3,
        linestyle = :solid,
        color = :black,
        label = "SLS")
    axislegend("Function",
        position = :rt)
    save(string(myDirPath, "memoryFnNonSingular.png"), fig1)

    # Plot the weakly singular kernels.

    fig2 = Figure(; size = (809, 500)) # (500ϕ, 500), ϕ is golden ratio
    ax2 = Axis(fig2[1, 1];
        xlabel = "time ÷ characteristic time  (t / τ_ϵ)",
        ylabel = "characteristic time × memory function  (τ_ϵ × k)",
        title = "Weakly-singular Memory Functions for Elastin",
        titlesize = 24,
        xlabelsize = 20,
        ylabelsize = 20)
    lines!(ax2, t₁, ccm;
        linewidth = 3,
        linestyle = :solid,
        color = :blue,
        label = "CCM")
    lines!(ax2, t₁, fls;
        linewidth = 3,
        linestyle = :solid,
        color = :red,
        label = "FLS")
    lines!(ax2, t₁, kww;
        linewidth = 3,
        linestyle = :solid,
        color = :black,
        label = "KWW")
    axislegend("Function",
        position = :rt)
    save(string(myDirPath, "memoryFnWeaklySingular.png"), fig2)
end # memoryFns

function AbelKernel(systemOfUnits::String, time::PhysicalScalar, parameters::Tuple)::Tuple
    kernel = 1 / sqrt(time)
    tau = PhysicalScalar(1.0, CGS_DIMENSIONLESS)
    return ("Abel", kernel, tau)
end # AbelKernel

function Abel(myDirPath::String)
    CairoMakie.activate!(type = "png")

    # Solve an Abel integral equation: solution parameters.
    t = PhysicalScalar(10.0, CGS_DIMENSIONLESS) # upper limit of integration
    p = ()
    c = PhysicalScalar(1.0, CGS_DIMENSIONLESS)
    f₀ = PhysicalScalar(CGS_DIMENSIONLESS)
    g₀ = PhysicalScalar(CGS_DIMENSIONLESS)

    # Solve an Abel integral equation: 30 steps.
    println("Solving the Abel integral in 30 steps.")
    N₁ = 30
    dt₁ = t / N₁
    W₁ = VolterraIntegralEquations.normalizedQuadratureWeights("CGS", N₁, dt₁, AbelKernel, p)
    g′ₙ = PhysicalScalar(CGS_DIMENSIONLESS)
    VIE₁ = VolterraIntegralEquations.VolterraIntegralScalarEquation("CGS", N₁, dt₁, f₀, g₀, W₁)
    for n in 1:N₁
        tₙ = n * get(dt₁)
        set!(g′ₙ, π*tₙ/2 + sqrt(tₙ))
        VolterraIntegralEquations.advance!(VIE₁, g′ₙ, c)
    end

    # Create the arrays for plotting: 30 steps.
    x₁ = zeros(Float64, N₁+1)
    for n in 1:N₁+1
        x₁[n] = (n - 1) * get(dt₁)
    end
    y₁ = zeros(Float64, N₁+1)
    for n in 1:N₁+1
        soln = VIE₁.f[n]
        y₁[n] = get(soln)
    end
    z₁ = zeros(Float64, N₁+1)
    for n in 1:N₁+1
        z₁[n] = (2/3)*x₁[n]^(3/2)
    end
    e₁ = zeros(Float64, N₁)
    ϵ₁ = zeros(Float64, N₁)
    for n in 1:N₁
        e₁[n] = x₁[n+1]
        if abs(z₁[n+1]) < 1.0
            ϵ₁[n] = abs(z₁[n+1] - y₁[n+1])
        else
            ϵ₁[n] = abs(z₁[n+1] - y₁[n+1]) / abs(z₁[n+1])
        end
    end

    # Solve an Abel integral equation: 100 steps.
    println("Solving the Abel integral in 100 steps.")
    N₂ = 100
    dt₂ = t / N₂
    W₂ = VolterraIntegralEquations.normalizedQuadratureWeights("CGS", N₂, dt₂, AbelKernel, p)
    VIE₂ = VolterraIntegralEquations.VolterraIntegralScalarEquation("CGS", N₂, dt₂, f₀, g₀, W₂)
    for n in 1:N₂
        t₂ = n * get(dt₂)
        set!(g′ₙ, π*t₂/2 + sqrt(t₂))
        VolterraIntegralEquations.advance!(VIE₂, g′ₙ, c)
    end

    # Create the arrays for plotting: 100 steps.
    x₂ = zeros(Float64, N₂+1)
    for n in 1:N₂+1
        x₂[n] = (n - 1) * get(dt₂)
    end
    y₂ = zeros(Float64, N₂+1)
    for n in 1:N₂+1
        soln = VIE₂.f[n]
        y₂[n] = get(soln)
    end
    z₂ = zeros(Float64, N₂+1)
    for n in 1:N₂+1
        z₂[n] = (2/3)*x₂[n]^(3/2)
    end
    e₂ = zeros(Float64, N₂)
    ϵ₂ = zeros(Float64, N₂)
    for n in 1:N₂
        e₂[n] = x₂[n+1]
        if abs(z₂[n+1]) < 1.0
            ϵ₂[n] = abs(z₂[n+1] - y₂[n+1])
        else
            ϵ₂[n] = abs(z₂[n+1] - y₂[n+1]) / abs(z₂[n+1])
        end
    end

    # Solve an Abel integral equation: 300 steps.
    println("Solving the Abel integral in 300 steps.")
    N₃ = 300
    dt₃ = t / N₃
    W₃ = VolterraIntegralEquations.normalizedQuadratureWeights("CGS", N₃, dt₃, AbelKernel, p)
    VIE₃ = VolterraIntegralEquations.VolterraIntegralScalarEquation("CGS", N₃, dt₃, f₀, g₀, W₃)
    for n in 1:N₃
        t₃ = n * get(dt₃)
        set!(g′ₙ, π*t₃/2 + sqrt(t₃))
        VolterraIntegralEquations.advance!(VIE₃, g′ₙ, c)
    end

    # Create the arrays for plotting: 300 steps.
    x₃ = zeros(Float64, N₃+1)
    for n in 1:N₃+1
        x₃[n] = (n - 1) * get(dt₃)
    end
    y₃ = zeros(Float64, N₃+1)
    for n in 1:N₃+1
        soln = VIE₃.f[n]
        y₃[n] = get(soln)
    end
    z₃ = zeros(Float64, N₃+1)
    for n in 1:N₃+1
        z₃[n] = (2/3)*x₃[n]^(3/2)
    end
    e₃ = zeros(Float64, N₃)
    ϵ₃ = zeros(Float64, N₃)
    for n in 1:N₃
        e₃[n] = x₃[n+1]
        if abs(z₃[n+1]) < 1.0
            ϵ₃[n] = abs(z₃[n+1] - y₃[n+1])
        else
            ϵ₃[n] = abs(z₃[n+1] - y₃[n+1]) / abs(z₃[n+1])
        end
    end

    # Solve an Abel integral equation: 1000 steps.
    println("Solving the Abel integral in 1000 steps.")
    N₄ = 1000
    dt₄ = t / N₄
    W₄ = VolterraIntegralEquations.normalizedQuadratureWeights("CGS", N₄, dt₄, AbelKernel, p)
    VIE₄ = VolterraIntegralEquations.VolterraIntegralScalarEquation("CGS", N₄, dt₄, f₀, g₀, W₄)
    for n in 1:N₄
        t₄ = n * get(dt₄)
        set!(g′ₙ, π*t₄/2 + sqrt(t₄))
        VolterraIntegralEquations.advance!(VIE₄, g′ₙ, c)
    end

    # Create the arrays for plotting: 1000 steps.
    x₄ = zeros(Float64, N₄+1)
    for n in 1:N₄+1
        x₄[n] = (n - 1) * get(dt₄)
    end
    y₄ = zeros(Float64, N₄+1)
    for n in 1:N₄+1
        soln = VIE₄.f[n]
        y₄[n] = get(soln)
    end
    z₄ = zeros(Float64, N₄+1)
    for n in 1:N₄+1
        z₄[n] = (2/3)*x₄[n]^(3/2)
    end
    e₄ = zeros(Float64, N₄)
    ϵ₄ = zeros(Float64, N₄)
    for n in 1:N₄
        e₄[n] = x₄[n+1]
        if abs(z₄[n+1]) < 1.0
            ϵ₄[n] = abs(z₄[n+1] - y₄[n+1])
        else
            ϵ₄[n] = abs(z₄[n+1] - y₄[n+1]) / abs(z₄[n+1])
        end
    end

    println("Constructing the figures.")
    fig1 = Figure(; size = (809, 500)) # (500ϕ, 500), ϕ is golden ratio
    ax1 = Axis(fig1[1, 1];
        xlabel = "Upper Limit of Integration, x",
        ylabel = "Logarithm of Solution Error, log(ϵ)",
        titlesize = 24,
        xlabelsize = 20,
        ylabelsize = 20,
        yscale = log10)
    lines!(ax1, e₁, ϵ₁;
        linewidth = 3,
        linestyle = :solid,
        color = :red,
        label = "30")
    lines!(ax1, e₂, ϵ₂;
        linewidth = 3,
        linestyle = :solid,
        color = :blue,
        label = "100")
    lines!(ax1, e₃, ϵ₃;
        linewidth = 3,
        linestyle = :solid,
        color = :orange,
        label = "300")
    lines!(ax1, e₄, ϵ₄;
        linewidth = 3,
        linestyle = :solid,
        color = :cyan,
        label = "1000")
    axislegend("N =",
        position = :rb)
    fileName = string("Abel_VIE_error.png")
    save(string(myDirPath, fileName), fig1)

    fig2 = Figure(; size = (809, 500)) # (500ϕ, 500), ϕ is golden ratio
    ax2 = Axis(fig2[1, 1];
        xlabel = "Upper Limit of Integration, x",
        ylabel = "Solution, ϕ(x)",
        titlesize = 24,
        xlabelsize = 20,
        ylabelsize = 20)
    lines!(ax2, x₁, y₁;
        linewidth = 3,
        linestyle = :solid,
        color = :red,
        label = "30")
    lines!(ax2, x₂, y₂;
        linewidth = 3,
        linestyle = :solid,
        color = :blue,
        label = "100")
    lines!(ax2, x₃, y₃;
        linewidth = 3,
        linestyle = :solid,
        color = :orange,
        label = "300")
    lines!(ax2, x₄, y₄;
        linewidth = 3,
        linestyle = :solid,
        color = :cyan,
        label = "1000")
    lines!(ax2, x₂, z₂;
        linewidth = 3,
        linestyle = :solid,
        color = :black,
        label = "exact")
    axislegend("N =",
        position = :lt)
    fileName = string("Abel_VIE_soln.png")
    save(string(myDirPath, fileName), fig2)
end # Abel

end # testVolterraIntEqs