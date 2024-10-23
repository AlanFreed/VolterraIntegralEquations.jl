module testVolterraIntEqs

using
    CairoMakie,
    PhysicalFields,
    VolterraIntegralEquations
 
import
    CairoMakie as CM,
    PhysicalFields as PF,
    VolterraIntegralEquations as VIE

export
    Abel,
    memoryFns

function memoryFns(my_dir_path::String)
    N = 150

    # viscoelastic material constants for elastin
    ν   = PF.PhysicalScalar(0.59, PF.CGS_DIMENSIONLESS)
    τ_ε = PF.PhysicalScalar(9.0E-4, PF.CGS_SECOND)
    τ_σ = PF.PhysicalScalar(7.0E-8, PF.CGS_SECOND)
    δ   = 0.1219 * τ_ε
    E_∞ = PF.PhysicalScalar(2.3E6, PF.BARYE)
    E_0 = (τ_ε/τ_σ)^get(ν) * E_∞

    # material constants for our model
    α = ν
    τ = τ_ε
    E = E_0
    c = (E_0 - E_∞) / E_∞
    println("The following viscoelastic model parameters are considered,")
    println("which are those of elastin fit to the FLS memory function.")
    println("  α  = ", PF.toString(α))
    println("  τ  = ", PF.toString(τ))
    println("  δ  = ", PF.toString(δ))
    println("  E₀ = ", PF.toString(E_0))
    println("  E∞ = ", PF.toString(E_∞))
    println("  c  = ", PF.toString(c), " = (E₀-E∞)/E∞\n")

    # assumed material constants for the Maxwell Chain Model.
    c₁ = PF.PhysicalScalar(0.25, PF.CGS_DIMENSIONLESS)
    c₂ = PF.PhysicalScalar(0.25, PF.CGS_DIMENSIONLESS)
    c₃ = PF.PhysicalScalar(0.25, PF.CGS_DIMENSIONLESS)
    c₄ = PF.PhysicalScalar(0.25, PF.CGS_DIMENSIONLESS)
    τ₁ = τ_σ
    τ  = exp10(log10(get(τ_σ)) + (log10(get(τ_ε)) - log10(get(τ_σ)))/3)
    τ₂ = PF.PhysicalScalar(τ, PF.CGS_SECOND)
    τ  = exp10(log10(get(τ_σ)) + 2(log10(get(τ_ε)) - log10(get(τ_σ)))/3)
    τ₃ = PF.PhysicalScalar(τ, PF.CGS_SECOND)
    τ₄ = τ_ε
    println("For the Prony series, let c_i = 0.25, i = 1,...,4, and")
    println("   τ₁ = ", PF.toString(τ₁))
    println("   τ₂ = ", PF.toString(τ₂))
    println("   τ₃ = ", PF.toString(τ₃))
    println("   τ₄ = ", PF.toString(τ₄))
    println()

    # Create the arrays to hold values for the memory function.
    reciprocalTime = PF.PhysicalUnits("CGS", 0, 0, 0, -1, 0, 0, 0)
    # weakly singular kernels
    arrayCCM = PF.ArrayOfPhysicalScalars(N, reciprocalTime)
    arrayFLS = PF.ArrayOfPhysicalScalars(N, reciprocalTime)
    arrayKWW = PF.ArrayOfPhysicalScalars(N, reciprocalTime)
    # non-singular kernels
    arrayBOX = PF.ArrayOfPhysicalScalars(N+1, reciprocalTime)
    arrayMCM = PF.ArrayOfPhysicalScalars(N, reciprocalTime)
    arrayMPL = PF.ArrayOfPhysicalScalars(N+1, reciprocalTime)
    arrayRFS = PF.ArrayOfPhysicalScalars(N+1, reciprocalTime)
    arraySLS = PF.ArrayOfPhysicalScalars(N+1, reciprocalTime)

    # Populate the memory function arrays
    time  = PF.PhysicalScalar(PF.CGS_SECOND)
    dtime = 3 * τ_ε / N
    (name, k, tau) = VIE.BOX("CGS", time, (τ_σ, τ_ε))
    arrayBOX[1] = k
    (name, k, tau) = VIE.MPL("CGS", time, (α, τ_ε))
    arrayMPL[1] = k
    (name, k, tau) = VIE.RFS("CGS", time, (α, δ, τ_ε))
    arrayRFS[1] = k
    (name, k, tau) = VIE.SLS("CGS", time, (τ_ε,))
    arraySLS[1] = k
    for n in 1:N
        time = time + dtime
        # weakly singular kernels
        (name, k, tau) = VIE.CCM("CGS", time, (α, τ_ε))
        arrayCCM[n] = k
        (name, k, tau) = VIE.FLS("CGS", time, (α, τ_ε))
        arrayFLS[n] = k
        (name, k, tau) = VIE.KWW("CGS", time, (α, τ_ε))
        arrayKWW[n] = k
        # non-singular kernels
        (name, k, tau) = VIE.BOX("CGS", time, (τ_σ, τ_ε))
        arrayBOX[n] = k
        (name, k, tau) = VIE.MCM("CGS", time, (c₁, c₂, c₃, c₄, τ₁, τ₂, τ₃, τ₄))
        arrayMCM[n] = k
        (name, k, tau) = VIE.MPL("CGS", time, (α, τ_ε))
        arrayMPL[n+1] = k
        (name, k, tau) = VIE.RFS("CGS", time, (α, δ, τ_ε))
        arrayRFS[n+1] = k
        (name, k, tau) = VIE.SLS("CGS", time, (τ_ε,))
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
    PF.set!(time, 0.0)
    t₂[1]  = get(time) / get(τ_ε)
    box[1] = get(τ_ε) * get(arrayBOX[1])
    mcm[1] = get(τ_ε) * get(arrayMCM[1])
    mpl[1] = get(τ_ε) * get(arrayMPL[1])
    rfs[1] = get(τ_ε) * get(arrayRFS[1])
    sls[1] = get(τ_ε) * get(arraySLS[1])
    for n in 1:N
        # weakly singular kernels
        time   = time + dtime
        t₁[n]  = get(time) / get(τ_ε)
        ccm[n] = get(τ_ε) * get(arrayCCM[n])
        fls[n] = get(τ_ε) * get(arrayFLS[n])
        kww[n] = get(τ_ε) * get(arrayKWW[n])
        # non-singular kernels
        t₂[n+1]  = t₁[n]
        box[n]   = get(τ_ε) * get(arrayBOX[n+1])
        mcm[n]   = get(τ_ε) * get(arrayMCM[n])
        mpl[n+1] = get(τ_ε) * get(arrayMPL[n+1])
        rfs[n+1] = get(τ_ε) * get(arrayRFS[n+1])
        sls[n+1] = get(τ_ε) * get(arraySLS[n+1])
    end

    # Plot the non-singular kernels.

    fig1 = CM.Figure(size = (809, 500)) # (500ϕ, 500), ϕ is golden ratio
    ax1 = CM.Axis(fig1[1, 1];
        xlabel = "time ÷ characteristic time  (t / τ_ε)",
        ylabel = "characteristic time × memory function  (τ_ε × k)",
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
    save(string(my_dir_path, "files/memoryFnNonSingular.png"), fig1)

    # Plot the weakly singular kernels.

    fig2 = CM.Figure(size = (809, 500)) # (500ϕ, 500), ϕ is golden ratio
    ax2 = CM.Axis(fig2[1, 1];
        xlabel = "time ÷ characteristic time  (t / τ_ε)",
        ylabel = "characteristic time × memory function  (τ_ε × k)",
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
    save(string(my_dir_path, "files/memoryFnWeaklySingular.png"), fig2)
end # memoryFns

function AbelKernel(system_of_units::String, time::PF.PhysicalScalar, parameters::Tuple)::Tuple
    kernel = 1 / sqrt(time)
    tau = PF.PhysicalScalar(1.0, PF.CGS_DIMENSIONLESS)
    return ("Abel", kernel, tau)
end # AbelKernel

function Abel(my_dir_path::String)
    CairoMakie.activate!(type = "png")

    # Solve an Abel integral equation: solution parameters.
    t = PF.PhysicalScalar(10.0, PF.CGS_DIMENSIONLESS) # upper limit of integration
    p = ()
    c = PF.PhysicalScalar(1.0, PF.CGS_DIMENSIONLESS)
    f₀ = PF.PhysicalScalar(PF.CGS_DIMENSIONLESS)
    g₀ = PF.PhysicalScalar(PF.CGS_DIMENSIONLESS)

    # Solve an Abel integral equation: 10 steps.
    println("Solving the Abel integral in 10 steps.")
    N₀ = 10
    dt₀ = t / N₀
    W₀ = VIE.normalizedQuadratureWeights("CGS", N₀, dt₀, AbelKernel, p)
    g′ₙ = PF.PhysicalScalar(PF.CGS_DIMENSIONLESS)
    VIE₀ = VIE.VolterraIntegralScalarEquation("CGS", N₀, dt₀, f₀, g₀, W₀)
    for n in 1:N₀
        tₙ = n * get(dt₀)
        PF.set!(g′ₙ, π*tₙ/2 + sqrt(tₙ))
        VIE.advance!(VIE₀, g′ₙ, c)
    end

    # Create the arrays for plotting: 10 steps.
    x₀ = zeros(Float64, N₀+1)
    for n in 1:N₀+1
        x₀[n] = (n - 1) * get(dt₀)
    end
    y₀ = zeros(Float64, N₀+1)
    for n in 1:N₀+1
        soln = VIE₀.f[n]
        y₀[n] = get(soln)
    end
    z₀ = zeros(Float64, N₀+1)
    for n in 1:N₀+1
        z₀[n] = (2/3)*x₀[n]^(3/2)
    end
    e₀ = zeros(Float64, N₀)
    ε₀ = zeros(Float64, N₀)
    for n in 1:N₀
        e₀[n] = x₀[n+1]
        if abs(z₀[n+1]) < 1.0
            ε₀[n] = abs(z₀[n+1] - y₀[n+1])
        else
            ε₀[n] = abs(z₀[n+1] - y₀[n+1]) / abs(z₀[n+1])
        end
    end

    # Solve an Abel integral equation: 30 steps.
    println("Solving the Abel integral in 30 steps.")
    N₁ = 30
    dt₁ = t / N₁
    W₁ = VIE.normalizedQuadratureWeights("CGS", N₁, dt₁, AbelKernel, p)
    g′ₙ = PF.PhysicalScalar(PF.CGS_DIMENSIONLESS)
    VIE₁ = VIE.VolterraIntegralScalarEquation("CGS", N₁, dt₁, f₀, g₀, W₁)
    for n in 1:N₁
        tₙ = n * get(dt₁)
        PF.set!(g′ₙ, π*tₙ/2 + sqrt(tₙ))
        VIE.advance!(VIE₁, g′ₙ, c)
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
    ε₁ = zeros(Float64, N₁)
    for n in 1:N₁
        e₁[n] = x₁[n+1]
        if abs(z₁[n+1]) < 1.0
            ε₁[n] = abs(z₁[n+1] - y₁[n+1])
        else
            ε₁[n] = abs(z₁[n+1] - y₁[n+1]) / abs(z₁[n+1])
        end
    end

    # Solve an Abel integral equation: 100 steps.
    println("Solving the Abel integral in 100 steps.")
    N₂ = 100
    dt₂ = t / N₂
    W₂ = VIE.normalizedQuadratureWeights("CGS", N₂, dt₂, AbelKernel, p)
    VIE₂ = VIE.VolterraIntegralScalarEquation("CGS", N₂, dt₂, f₀, g₀, W₂)
    for n in 1:N₂
        t₂ = n * get(dt₂)
        PF.set!(g′ₙ, π*t₂/2 + sqrt(t₂))
        VIE.advance!(VIE₂, g′ₙ, c)
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
    ε₂ = zeros(Float64, N₂)
    for n in 1:N₂
        e₂[n] = x₂[n+1]
        if abs(z₂[n+1]) < 1.0
            ε₂[n] = abs(z₂[n+1] - y₂[n+1])
        else
            ε₂[n] = abs(z₂[n+1] - y₂[n+1]) / abs(z₂[n+1])
        end
    end

    # Solve an Abel integral equation: 300 steps.
    println("Solving the Abel integral in 300 steps.")
    N₃ = 300
    dt₃ = t / N₃
    W₃ = VIE.normalizedQuadratureWeights("CGS", N₃, dt₃, AbelKernel, p)
    VIE₃ = VIE.VolterraIntegralScalarEquation("CGS", N₃, dt₃, f₀, g₀, W₃)
    for n in 1:N₃
        t₃ = n * get(dt₃)
        PF.set!(g′ₙ, π*t₃/2 + sqrt(t₃))
        VIE.advance!(VIE₃, g′ₙ, c)
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
    ε₃ = zeros(Float64, N₃)
    for n in 1:N₃
        e₃[n] = x₃[n+1]
        if abs(z₃[n+1]) < 1.0
            ε₃[n] = abs(z₃[n+1] - y₃[n+1])
        else
            ε₃[n] = abs(z₃[n+1] - y₃[n+1]) / abs(z₃[n+1])
        end
    end

    # Solve an Abel integral equation: 1000 steps.
    println("Solving the Abel integral in 1000 steps.")
    N₄ = 1000
    dt₄ = t / N₄
    W₄ = VIE.normalizedQuadratureWeights("CGS", N₄, dt₄, AbelKernel, p)
    VIE₄ = VIE.VolterraIntegralScalarEquation("CGS", N₄, dt₄, f₀, g₀, W₄)
    for n in 1:N₄
        t₄ = n * get(dt₄)
        PF.set!(g′ₙ, π*t₄/2 + sqrt(t₄))
        VIE.advance!(VIE₄, g′ₙ, c)
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
    ε₄ = zeros(Float64, N₄)
    for n in 1:N₄
        e₄[n] = x₄[n+1]
        if abs(z₄[n+1]) < 1.0
            ε₄[n] = abs(z₄[n+1] - y₄[n+1])
        else
            ε₄[n] = abs(z₄[n+1] - y₄[n+1]) / abs(z₄[n+1])
        end
    end

    println("Constructing the figures.")
    fig1 = CM.Figure(size = (809, 500)) # (500ϕ, 500), ϕ is golden ratio
    ax1 = CM.Axis(fig1[1, 1];
        xlabel = "Upper Limit of Integration, x",
        ylabel = "Logarithm of Solution Error, log(ε)",
        titlesize = 24,
        xlabelsize = 20,
        ylabelsize = 20,
        yscale = log10)
    lines!(ax1, e₀, ε₀;
        linewidth = 3,
        linestyle = :solid,
        color = :red,
        label = "10")
    lines!(ax1, e₁, ε₁;
        linewidth = 3,
        linestyle = :solid,
        color = :magenta,
        label = "30")
    lines!(ax1, e₂, ε₂;
        linewidth = 3,
        linestyle = :solid,
        color = :orange,
        label = "100")
    lines!(ax1, e₃, ε₃;
        linewidth = 3,
        linestyle = :solid,
        color = :cyan,
        label = "300")
    lines!(ax1, e₄, ε₄;
        linewidth = 3,
        linestyle = :solid,
        color = :blue,
        label = "1000")
    axislegend("N =",
        position = :rb)
    file_name = string("files/Abel_VIE_error.png")
    save(string(my_dir_path, file_name), fig1)

    fig2 = CM.Figure(size = (809, 500)) # (500ϕ, 500), ϕ is golden ratio
    ax2 = CM.Axis(fig2[1, 1];
        xlabel = "Upper Limit of Integration, x",
        ylabel = "Solution, ϕ(x)",
        titlesize = 24,
        xlabelsize = 20,
        ylabelsize = 20)
    lines!(ax2, x₀, y₀;
        linewidth = 1,
        linestyle = :solid,
        color = :red,
        label = "10")
    lines!(ax2, x₁, y₁;
        linewidth = 1,
        linestyle = :solid,
        color = :magenta,
        label = "30")
    lines!(ax2, x₂, y₂;
        linewidth = 1,
        linestyle = :solid,
        color = :orange,
        label = "100")
    lines!(ax2, x₃, y₃;
        linewidth = 1,
        linestyle = :solid,
        color = :cyan,
        label = "300")
    lines!(ax2, x₄, y₄;
        linewidth = 1,
        linestyle = :solid,
        color = :blue,
        label = "1000")
    lines!(ax2, x₂, z₂;
        linewidth = 1,
        linestyle = :solid,
        color = :black,
        label = "exact")
    axislegend("N =",
        position = :lt)
    file_name = string("files/Abel_VIE_soln.png")
    save(string(my_dir_path, file_name), fig2)
end # Abel

end # testVolterraIntEqs
