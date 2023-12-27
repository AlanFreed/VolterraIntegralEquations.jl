module testVolterraIntEqs

using
    CairoMakie,       # Pixel based figure construction.
    PhysicalFields,
    ..VolterraIntegralEquations

export
    Abel,
    memoryFns,
    persistence

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
    (k, tau) = BOX("CGS", time, (τ_σ, τ_ϵ))
    arrayBOX[1] = k
    (k, tau) = MPL("CGS", time, (α, τ_ϵ))
    arrayMPL[1] = k
    (k, tau) = RFS("CGS", time, (α, δ, τ_ϵ))
    arrayRFS[1] = k
    (k, tau) = SLS("CGS", time, (τ_ϵ,))
    arraySLS[1] = k
    for n in 1:N
        time = time + dTime
        # weakly singular kernels
        (k, tau) = CCM("CGS", time, (α, τ_ϵ))
        arrayCCM[n] = k
        (k, tau) = FLS("CGS", time, (α, τ_ϵ))
        arrayFLS[n] = k
        (k, tau) = KWW("CGS", time, (α, τ_ϵ))
        arrayKWW[n] = k
        # non-singular kernels
        (k, tau) = BOX("CGS", time, (τ_σ, τ_ϵ))
        arrayBOX[n] = k
        (k, tau) = MCM("CGS", time, (c₁, c₂, c₃, c₄, τ₁, τ₂, τ₃, τ₄))
        arrayMCM[n] = k
        (k, tau) = MPL("CGS", time, (α, τ_ϵ))
        arrayMPL[n+1] = k
        (k, tau) = RFS("CGS", time, (α, δ, τ_ϵ))
        arrayRFS[n+1] = k
        (k, tau) = SLS("CGS", time, (τ_ϵ,))
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

    fig = Figure(resolution = (809, 500)) # (500ϕ, 500), ϕ is golden ratio
    ax = Axis(fig[1, 1];
        xlabel = "time ÷ characteristic time  (t / τ_ϵ)",
        ylabel = "characteristic time × memory function  (τ_ϵ × k)",
        title = "Non-singular Memory Functions for Elastin",
        titlesize = 24,
        xlabelsize = 20,
        ylabelsize = 20)
    lines!(ax, t₂, box;
        linewidth = 3,
        linestyle = :solid,
        color = :green,
        label = "BOX")
    lines!(ax, t₁, mcm;
        linewidth = 3,
        linestyle = :solid,
        color = :cyan,
        label = "MCM")
    lines!(ax, t₂, mpl;
        linewidth = 3,
        linestyle = :solid,
        color = :blue,
        label = "MPL")
    lines!(ax, t₂, rfs;
        linewidth = 3,
        linestyle = :solid,
        color = :red,
        label = "RFS")
    lines!(ax, t₂, sls;
        linewidth = 3,
        linestyle = :solid,
        color = :black,
        label = "SLS")
    axislegend("Function",
        position = :rt)
    save(string(myDirPath, "memoryFnNonSingular.png"), fig)

    # Plot the weakly singular kernels.

    fig = Figure(resolution = (809, 500)) # (500ϕ, 500), ϕ is golden ratio
    ax = Axis(fig[1, 1];
        xlabel = "time ÷ characteristic time  (t / τ_ϵ)",
        ylabel = "characteristic time × memory function  (τ_ϵ × k)",
        title = "Weakly-singular Memory Functions for Elastin",
        titlesize = 24,
        xlabelsize = 20,
        ylabelsize = 20)
    lines!(ax, t₁, ccm;
        linewidth = 3,
        linestyle = :solid,
        color = :blue,
        label = "CCM")
    lines!(ax, t₁, fls;
        linewidth = 3,
        linestyle = :solid,
        color = :red,
        label = "FLS")
    lines!(ax, t₁, kww;
        linewidth = 3,
        linestyle = :solid,
        color = :black,
        label = "KWW")
    axislegend("Function",
        position = :rt)
    save(string(myDirPath, "memoryFnWeaklySingular.png"), fig)
end # memoryFns

function AbelKernel(systemOfUnits::String, time::PhysicalScalar, parameters::Tuple)::Tuple
    kernel = 1 / sqrt(time)
    tau = PhysicalScalar(1.0, CGS_DIMENSIONLESS)
    return (kernel, tau)
end # AbelKernel

function Abel(myDirPath::String)
    CairoMakie.activate!(type = "png")

    # Solve an Abel integral equation: solution parameters.
    t  = PhysicalScalar(6.082201995573399, CGS_DIMENSIONLESS) # upper limit of integration
    p  = ()
    c  = PhysicalScalar(1.0, CGS_DIMENSIONLESS)
    f₀ = PhysicalScalar(CGS_DIMENSIONLESS)
    g₀ = PhysicalScalar(CGS_DIMENSIONLESS)

    # Solve an Abel integral equation: 10 steps.
    N₁ = 10
    Nₘₐₓ = 1000
    sigFigs = 3
    dt₁ = t / N₁
    W₁ = normalizedQuadratureWeights("SI", dt₁, p, AbelKernel, Nₘₐₓ, sigFigs)
    g′ₙ = ArrayOfPhysicalScalars(3, CGS_DIMENSIONLESS)
    VIE₁ = VolterraIntegralScalarEquation("SI", N₁, dt₁, f₀, g₀, W₁)
    for n in 1:N₁
        tₙ₁ = (1/6 + n - 1)*dt₁
        tₙ₂ = (1/2 + n - 1)*dt₁
        tₙ₃ = (5/6 + n - 1)*dt₁
        g′ₙ[1] = π*tₙ₁/2 + sqrt(tₙ₁)
        g′ₙ[2] = π*tₙ₂/2 + sqrt(tₙ₂)
        g′ₙ[3] = π*tₙ₃/2 + sqrt(tₙ₃)
        advance!(VIE₁, g′ₙ, c)
    end

    # Create the arrays for plotting: 10 steps.
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
    N₂ = 100
    dt₂ = t / N₂
    W₂ = normalizedQuadratureWeights("SI", dt₂, (), AbelKernel, Nₘₐₓ, sigFigs)
    VIE₂ = VolterraIntegralScalarEquation("SI", N₂, dt₂, f₀, g₀, W₂)
    for n in 1:N₂
        tₙ₁ = (1/6 + n - 1)*dt₂
        tₙ₂ = (1/2 + n - 1)*dt₂
        tₙ₃ = (5/6 + n - 1)*dt₂
        g′ₙ[1] = π*tₙ₁/2 + sqrt(tₙ₁)
        g′ₙ[2] = π*tₙ₂/2 + sqrt(tₙ₂)
        g′ₙ[3] = π*tₙ₃/2 + sqrt(tₙ₃)
        advance!(VIE₂, g′ₙ, c)
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

    # Solve an Abel integral equation: 1000 steps.
    N₃ = 1000
    dt₃ = t / N₃
    W₃ = normalizedQuadratureWeights("SI", dt₃, (), AbelKernel, Nₘₐₓ, sigFigs)
    VIE₃ = VolterraIntegralScalarEquation("SI", N₃, dt₃, f₀, g₀, W₃)
    for n in 1:N₃
        tₙ₁ = (1/6 + n - 1)*dt₃
        tₙ₂ = (1/2 + n - 1)*dt₃
        tₙ₃ = (5/6 + n - 1)*dt₃
        g′ₙ[1] = π*tₙ₁/2 + sqrt(tₙ₁)
        g′ₙ[2] = π*tₙ₂/2 + sqrt(tₙ₂)
        g′ₙ[3] = π*tₙ₃/2 + sqrt(tₙ₃)
        advance!(VIE₃, g′ₙ, c)
    end

    # Create the arrays for plotting: 1000 steps.
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

    fig = Figure(resolution = (809, 500)) # (500ϕ, 500), ϕ is golden ratio
    ax = Axis(fig[1, 1];
        xlabel = "Upper Limit of Integration, 𝑥",
        ylabel = "Logarithm of Solution Error, ϵ",
        title = "Accuracy of Young's Algorithm: An Abel Kernel",
        titlesize = 24,
        xlabelsize = 20,
        ylabelsize = 20,
        yscale = log10)
    lines!(ax, e₁, ϵ₁;
        linewidth = 3,
        linestyle = :solid,
        color = :red,
        label = "10")
    lines!(ax, e₂, ϵ₂;
        linewidth = 3,
        linestyle = :solid,
        color = :blue,
        label = "100")
    lines!(ax, e₃, ϵ₃;
        linewidth = 3,
        linestyle = :solid,
        color = :black,
        label = "1000")
    axislegend("N =",
        position = :rc)
    save(string(myDirPath, "AbelVIEerror.png"), fig)

    fig = Figure(resolution = (809, 500)) # (500ϕ, 500), ϕ is golden ratio
    ax = Axis(fig[1, 1];
        xlabel = "Upper Limit of Integration, 𝑥",
        ylabel = "Solution, 𝑓(𝑥)",
        title = "Accuracy of Young's Algorithm: An Abel Kernel",
        titlesize = 24,
        xlabelsize = 20,
        ylabelsize = 20)
    lines!(ax, x₁, y₁;
        linewidth = 3,
        linestyle = :solid,
        color = :red,
        label = "10")
    lines!(ax, x₂, y₂;
        linewidth = 3,
        linestyle = :solid,
        color = :blue,
        label = "100")
    lines!(ax, x₃, y₃;
        linewidth = 3,
        linestyle = :solid,
        color = :cyan,
        label = "1000")
    lines!(ax, x₂, z₂;
        linewidth = 3,
        linestyle = :solid,
        color = :black,
        label = "exact")
    axislegend("N =",
        position = :lt)
    save(string(myDirPath, "AbelVIEsoltion.png"), fig)
end # Abel

function persistence(myDirPath::String)
    CairoMakie.activate!(type = "png")

    # Solve an Abel integral equation: solution parameters.
    t  = PhysicalScalar(6.082201995573399, CGS_DIMENSIONLESS) # upper limit of integration
    c  = PhysicalScalar(1.0, CGS_DIMENSIONLESS)
    f₀ = PhysicalScalar(CGS_DIMENSIONLESS)
    g₀ = PhysicalScalar(CGS_DIMENSIONLESS)

    # Solve an Abel integral equation: 10 steps.
    p = ()
    N₁ = 10
    Nₘₐₓ = 1000
    dt₁ = t / N₁
    sigFigs = 2
    W₁ = normalizedQuadratureWeights("CGS", dt₁, p, AbelKernel, Nₘₐₓ, sigFigs)
    g′ₙ = ArrayOfPhysicalScalars(3, CGS_DIMENSIONLESS)
    VIE₁ = VolterraIntegralScalarEquation("CGS", N₁, dt₁, f₀, g₀, W₁)
    for n in 1:N₁
        tₙ₁ = (1/6 + n - 1)*dt₁
        tₙ₂ = (1/2 + n - 1)*dt₁
        tₙ₃ = (5/6 + n - 1)*dt₁
        g′ₙ[1] = π*tₙ₁/2 + sqrt(tₙ₁)
        g′ₙ[2] = π*tₙ₂/2 + sqrt(tₙ₂)
        g′ₙ[3] = π*tₙ₃/2 + sqrt(tₙ₃)
        advance!(VIE₁, g′ₙ, c)
    end

    # Save this solution to file.
    my_dir_path = myDirPath
    json_stream = openJSONWriter(my_dir_path, "testAbelIntEq.json")
    VolterraIntegralEquations.toFile(VIE₁, json_stream)
    close(json_stream)
    json_stream = openJSONReader(my_dir_path, "testAbelIntEq.json")
    VIE₂ = VolterraIntegralEquations.fromFile(VolterraIntegralScalarEquation, json_stream)
    close(json_stream)

    # Create arrays for plotting.
    x₁ = zeros(Float64, N₁+1)
    for n in 1:N₁+1
        x₁[n] = (n - 1) * get(dt₁)
    end
    y₁ = zeros(Float64, N₁+1)
    for n in 1:N₁+1
        soln = VIE₁.f[n]
        y₁[n] = get(soln)
    end
    y₂ = zeros(Float64, N₁+1)
    for n in 1:N₁+1
        soln = VIE₂.f[n]
        y₂[n] = get(soln)
    end

    fig = Figure(resolution = (809, 500)) # (500ϕ, 500), ϕ is golden ratio
    ax = Axis(fig[1, 1];
        xlabel = "Upper Limit of Integration, 𝑥",
        ylabel = "Solution, 𝑓(𝑥)",
        title = "Verify Stored JSON File",
        titlesize = 24,
        xlabelsize = 20,
        ylabelsize = 20)
    lines!(ax, x₁, y₁;
        linewidth = 3,
        linestyle = :solid,
        color = :red,
        label = "To File")
    lines!(ax, x₁, y₂;
        linewidth = 3,
        linestyle = :dash,
        color = :black,
        label = "From File")
    axislegend("JSON",
        position = :lt)
    save(string(my_dir_path, "VIEsoltionJSON.png"), fig)
end # persistence

end # testVolterraIntEqs