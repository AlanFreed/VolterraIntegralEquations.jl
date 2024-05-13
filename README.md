# VolterraIntegralEquations

This module provides a solver for Volterra integral equations of the second kind [13]. Specifically, this solver solves an integro-differential equation of the form
$$
\boldsymbol{f}^{\prime}(t) = \boldsymbol{g}^{\prime}(t) - c(t) \int_0^t K(t-τ) \, \boldsymbol{f}^{\prime}(τ) \, \mathrm{d}τ
$$
where the lower limit of integration is an initial time, taken to be 0, while the upper limit of integration is some current time $t$, with $\tau$ denoting a dummy variable of integration. For illustrative purposes, the independent variable in this document is taken to be time, but this is not necessary. 

This equation becomes a Volterra integral equation of the first kind whenever $\boldsymbol{g}^{\prime}$ is $\boldsymbol{0}$.

This Volterra integral equation describes an ordinary differential equation that, in turn, must be solved, i.e., it is an integro-differential equation. Here derivative $\boldsymbol{g}^{\prime} = \mathrm{d}\boldsymbol{g}/\mathrm{d}t$ is a known rate for some control function $\boldsymbol{g}$, while derivative $\boldsymbol{f}^{\prime} = \mathrm{d}\boldsymbol{f}/\mathrm{d}t$ is an unknown rate for a response function $\boldsymbol{f}$ that is to be determined, whose evolution is characterized by a Volterra integral equation. Coefficient $c(t)$ and kernel $K(t-τ)$ are also taken to be known quantities.

Functions $\boldsymbol{f}$ and $\boldsymbol{g}$ have the same physical units. They may be scalars, vectors or tensors in stature. Coefficient $c$ is a dimensionless, scalar-valued function. Kernel $K$ is a memory function, which is the derivative of a generalized creep function. It is a scalar-valued function with physical units of $t^{-1}$, e.g., reciprocal time whenever $t$ has units of time.

To use this module, you will need to add the following repositories to your project:

```
using Pkg
Pkg.add(url = "https://github.com/AlanFreed/PhysicalFields.jl")
Pkg.add(url = "https://github.com/AlanFreed/VolterraIntegralEquations.jl")
```

This software was written to accompany a book being written by the authors of this software [7].

## Example: A Viscoelastic Fiber

A linear viscoelastic fiber can be modeled as a Volterra integral equation of the second kind, viz., [6]
$$
\frac{\mathrm{d} \sigma}{\mathrm{d} t} = \mathcal{E}(t) \left(
	\frac{\mathrm{d} \epsilon}{\mathrm{d} t} -
	\frac{E_0 - E_{\infty}}{E_0 \, E_{\infty}}
	\int_0^t K (t - \tau) \, \frac{\mathrm{d} \sigma}{\mathrm{d}\tau} \,
	\mathrm{d}\tau \right)
$$
where $\sigma$ is stress and $\epsilon$ is strain, and where $\mathcal{E}$,  $E_{\infty}$, and $E_0$ are its tangent, rubbery, and glassy elastic moduli, respectively, the former being a possible function of time, while the latter two are material constants. Kernel $K$ is a positive, monotonic-decreasing function known as a memory function whose units are reciprocal time.  The first term on the right-hand side of the above equation provides for a glassy elastic change in stress that is attenuated by the second term, which introduces a viscous loss to this change in stress.

The elastic tangent modulus $\mathcal{E}$ for a linear Hookean fiber is
$$
\mathcal{E} = E_0
$$
wherein $E_0$ denotes its glassy modulus.

In contrast, for a nonlinear biologic fiber, its tangent modulus associates with a nonlinear compliance of
$$
\frac{1}{\mathcal{E}} = \frac{1}{E_0} + \frac{\beta}{E_r \beta + 2(\sigma - \sigma_r)} \sqrt{\frac{\beta E_r}{\beta E_r + 2(\sigma - \sigma_r)}}
$$
wherein $\sigma_r$ denotes a residual stress, and $\beta$ designates a limiting state for an internal strain, a strain caused by molecular reconfiguration.  The elastic tangent modulus associated with a fiber's strain-free reference configuration is designated as $E_r$ ($> 0$). Within a biologic fiber's linear region of response, wherein strain is caused by molecular stretching, moduli $E_\infty$ ($> E_r$) and $E_0$ ($> E_{\infty}$) denote its rubbery and glassy elastic moduli, respectively. The resulting elastic tangent modulus $\mathcal{E}$ is thereby bounded by $E_r$ from below and $E_0$ from above.

### Volterra Functions

For a viscoelastic fiber, the unknown forcing function to be ascertained is a stress rate, viz.,
$$
f^{\prime}(t) = \frac{\mathrm{d} \sigma}{\mathrm{d} t}
$$
that is to be solved for in terms of:
1) a known kernel function and its coefficient
$$
c(t) \, K(t-\tau) 
\quad \text{with coefficient} \quad 
c(t) = \mathcal{E}(t) \, \frac{E_0 - E_{\infty}}{E_0 \, E_{\infty}}
$$
resulting in a modulus-scaled memory function, and 
2) a known control function
$$
g^{\prime}(t) = \mathcal{E}(t) \, \frac{\mathrm{d}\epsilon}{\mathrm{d}t}
$$
which is a modulus-scaled strain rate.

## Memory Functions

A selection of positive, monotonic-decreasing, viscoelastic, memory functions $K(t)$, whose units are reciprocal time, have been preprogrammed into this software. These are representative of the many kernel functions that have been proposed in the literature. These memory functions are the derivatives of creep functions, the latter of which are more commonly found in the literature. Consequently, their characteristic times associate with creep (not stress relaxation). See Freed [6,8] for a discussion of these functions.

1) **BOX**: the *𝑏𝑜𝑥* energy dissipation model of Neuber [12], a.k.a. Fung's [9] **Q**uasi-**L**inear **V**iscoelastic (QLV) kernel, has a memory function of
$$
K(t) = \frac{\exp(-t/\tau_2) - \exp(-t/\tau_1)}{t \, \ln(\tau_2/\tau_1)}
\quad \text{with} \quad
K(0) = \frac{1/\tau_1 - 1/\tau_2}{\ln(\tau_2/\tau_1)}
$$
wherein $\tau_1$ and $\tau_2$ are characteristic times, ordered so that $0 < \tau_1 < \tau_2$. It is at rates between $1/\tau_2$ and $1/\tau_1$ where dissipation is considered to occur.

2) **CCM**: **C**ole and **C**ole's [4,5] power-law **M**odel has a memory function of
$$
K(t) = \frac{\alpha}{t} \, \left( \frac{t}{\tau} \right)^{\alpha} 
\frac{1}{\bigl( 1 + (t / \tau)^{\alpha} \bigr)^2}
\quad \text{with} \quad
K(0) = \infty
$$
wherein $\tau$ is a characteristic time and $\alpha \in (0,1]$ is the exponent of a power law. The CCM memory kernel is weakly singular whenever $\alpha \in (0,1)$.

3) **FLS**: Caputo and Mainardi's [2,3] **F**ractional **L**inear **S**olid has a memory function of
$$
K(t) = -\frac{E_{\alpha,0} \left( - \left( t / \tau \right)^{\alpha} \right)}{t}
\quad \text{with} \quad
K(0) = \infty
$$
wherein $\tau$ is a characteristic time, and $\alpha \in (0, 1]$ is a fractional order of evolution, with $E_{\alpha,\beta}(t)$ being the two-parameter Mittag-Leffler function. The FLS model contains the SLS model as a special case; specifically, they become equivalent whenever $\alpha = 1$. Mainardi's memory kernel is weakly singular whenever $\alpha \in (0,1)$.

4) **KWW**: **K**ohlrausch's [10] and **W**illiams & **W**atts' [14] stretched exponential has a memory function of
$$
K(t) = \alpha \, \left( \frac{t}{\tau} \right)^{\alpha} \;
\frac{\exp \bigl( -(t/\tau)^{\alpha} \bigr)}{t}
\quad \text{with} \quad
K(0) = \infty
$$
wherein $\tau$ is a characteristic time and $\alpha \in (0,1]$ is an exponent for the power of the argument in the exponential. The KWW model contains the SLS model as a special case; specifically, they become equivalent whenever $\alpha = 1$. The KWW memory kernel is weakly singular whenever $\alpha \in (0,1)$.

5) **MCM**: **M**axwell's **C**hain **M**odel, a.k.a. the Prony series model, has a memory function of
$$
K(t) = \sum_{\ell=1}^L \frac{c_{\ell}}{\tau_{\ell}} \, \exp(-t/\tau_{\ell})
\quad \text{with} \quad
K(0) = \sum_{\ell=1}^L \frac{c_{\ell}}{\tau_{\ell}}
$$
whose coefficients $c_{\ell}$ are positive and sum as $\sum_{\ell=1}^L c_{\ell} = 1$, and whose characteristic times $\tau_{\ell}$, of which there are $L$, are ordered such that $0 < \tau_1 < \tau_2 < \cdots < \tau_L$. This is the popular Prony series that pervades the viscoelastic literature.  Its parameters are not unique, and therefore, they lack physical interpretation.

6) **MPL**: Williams' [14] **M**odified **P**ower-**L**aw model has a memory function of
$$
K(t) = \frac{\alpha}{\tau} \frac{1}{(1 + t / \tau)^{1+\alpha}}
\quad \text{with} \quad
K(0) = \frac{\alpha}{\tau}
$$
wherein $\tau$ is a characteristic time and $\alpha$ ($>0$) is in the exponent of a power law. The MPL kernel equates with the SLS kernel at time 0 whenever $\alpha = 1$.

7) **RFS**: Freed and Rajagopal's [8] **R**egularized **F**ractional linear **S**olid has a memory function of
$$
K(t) = \frac{-1}{E_{\alpha,1} \bigl( - ( \delta / \tau )^{\alpha} \bigr)} 
\frac{E_{\alpha,0} \bigl( - \bigl( ( \delta + t ) / \tau \bigr)^{\alpha} \bigr)}{\delta + t}
\quad \text{with} \quad
K(0) = - \frac{E_{\alpha,0} \bigl( - ( \delta / \tau )^{\alpha} \bigr)}{\delta E_{\alpha,1} \bigl( - ( \delta / \tau )^{\alpha} \bigr)}
$$
wherein $\tau$ is a characteristic time, $\alpha \in (0, 1)$ is a fractional order of evolution, and $\delta$ is a regularization parameter satisfying the implicit equation
$$
\frac{\delta}{\tau} = -E_{\alpha,0} \left( - \left( 
\frac{\delta}{\tau} \right)^{\alpha}\right)
$$
which effectively shifts the singularity of an FLS kenel a short distance into negative time so that the RFS creep compliance equates with the SLS creep compliance at time 0. As a consequence, the RFS memory kernel is finite valued throughout its domain of application.

8) **SLS**: Zener's [17] **S**tandard **L**inear **S**olid has a Maxwell-Debye kernel resulting in a memory function of
$$
K(t) = \frac{\mathrm{exp}(-t / \tau)}{\tau}
\quad \text{with} \quad
K(0) = \frac{1}{\tau}
$$
wherein $\tau$ is a characteristic time.

A capability is also provided in this software for the user to be able to define and use their own memory function, if it is other than one of those listed above.

## Numerical Method

This solver implements a numerical method developed in an appendix of a book that is currently being written by the authors of this software [6].

In practice, global solutions to such integral equations are gotten at sequential instants in time.  In our case, given an initial time $t_0 = 0$, solutions are sought at times $t_1, t_2, \ldots, t_N$ \mbox{($t_0 < t_1 < t_2 < \cdots < t_N$)} of which there are $N$ solutions to be gotten.  These global instants are taken to be uniformly spaced in time, separated by a time increment of $\mathrm{d} t$ such that $t_n = t_{n-1} + \mathrm{d}t$ for all $n = 1, 2, \ldots, N$.

It is sufficient for our purposes to consider quadrature nodes $t_n$ that are of adequate density so as to minimize the method's approximation error. To assess if a given nodal density is adequate in this regard, multiple solutions with different nodal densities ought to be acquired, whose errors of approximation converge with increasing nodal density.

Solutions to Volterra integrals of the second kind can be advanced through a block-by-block algorithmic strategy of
$$
\boldsymbol{f}^{\prime}(t_1) = \boldsymbol{g}^{\prime}(t_1) - c(t_1) \int_{0}^{t_1} K(t_1 - \tau) \, \boldsymbol{f}^{\prime}(\tau) \, \mathrm{d}\tau
$$
$$
\boldsymbol{f}^{\prime}(t_2) = \boldsymbol{g}^{\prime}(t_2) - c(t_2) \int_{0}^{t_2} K(t_2 - \tau) \, \boldsymbol{f}^{\prime}(\tau) \, \mathrm{d}\tau 
$$
$$
\boldsymbol{f}^{\prime}(t_3) = \boldsymbol{g}^{\prime}(t_3) - c(t_3) \int_{0}^{t_3} K(t_3 - \tau) \, \boldsymbol{f}^{\prime}(\tau) \, \mathrm{d}\tau 
$$
$$
\vdots
$$
$$
\boldsymbol{f}^{\prime}(t_N) = \boldsymbol{g}^{\prime}(t_N) - c(t_N) \int_{0}^{t_N} K(t_N - \tau) \, \boldsymbol{f}^{\prime}(\tau) \, \mathrm{d}\tau
$$
subjected to initial conditions of
$$
\boldsymbol{f}(t_0) = \boldsymbol{f}_0 
\quad \text{and} \quad
\boldsymbol{g}(t_0) = \boldsymbol{g}_0
$$
that, from the fundamental theorem of calculus, can be rewritten as
$$
\boldsymbol{f}^{\prime}(t_1) = \boldsymbol{g}^{\prime}(t_1) - c(t_1) \int_{0}^{t_1} K(t_1 - \tau) \, \boldsymbol{f}^{\prime}(\tau) \, \mathrm{d}\tau
$$
$$
\begin{aligned}
\boldsymbol{f}^{\prime}(t_2) & = \boldsymbol{g}^{\prime}(t_2) - c(t_2) \int_{0}^{t_1} K(t_2 - \tau) \, \boldsymbol{f}^{\prime}(\tau) \, \mathrm{d}\tau \\
& \phantom{= \boldsymbol{g}^{\prime}(t_3)} \; - c(t_2)
\int_{t_1}^{t_2} K(t_2 - \tau) \, \boldsymbol{f}^{\prime}(\tau) \, \mathrm{d}\tau
\end{aligned}
$$
$$
\begin{aligned}
\boldsymbol{f}^{\prime}(t_3) & = \boldsymbol{g}^{\prime}(t_3) - c(t_3) \int_{0}^{t_1} K(t_3 - \tau) \, \boldsymbol{f}^{\prime}(\tau) \, \mathrm{d}\tau \\
& \phantom{= \boldsymbol{g}^{\prime}(t_3)} \; - c(t_3)
		\int_{t_1}^{t_2} K(t_3 - \tau) \, \boldsymbol{f}^{\prime}(\tau) \, \mathrm{d}\tau \\
& \phantom{= \boldsymbol{g}^{\prime}(t_3)} \; - c(t_3)
		\int_{t_2}^{t_3} K(t_3 - \tau) \, \boldsymbol{f}^{\prime}(\tau) \, \mathrm{d}\tau
\end{aligned}
$$
$$
\vdots
$$
$$
\begin{aligned}
\boldsymbol{f}^{\prime}(t_N) & = \boldsymbol{g}^{\prime}(t_N) - c(t_N) \int_{0}^{t_1} K(t_N - \tau) \, \boldsymbol{f}^{\prime}(\tau) \, \mathrm{d}\tau \\
& \phantom{\= \boldsymbol{g}^{\prime}(t_N)} \; - c(t_N)
        \int_{t_1}^{t_2} K(t_N - \tau) \, \boldsymbol{f}^{\prime}(\tau) \, \mathrm{d}\tau \\
& \phantom{= \boldsymbol{g}^{\prime}(t_N) -} \vdots \\
& \phantom{= \boldsymbol{g}^{\prime}(t_N)} \; - c(t_N)
        \int_{t_{N-1}}^{t_N} K(t_N - \tau) \, \boldsymbol{f}^{\prime}(\tau) \, \mathrm{d}\tau
\end{aligned}
$$
wherein, e.g., the product integral $\int_{0}^{t_1} K(t_2 - \tau) \, \boldsymbol{f}^{\prime}(\tau) \, \mathrm{d}\tau$ has a kernel $K$ whose fixed time $t_2$ lies outside the interval of integration, in this case $[0, t_1]$, over which time $\tau$ spans, while the product integral $\int_{t_1}^{t_2} K(t_2 - \tau) \, \boldsymbol{f}^{\prime}(\tau) \, \mathrm{d}\tau$ has a kernel $K$ whose fixed time, in this case $t_2$, is now the upper limit of integration. The first of these two integrals is called a Fredholm integral, while the latter is called a Volterra integral.

### Quadrature Rule

The convolution integrals in the above Volterra integral equations are treated as product integrals [16], a product between a forcing function and a kernel function, that are considered to obey the following system of implicit equations
$$
\begin{aligned}
& \boldsymbol{f}^{\prime}_1 = \boldsymbol{g}^{\prime}_1 - c_1 \left( W_N \boldsymbol{f}^{\prime}_1 + \boldsymbol{\epsilon}_1 \right) \\
& \boldsymbol{f}^{\prime}_2 = \boldsymbol{g}^{\prime}_2 - c_2 \left( W_N \boldsymbol{f}^{\prime}_2 + W_{N-1} \boldsymbol{f}^{\prime}_1 + \boldsymbol{\epsilon}_2 \right) \\
& \boldsymbol{f}^{\prime}_3 = \boldsymbol{g}^{\prime}_3 - c_3 \left( W_N \boldsymbol{f}^{\prime}_3 + W_{N-1} \boldsymbol{f}^{\prime}_2 + W_{N-2} \boldsymbol{f}^{\prime}_1 + \boldsymbol{\epsilon}_3 \right) \\
& \phantom{= \boldsymbol{f}^{\prime}_3} \!\vdots \\
& \boldsymbol{f}^{\prime}_N = \boldsymbol{g}^{\prime}_N - c_N \left( W_N \boldsymbol{f}^{\prime}_N + W_{N-1} \boldsymbol{f}^{\prime}_{N-1} + \cdots + W_1 \boldsymbol{f}^{\prime}_1 + \boldsymbol{\epsilon}_N \right)
\end{aligned}
$$
where $W_n$ is the weight of quadrature at step $n$, $n = 1, 2, \ldots , N$, in the last step of integration $N$, whose approximation error is $\boldsymbol{\epsilon}_n$. 

Assuming the nodal density to be sufficiently large so that the errors of approximation $\boldsymbol{\epsilon}_n$ can be safely neglected, then upon rearranging these expressions one gets
$$
\begin{aligned}
& \boldsymbol{f}^{\prime}_1 = (1 + c_1 W_N)^{-1} \boldsymbol{g}^{\prime}_1 \\
& \boldsymbol{f}^{\prime}_2 = (1 + c_2 W_N)^{-1} \left( \boldsymbol{g}^{\prime}_2 - c_2 W_{N-1} \boldsymbol{f}^{\prime}_1 \right) \\
& \boldsymbol{f}^{\prime}_3 = (1 + c_3 W_N)^{-1} \left( \boldsymbol{g}^{\prime}_3 - c_3 \left( W_{N-1} \boldsymbol{f}^{\prime}_2 + W_{N-2} \boldsymbol{f}^{\prime}_1 \right) \right) \\
& \phantom{= \boldsymbol{f}^{\prime}_3} \!\vdots \\
& \boldsymbol{f}^{\prime}_N = (1 + c_N W_N)^{-1} \left(\boldsymbol{g}^{\prime}_N - c_N \left( W_{N-1} \boldsymbol{f}^{\prime}_{N-1} + W_{N-2} \boldsymbol{f}^{\prime}_{N-2} + \cdots + W_1 \boldsymbol{f}^{\prime}_1 \right) \right)
\end{aligned}
$$
where each integral above is taken to obey a product quadrature rule [16] of
$$
\int_{t_{n-1}}^{t_n} K (t_N - \tau) \, \boldsymbol{f}^{\prime} ( \tau ) \, \mathrm{d} \tau = W_n f^{\prime}_n + \boldsymbol{\epsilon}_n
$$
with an approximation error of $\boldsymbol{\epsilon}_n$. All functions are to be evaluated at the end point $t_n = n \, \mathrm{d}t$ of each local integration in that
$$
c_n = c(t_n) , \quad \boldsymbol{f}^{\prime}_n = \boldsymbol{f}^{\prime} (t_n) \quad \text{and} \quad \boldsymbol{g}^{\prime}_n = \boldsymbol{g}^{\prime}(t_n)
$$
which is consistent with the fact that the solution $\boldsymbol{f}^{\prime}_N$ is evaluated at the upper limit of integration. 

The weights of quadrature are described by integrals, specifically
$$
W_n = \int_{t_{n-1}}^{t_n} K ( t_N - \tau ) \, \mathrm{d} \tau
$$
that upon adopting the transformation
$$
\int_a^b \phi (x) \, \mathrm{d}x = \frac{b-a}{2} \int_{-1}^1 \phi \left( \tfrac{b-a}{2} \xi + \tfrac{a+b}{2} \right) \mathrm{d} \xi
$$
place the above integrals for quadrature into Gauss' form $\int_{-1}^1 f(x) \, \mathrm{d}x = \sum_{s=1}^S w_s f(\xi_s) + \epsilon$, thereby resulting in the following formula for our weights of quadrature
$$
W_n = \frac{\mathrm{d}t}{2} \sum_{s=1}^S w_s K \Bigl( \bigl( N - \tfrac{1}{2} ( 2n - 1 + \xi_s ) \bigr) \mathrm{d} t \Bigr)
$$
where $w_s$ and $\xi_s$ are the weights and nodes of Gauss' quadrature rule, see the table below. For our application, we choose $S=4$ whenever $t_n < \tfrac{1}{2} \tau$, $S=3$ whenever $\tfrac{1}{2} \tau \le t_n < 2\tau$, $S=2$ whenever $2\tau \le t_n , 5\tau$, and $S=1$ whenever $t_n \ge 5\tau$, where $\tau$ is the characteristic time for kernel $K$.

| $S=$ | $w_s=$ | $\xi_s =$ |
| --- | --- | --- |
| 1 | $\{ 2 \}$ | $\{ 0 \}$ |
| 2 | $\{ 1, 1 \}$ | $\{ -\sqrt{1/3} , \sqrt{1/3} \}$ |
| 3 | $\{ 5/9 , 8/9 , 5/9 \}$ | $\{ -\sqrt{3/5} , 0, \sqrt{3/5} \}$ |

Table: Gauss' weights $w_s$ and nodes $\xi_s$ of quadrature for the numerical integration of integral $\int_{-1}^1 f(x) \, \mathrm{d}x = \sum_{s=1}^S w_s f(\xi_s) + \epsilon$. 

It is noteworthy to point out that Gauss' method of integration does not incorporate an evaluation of its kernel $K$ at the upper limit of integration, viz., at $t_N = N \, \mathrm{d}t$, where some kernels of interest become singular. Consequently, the error of approximation for this integrator will converge to a finite value with decreasing step size. [1]

## Solving the Resulting Differential Equations

Here a Volterra integral equation of the second kind, whose method of solution has just been laid out, describes a differential equation whose solution is sought. Specifically, the control $\boldsymbol{g}^{\prime}$ and response $\boldsymbol{f}^{\prime}$ functions are differential equations that are integrated here using a Backward Difference Formula (BDF), e.g., given an initial condition of $\boldsymbol{f}_0 = \boldsymbol{f}(t_0)$, one starts off with
$$
\boldsymbol{f}_1 = \boldsymbol{f}_0 + \tfrac{1}{2} \boldsymbol{f}^{\prime}_1 \, \mathrm{d} t \quad \text{for} \; n = 1
$$
and then, for the second step, one advances the solution with
$$
\boldsymbol{f}_2 = \tfrac{4}{3} \boldsymbol{f}_1 - \tfrac{1}{3} \boldsymbol{f}_0 + \tfrac{2}{3} \boldsymbol{f}^{\prime}_2 \, \mathrm{d}t \quad \text{for} \; n = 2
$$
while continuing thereafter with
$$
\boldsymbol{f}_n = \tfrac{18}{11} \boldsymbol{f}_{n-1} - \tfrac{9}{11} \boldsymbol{f}_{n-2} + \tfrac{2}{11} \boldsymbol{f}_{n-3} + \tfrac{6}{11} \boldsymbol{f}^{\prime}_n \, \mathrm{d}t \quad \text{for all} \; n > 2
$$
whose rates $\boldsymbol{f}^{\prime}_n$ are solutions to a Volterra integral equation of the second kind. A like strategy is used to integrate the control function $\boldsymbol{g}^{\prime}$.

## Software

This software uses the `PhysicalFields.jl` package.

### Memory Functions

Normalized memory functions $K(t)$, like those introduced above, are callable in software through the following functions. For those kernels that are not weakly singular
```
function BOX(systemOfUnits::String, time::PhysicalScalar, parameters::Tuple)::Tuple
function MCM(systemOfUnits::String, time::PhysicalScalar, parameters::Tuple)::Tuple
function MPL(systemOfUnits::String, time::PhysicalScalar, parameters::Tuple)::Tuple
function RFS(systemOfUnits::String, time::PhysicalScalar, parameters::Tuple)::Tuple
function SLS(systemOfUnits::String, time::PhysicalScalar, parameters::Tuple)::Tuple
```
and for those kernels that are weakly singular
```
function CCM(systemOfUnits::String, time::PhysicalScalar, parameters::Tuple)::Tuple
function FLS(systemOfUnits::String, time::PhysicalScalar, parameters::Tuple)::Tuple
function KWW(systemOfUnits::String, time::PhysicalScalar, parameters::Tuple)::Tuple
```
while other memory functions of one's own design can be introduced, too, provided they have an interface of:
```
function <myMemoryFunction>(systemOfUnits::String, time::PhysicalScalar, parameters::Tuple)::Tuple
```
Herein string `systemOfUnits` is (at present) either "SI" or "CGS", with `time` specifying the argument of `K(t)` sent to the generalized memory function for evaluation, and whose `parameters,` or material constants, are supplied via a Tuple.

The tuples that are to be supplied for the incorporated memory functions for creep listed above, i.e. the model's material constants, are:

| $K$ | parameters |
| --- | --- |
| BOX | $(τ_1, τ_2)$ |
| CCM | $(α, τ)$ |
| FLS | $(α, τ)$ |
| KWW | $(α, τ)$ |
| MCM | $(c_1, c_2, …, c_L, τ_1, τ_2, …, τ_L)$ |
| MPL | $(α, τ)$ |
| RFS | $(α, δ, τ)$ |
| SLS | $(τ,)$ |

All memory functions return a tuple. Specifically, they return tuple $(name, k, \tau)$ wherein $name$ is string specifying the name of the kernel, e.g. "BOX", while $k$ is the value of $K(t)$ and $\tau$ is its rate controlling characteristic time, both of which are instances of type `PhysicalFields.PhysicalScalar.`

### Constructing Weights of Quadrature for Volterra Integral Equations

The weights of quadrature for this, a solver for Volterra integral equationss of the second kind, can be created with a call to the following function
```
function normalizedQuadratureWeights(systemOfUnits::String,
                                     N::Integer,
                                     dTime::PhysicalScalar,
                                     kernel::Function,
                                     parameters::Tuple)::ArrayOfPhysicalScalars
```
where the `kernel` is any of the eight memory functions addressed above, or one of your own design. The supplied `kernel` function is called internally, whose arguments include `systemOfUnits` and `parameters,` with its argument for `time` being an integer multiple of step size `dTime,` which is the time-step size used by the solver, up to a time of `N dTime,` where `N` is the number of solution steps to be taken. The returned quadrature weights are contained within an instance of type `PhysicalFields.ArrayOfPhysicalScalars` whose length is `N.`

The returned quadrature weights are of type `PhysicalFields.ArrayOfPhysicalScalars,` and as such, they can be stored to a JSON file for a later retrieval. In fact, a call to `normalizedQuadratureWeights` will first check to see if they exist in a file, and if so they will be read in and returned; otherwise, they are computed and written to file before this function returns them.

# Solver for Volterra Integral Equations of the Second Kind

The exported solvers for Volterra integral equations are implementations of the abstract class
```
abstract type VolterraIntegralEquation end
```
and they come in three flavors: for scalar equations, for vector equations, and for tensor equations.

## For Scalar Equations:

```
struct VolterraIntegralScalarEquation <: VolterraIntegralEquation
    # Dimensioning fields
    dt::PhysicalScalar          # distance separating global integration nodes
    N::Int64                    # number of integration nodes in a solution path
    n::MInteger                 # current node along a solution path
    # Arrays of length N+1 holding control and response fields, and their rates
    f::ArrayOfPhysicalScalars   # array of integrated response function values
    f′::ArrayOfPhysicalScalars  # array of response function rates
    g::ArrayOfPhysicalScalars   # array of integrated control function values
    g′::ArrayOfPhysicalScalars  # array of control function rates
    t::ArrayOfPhysicalScalars   # array of times, the independent variable
    # Array of N normalized weights of quadrature for a product integral
    W::ArrayOfPhysicalScalars   # array holding quadrature weights
end
```

### Internal constructors:

There are two such constructors, i.e., the first constructor is
```
function VolterraIntegralScalarEquation(systemOfUnits::String,
                                        N::Int64,
                                        dt::PhysicalScalar,
                                        f₀::PhysicalScalar,
                                        g₀::PhysicalScalar,
                                        W::ArrayOfPhysicalScalars)
```
wherein `f₀` and `g₀` are initial values for the response variable and the control variable, respectively, while the second constructor is
```
function VolterraIntegralScalarEquation(dt::PhysicalScalar,
                                        N::Int64,
                                        n::MInteger, 
                                        f::ArrayOfPhysicalScalars,
                                        f′::ArrayOfPhysicalScalars,
                                        g::ArrayOfPhysicalScalars,
                                        g′::ArrayOfPhysicalScalars,
                                        t::ArrayOfPhysicalScalars,
                                        W::ArrayOfPhysicalTensors)
```
The former is usually used by a user to construct such a data structure, while the latter is used by `JSON3.jl` whenever a data structure of this type is to be recreated from a JSON file.

### Methods

To create a copy of an object of this type, call
```
function Base.:(copy)(vie::VolterraIntegralScalarEquation)::VolterraIntegralScalarEquation
```
or to create a deep copy of the object, call
```
function Base.:(deepcopy)(vie::VolterraIntegralScalarEquation)::VolterraIntegralScalarEquation
```

### Working with JSON files.

These objects are persistent in the sense that they can be stored-to and retrieved-from a file; specifically, one can call the following method to write an object to a JSON file.
```
function toFile(y::VolterraIntegralScalarEquation, json_stream::IOStream)
```
Or one can call the following method to read an object from a JSON file.
```
function fromFile(::Type{VolterraIntegralScalarEquation}, json_stream::IOStream)::VolterraIntegralScalarEquation
```
where an appropriate `json_stream` comes from functions
```
function PhysicalFields.openJSONReader(my_dir_path::String, my_file_name::String)::IOStream
function PhysicalFields.openJSONWriter(my_dir_path::String, my_file_name::String)::IOStream
function PhysicalFields.closeJSONStream(json_stream::IOStream)
```

### Calling the solver to advance a solution.

After a Volterra integral's data structure has been created, thereby initializing a problem, one can advance a solution along its path from node *n-1* to node *n* with *n* advancing to *N*. This is done by sequentially calling the following function:
```
function advance!(vie::VolterraIntegralScalarEquation, g′ₙ::PhysicalScalar, cₙ::PhysicalScalar)
```
wherein `vie` is an instance of type `VolterraIntegralScalarEquation.` Argument `g′ₙ` is a scalar that contains the known rate for the control function at time $t_n = t_{n-1} + \mathrm{d}t$, while argument `cₙ` is a known scalar field also evaluated at time $t_n$.

There are applications where a solution is to be secured iteratively, e.g., during an optimization problem. In such cases one can call the following method as many times as needed before advancing a solution to its next node along its path via a call to `advance!.` Such iterative refinements result from a call to
```
function update!(vie::VolterraIntegralScalarEquation, g′ₙ::PhysicalScalar, cₙ::PhysicalScalar)
```
where `g′ₙ` and `cₙ` are scalar fields that contain updated estimates for the control function and coefficient, e.g., as supplied by a finite element solver. One need not call method `update!` unless either the coefficient `cₙ`  or the control rate `g′ₙ` is not known explicitly, and therefore must be obtained iteratively.

## For Vector Equations:

```
struct VolterraIntegralVectorEquation <: VolterraIntegralEquation
    # Dimensioning fields
    dt::PhysicalScalar          # distance separating global integration nodes
    N::Int64                    # number of integration nodes in a solution path
    n::MInteger                 # current node along a solution path
    # Arrays of length N+1 holding control and response fields, and their rates
    f::ArrayOfPhysicalVectors   # array of integrated response function values
    f′::ArrayOfPhysicalVectors  # array of response function rates
    g::ArrayOfPhysicalVectors   # array of integrated control function values
    g′::ArrayOfPhysicalVectors  # array of control function rates
    t::ArrayOfPhysicalScalars   # array of times, the independent variable
    # Array of N normalized weights of quadrature for a product integral
    W::ArrayOfPhysicalScalars   # array holding quadrature weights
end
```

### Internal constructors:

There are two such constructors, i.e., the first constructor is
```

    function VolterraIntegralVectorEquation(systemOfUnits::String,
                                            N::Int64,
                                            dt::PhysicalScalar,
                                            f₀::PhysicalVector,
                                            g₀::PhysicalVector,
                                            W::ArrayOfPhysicalScalars)
```
wherein `f₀` and `g₀` are initial values for the response variable and the control variable, respectively, while the second constructor is
```
function VolterraIntegralVectorEquation(dt::PhysicalScalar,
                                        N::Int64,
                                        n::MInteger,
                                        f::ArrayOfPhysicalVectors,
                                        f′::ArrayOfPhysicalVectors,
                                        g::ArrayOfPhysicalVectors,
                                        g′::ArrayOfPhysicalVectors, 
                                        t::ArrayOfPhysicalScalars,
                                        W::ArrayOfPhysicalScalars)
```
The former is used to initially construct such a data structure, while the latter is used by `JSON3.jl` whenever a data structure of this type is to be recreated from a JSON file.

### Methods

To create a copy of an object of this type, call
```
function Base.:(copy)(vie::VolterraIntegralVectorEquation)::VolterraIntegralVectorEquation
```
or to create a deep copy of the object, call
```
function Base.:(deepcopy)(vie::VolterraIntegralVectorEquation)::VolterraIntegralVectorEquation
```

### Working with JSON files.

These objects are persistent in the sense that they can be stored-to and retrieved-from a file; specifically, one can call the following method to write an object to a JSON file.
```
function toFile(y::VolterraIntegralVectorEquation, json_stream::IOStream)
```
Or one can call the following method to read an object from a JSON file.
```
function fromFile(::Type{VolterraIntegralVectorEquation}, json_stream::IOStream)::VolterraIntegralVectorEquation
```
where an appropriate `json_stream` comes from functions
```
function PhysicalFields.openJSONReader(my_dir_path::String, my_file_name::String)::IOStream
function PhysicalFields.openJSONWriter(my_dir_path::String, my_file_name::String)::IOStream
function PhysicalFields.closeJSONStream(json_stream::IOStream)
```

### Calling the solver to advance a solution.

After a Volterra integral's data structure has been created, thereby initializing a problem, one can advance a solution along its path from node *n-1* to node *n* by sequentially calling the following function:
```
function advance!(vie::VolterraIntegralVectorEquation, g′ₙ::PhysicalVector, cₙ::PhysicalScalar)
```
wherein `vie` is an instance of type `VolterraIntegralVectorEquation.` Argument `g′ₙ` is a vector that contains the known rate for the control function at time $t_n = t_{n-1} + \mathrm{d}t$, while argument `cₙ` is a known scalar field also evaluated at time $t_n$.

There are applications where a solution is to be secured iteratively, e.g., during an optimization problem. In such cases one can call the following method as many times as needed before advancing a solution to its next node along its path via a call to `advance!.` Such iterative refinements result from a call to
```
function update!(vie::VolterraIntegralVectorEquation, g′ₙ::PhysicalVector, cₙ::PhysicalScalar)
```
where `g′ₙ` is a vector and `cₙ` is a scalar. These fields contain updated estimates for the control function and/or the coefficient, e.g., as supplied by a finite element solver. One need not call method `update!` unless either the coefficient `cₙ`  or the control rate `g′ₙ` is not known explicitly.

## For Tensor Equations:

```
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
    W::ArrayOfPhysicalScalars   # array holding the quadrature weights weights
end
```

### Internal constructors:

There are two such constructors, i.e., the first constructor is
```
function VolterraIntegralTensorEquation(systemOfUnits::String, 
                                        N::Int64,
                                        dt::PhysicalScalar,
                                        f₀::PhysicalTensor,
                                        g₀::PhysicalTensor,
                                        W::ArrayOfPhysicalScalars)
```
wherein `f₀` and `g₀` are initial values for the response variable and the control variable, respectively, while the second constructor is
```
function VolterraIntegralTensorEquation(dt::PhysicalScalar,
                                        N::Int64,
                                        n::MInteger,
                                        f::ArrayOfPhysicalTensors,
                                        f′::ArrayOfPhysicalTensors,
                                        g::ArrayOfPhysicalTensors,
                                        g′::ArrayOfPhysicalTensors,
                                        t::ArrayOfPhysicalScalars,
                                        W::ArrayOfPhysicalScalars)
```
The former is used to initially construct such a data structure, while the latter is used by `JSON3.jl` whenever a data structure of this type is to be recreated from a JSON file.

### Methods

To create a copy of an object of this type, call
```
function Base.:(copy)(vie::VolterraIntegralTensorEquation)::VolterraIntegralTensorEquation
```
or to create a deep copy of the object, call
```
function Base.:(deepcopy)(vie::VolterraIntegralTensorEquation)::VolterraIntegralTensorEquation
```

### Working with JSON files.

These objects are persistent in the sense that they can be stored-to and retrieved-from a file; specifically, one can call the following method to write an object to a JSON file.
```
function toFile(y::VolterraIntegralTensorEquation, json_stream::IOStream)
```
Or one can call the following method to read an object from a JSON file.
```
function fromFile(::Type{VolterraIntegralTensorEquation}, json_stream::IOStream)::VolterraIntegralTensorEquation
```
where an appropriate `json_stream` comes from functions
```
function PhysicalFields.openJSONReader(my_dir_path::String, my_file_name::String)::IOStream
function PhysicalFields.openJSONWriter(my_dir_path::String, my_file_name::String)::IOStream
function PhysicalFields.closeJSONStream(json_stream::IOStream)
```

### Calling the solver to advance a solution.

After a Volterra integral's data structure has been created, thereby initializing a problem, one can advance a solution along its path from node *n-1* to node *n* by sequentially calling the following function:
```
function advance!(vie::VolterraIntegralTensorEquation, g′ₙ::PhysicalTensor, cₙ::PhysicalScalar)
```
wherein `vie` is an instance of type `VolterraIntegralTensorEquation.` Argument `g′ₙ` is a tensor that contains the known rate for the control function at time $t_n = t_{n-1} + \mathrm{d}t$, while argument `cₙ` is a known scalar field also evaluated at time $t_n$.

There are applications where a solution is to be secured iteratively, e.g., during an optimization problem. In such cases one can call the following method as many times as needed before advancing a solution to its next node along its path via a call to `advance!.` Such iterative refinements result from a call to
```
function update!(vie::VolterraIntegralTensorEquation, g′ₙ::PhysicalTensor, cₙ::PhysicalScalar)
```
where `g′ₙ` is a tensor field and `cₙ` is a scalar field that contain updated estimates for the control function and a coefficient, e.g., as supplied by a finite element solver. One need not call method `update!` unless either the coefficient `cₙ`  or the control rate `g′ₙ` is not known explicitly.

## References

1) Braß, H., "On the Principle of Avoiding the Singularity in Quadrature," Zeitschrift für angewandte Mathematik und Mechanik, 75 (1995), S617-S618.

2) Caputo, M. and Mainardi, F., "Linear models of dissipation in anelastic solids," *Rivista del Nuoro Cimento*, **1** (1971), 161-198.

3) Caputo, M. and Mainardi, F., "A new dissipation model based on memory mechanism," *Pure and Applied Geophysics*, **91** (1971), 134-147.

4) Cole, K.S. and Cole, R.H., "Dispersion and absorption in dielectrics I. Alternating current characteristics," *Journal of Chemical Physics*, **9** (1941), 342-351.

5) Cole, K.S. and Cole, R.H., "Dispersion and absorption in dielectrics II. Direct current characteristics," *Journal of Chemical Physics*, **10** (1942), 98-105.

6) Freed, A.D., *Soft Solids: A primer to the theoretical mechanics of materials*, Modeling and Simulation in Science, Engineering and Technology. Basel: Birkhäuser, 2014.

7) Freed, A.D. and Clayton, J.D., *Application of Laplace Stretch in Alveolar Mechanics: A Case Study of Blast and Blunt Trauma*, in preparation.

8) Freed, A.D. and Rajagopal, K.R. "A viscoelastic model for describing the response of biological fibers," *ACTA Mechanica*, **227** (2016), 3367-3380.

9) Fung, Y.-C., "Biorheology of Soft Tissues," Biorheology, 10 (1973), 139-155.

10) Kohlrausch, R., "Ueber das Dellmann'sche Elektrometer," *Annalen der Physik und Chemie*, **72** (1847), 353-405.

11) Maxwell, J.C., "On the dynamical theory of gases," *Philosophical Transactions of the Royal Society, London*, **157** (1867), 49-88.

12) Neubert, H.K., "A simple model representing internal damping in solid materials," *The Aeronautical Quarterly*, **14** (1963), 187-210.

13) Volterra, V. *Theory of functionals and of integral and integro-differential equations*. Glasgow: Blackie and Son, 1930.

14) Williams, G. and Watts, D.C., "Non-symmetrical dielectric relaxation behaviour arising from a simple empirical decay function," *Transactions of the Faraday Society*, **66** (1970), 80-85.

15) Williams, M.L., "Structural analysis of viscoelastic materials," *AIAA Journal*, **2** (1964), 785-808.

16) Young, A., "Approximate product integration," *Proceedings of the Royal Society, London*, **A-224** (1954), 552-561.

17) Zener, C., *Elasticity and Anelasticity of Metals*. Chicago: University of Chicago Press, 1948.

## Version History

### Version 0.2.0

A different and simpler solver is employed here. The previous algorithm expanded the kernel in a Taylor expansion to construct moments that were used as quadrature weights [16]. This resulted in an artifact--a numerical ringing--that would occur at the beginning of some solutions. The new solver does not include such a sophistication, and seems to be better behaved for it.

### Version 0.1.3

Bug fixes.

### Version 0.1.2

The coefficient c within a Volterra integral equation is now taken to be a function, assigned during construction. It is introduced via the solver, viz., as an argument in methods advance! and update!, instead of as a field in the data structure.

### Version 0.1.1

Rearranged order of fields in the three viscoelastic data structures so that they conform with the order of other data structures used for solvers in this family of software.

### Version 0.1.0

Initial release.