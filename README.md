# VolterraIntegralEquations

This module provides a solver for Volterra integral equations of the second kind [13]. Specifically, this solver solves an integro-differential equation of the form
$$
\boldsymbol{f}^{\prime}(t) = \boldsymbol{g}^{\prime}(t) - c(t) \int_0^t K(t-τ) \, \boldsymbol{f}^{\prime}(τ) \, \mathrm{d}τ
$$
where the lower limit of integration is an initial time, taken to be 0, while the upper limit of integration is some current time $t$, with $\tau$ denoting a dummy variable of integration. For illustrative purposes, the independent variable in this document is taken to be time, but this is not necessary. 

This equation becomes a Volterra integral equation of the first kind whenever $\boldsymbol{g}^{\prime}$ is $\boldsymbol{0}$.

This Volterra integral equation describes an ordinary differential equation that, in turn, must be solved, i.e., it is an integro-differential equation. Here derivative $\boldsymbol{g}^{\prime} = \mathrm{d}\boldsymbol{g}/\mathrm{d}t$ is a known rate of some control function $\boldsymbol{g}$, while derivative $\boldsymbol{f}^{\prime} = \mathrm{d}\boldsymbol{f}/\mathrm{d}t$ is an unknown rate of a response function $\boldsymbol{f}$ to be solved for, whose evolution is characterized by a Volterra integral equation. Coefficient $c(t)$ and kernel $K(t-τ)$ are also taken to be known quantities.

Functions $\boldsymbol{f}$ and $\boldsymbol{g}$ have the same physical units. They may be scalars, vectors or tensors in stature. Coefficient $c$ is a dimensionless, scalar-valued function. Kernel $K$ is a memory function, which is the derivative of a generalized creep function. It is a scalar-valued function with physical units of $t^{-1}$, e.g., reciprocal time whenever $t$ has units of time.

To use this module, you will need to add the following repositories to your project:

```
using Pkg
Pkg.add(url = "https://github.com/AlanFreed/PhysicalFields.jl")
Pkg.add(url = "https://github.com/AlanFreed/VolterraIntegralEquations.jl")
```

This software was written to accompany a book the authors are writing [7].

## Example: A Viscoelastic Fiber

A linear viscoelastic fiber can be modeled as a Volterra integral equation of the second kind, viz., [6]
$$
\frac{\mathrm{d} \sigma}{\mathrm{d} t} = \mathcal{E}(t) \left(
	\frac{\mathrm{d} \epsilon}{\mathrm{d} t} -
	\frac{E_0 - E_{\infty}}{E_0 \, E_{\infty}}
	\int_0^t K (t - \tau) \, \frac{\mathrm{d} \sigma}{\mathrm{d}\tau} \,
	\mathrm{d}\tau \right)
$$
where $\sigma$ is stress and $\epsilon$ is strain, and where $\mathcal{E}$,  $E_{\infty}$, and $E_0$ are its local tangent, rubbery, and glassy elastic moduli, respectively, the former being a possible function of time, while the latter two are material constants. Kernel $K$ is a positive, monotonic-decreasing function known as a memory function whose units are reciprocal time.  The first term on the right-hand side of the above equation provides for a glassy elastic change in stress that is attenuated by the second term, which introduces a viscous loss to this change in stress.

The elastic tangent modulus $\mathcal{E}$ for a linear Hookean fiber is
$$
\mathcal{E} = E_0
$$
wherein $E_0$ denotes its glassy modulus.

In contrast, for a nonlinear biologic fiber, its tangent modulus associates with a nonlinear compliance of
$$
\frac{1}{\mathcal{E}} = \frac{1}{E_0} + \frac{\beta}{E_r \beta + 2(\sigma - \sigma_r)} \sqrt{\frac{\beta E_r}{\beta E_r + 2(\sigma - \sigma_r)}}
$$
wherein $\sigma_r$ denotes a residual stress, and $\beta$ designates a limiting state for an internal strain, a strain caused by molecular reconfiguration.  The elastic tangent modulus associated with a fiber's strain-free reference configuration is designated as $E_r$ ($> 0$). Within a biologic fiber's linear region of response, wherein strain is caused by molecular stretching, moduli $E_\infty$ ($> E_r$) and $E_0$ ($> E_{\infty}$) denote its rubbery and glassy moduli, respectively. The resulting elastic tangent modulus $\mathcal{E}$ is thereby bounded by $E_r$ from below and $E_0$ from above.

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

A selection of positive, monotonic-decreasing, viscoelastic, memory functions $K(t)$, whose units are reciprocal time, have been preprogrammed into this software. These are representative of the many kernel functions that have been proposed in the literature. These memory functions are the derivatives of creep functions, the latter of which are more commonly found in the literature. Consequently, these characteristic times associate with creep (not stress relaxation). See Freed [6,8] for a discussion of these functions.

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
wherein $\tau$ is a characteristic time, and $\alpha \in (0, 1]$ is a fractional order of evolution, with $E_{\alpha,\beta}(t)$ being the two-parameter Mittag-Leffler function. The FLS model contains the SLS model as a special case; specifically, they are equivalent whenever $\alpha = 1$. Mainardi's memory kernel is weakly singular.

4) **KWW**: **K**ohlrausch's [10] and **W**illiams & **W**atts' [14] stretched exponential has a memory function of
$$
K(t) = \alpha \, \left( \frac{t}{\tau} \right)^{\alpha} \;
\frac{\exp \bigl( -(t/\tau)^{\alpha} \bigr)}{t}
\quad \text{with} \quad
K(0) = \infty
$$
wherein $\tau$ is a characteristic time and $\alpha \in (0,1]$ is an exponent for the power of the argument in the exponential. The KWW model contains the SLS model as a special case; specifically, they are equivalent whenever $\alpha = 1$. The KWW memory kernel is weakly singular whenever $\alpha \in (0,1)$.

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
which effectively shifts the singularity of an FLS kenel a short distance into negative time so the RFS creep compliance equates with the SLS creep compliance at time 0. As a consequence, the RFS memory kernel is finite valued throughout its domain of application.

8) **SLS**: Zener's [17] **S**tandard **L**inear **S**olid has a Maxwell-Debye kernel resulting in a memory function of
$$
K(t) = \frac{\mathrm{exp}(-t / \tau)}{\tau}
\quad \text{with} \quad
K(0) = \frac{1}{\tau}
$$
wherein $\tau$ is a characteristic time.

A capability is provided in this software for the user to define and use their own memory function, too, if it is other than one of those listed above.

## Numerical Method

This solver implements a numerical method developed in an appendix of a book that is currently being written by the authors of this software [6]. It builds upon another numerical method developed by Young [16] which comes from the actuarial sciences literature.

Solutions to Volterra integrals of the second kind can be advanced through a block-by-block algorithmic strategy of
$$
\boldsymbol{f}^{\prime}(t_1) = \boldsymbol{g}^{\prime}(t_1) - c(t_1) \int_{0}^{t_1} K(t_1 - \tau) \, \boldsymbol{f}^{\prime}(\tau) \, \mathrm{d}\tau
$$
$$
\boldsymbol{f}^{\prime}(t_2) = \boldsymbol{g}^{\prime}(t_2) - c(t_2) \int_{0}^{t_2} K(t_2 - \tau) \, \boldsymbol{f}^{\prime}(\tau) \, \mathrm{d}\tau 
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
wherein, e.g., the product integral $\int_{0}^{t_1} K(t_2 - \tau) \, \boldsymbol{f}^{\prime}(\tau) \, \mathrm{d}\tau$ has a kernel $K$ whose fixed time $t_2$ lies outside the interval of integration, in this case $[0, t_1]$, over which time $\tau$ spans, while the product integral $\int_{t_1}^{t_2} K(t_2 - \tau) \, \boldsymbol{f}^{\prime}(\tau) \, \mathrm{d}\tau$ has a kernel $K$ whose fixed time, in this case $t_2$, is now the upper limit of integration. These two integrals will have different rules of quadrature representing them.

### Quadrature Rule

Many of the memory kernels used in practice are weakly singular in that they become unbounded at the upper limit of integration. Braß [1] has shown that open quadrature methods (those that do not contain the interval's endpoints as nodes of quadrature, and therefore avoid its singularity) have error functions that converge; whereas, every closed quadrature method (those that contain the singular point) has an error function that diverges. Consequently, open quadrature methods must be used to approximate such integrals, like the method being considered here.

There are product integrals of the Fredholm type, viz.,
$$
\int_{t_{l-1}}^{t_l} K(t_n - \tau) \, \boldsymbol{f}^{\prime}(\tau) \, \mathrm{d}\tau =
	\sum_{j=1}^J W_j (K, \mathrm{d}t) \, \boldsymbol{f}^{\prime}(t_{l-1} + t_j)  + \boldsymbol{\varepsilon}
$$
where $l = 0,1,2, \ldots, n \! - \! 1,$ with $t_j \in (0, \mathrm{d}t)$ representing the $J$ local nodes of quadrature for a selected integration rule, which are taken to be an open set over the interval, while recalling that $\mathrm{d}t = t_l - t_{l-1}$. Here $W_j$ denotes a weight of quadrature for which $\boldsymbol{\varepsilon}$ is its local truncation error.

There are also product integrals of the Volterra type, viz.,
$$
\int_{t_{n-1}}^{t_n} K(t_n - \tau) \, \boldsymbol{f}^{\prime}(\tau) \, \mathrm{d} \tau = \sum_{j=1}^J W_j (K, \mathrm{d}t) \, \boldsymbol{f}^{\prime}(t_{n-1} + t_j) + \boldsymbol{\varepsilon}
$$
where the upper limit of integration $t_n$ now appears as an argument in the kernel function $K$.

These two quadrature rules result in similar, but different, formula for their weights of quadrature.

In practice, it is considered that the global interval of integration $\mathrm{d}t$ is sufficiently small so that its local truncation error $\boldsymbol{\varepsilon}$ can be safely neglected. Even so, the user ought to be cognizant of this consideration whenever one is assessing the validity of one's solution. A theorem by Braß [1] assures that the error converges with a refinement of nodal mesh density, and therefore, the solution converges, too.

### Block-by-Block Solution Strategy

Solutions to Volterra integral equations of the second kind take on the form of a sequentially solved linear equation
$$
\boldsymbol{f}^{\prime}_{n} = \Bigl( \boldsymbol{I} + c_n \boldsymbol{W}^{\mathsf{T}}_{\!1} \Bigr)^{-1} \left( \boldsymbol{g}^{\prime}_n - c_n \sum_{m=1}^{n-1} \boldsymbol{W}^{\mathsf{T}}_{\!n-m+1} \boldsymbol{f}^{\prime}_{m} \right) ,
\qquad n = 1, 2, \ldots, N
$$
where $N$ specifies the number of global integration steps, called nodes, that are required to traverse a solution's path, whose weights of quadrature are $\boldsymbol{W}_n$, with weight $\boldsymbol{W}_1$ being distinct in form from all the others. In order for this algorithm to work, it is necessary that the $J \times J$ matrix $\boldsymbol{I} + c_n \boldsymbol{W}_1^{\mathsf{T}}$ not be singular so that its inverse exists, with $c_n = c(t_n)$. 

For the method implemented here, there are three local nodes of integration per single global node, i.e., $J = 3$, and as such, the control function $\boldsymbol{g}^{\prime}_n$ and the forcing function $\boldsymbol{f}^{\prime}_n$ are each arrays of length 3, while the weights of quadrature $\boldsymbol{W}_n$ are $3 \times 3$ matrices. Specifically, an open, mid-point, quadrature rule is selected wherein
$$
t_j = \{ \tfrac{1}{6} \mathrm{d}t , \tfrac{1}{2} \mathrm{d}t , \tfrac{5}{6} \mathrm{d}t \}
$$
so that
$$
\boldsymbol{f}^{\prime}_n = \left\{ \begin{matrix}
\boldsymbol{f}^{\prime}(t_{n,1}) \\
\boldsymbol{f}^{\prime}(t_{n,2}) \\
\boldsymbol{f}^{\prime}(t_{n,3}) \\
\end{matrix} \right\}
= \left\{ \begin{matrix}
\boldsymbol{f}^{\prime}(t_{n-1} \! + \! \tfrac{1}{6} \mathrm{d}t) \\
\boldsymbol{f}^{\prime}(t_{n-1} \! + \! \tfrac{1}{2} \mathrm{d}t) \\
\boldsymbol{f}^{\prime}(t_{n-1} \! + \! \tfrac{5}{6} \mathrm{d}t) \\
\end{matrix} \right\}
$$
and
$$
\boldsymbol{g}^{\prime}_n = \left\{ \begin{matrix}
\boldsymbol{g}^{\prime}(t_{n,1}) \\
\boldsymbol{g}^{\prime}(t_{n,2}) \\
\boldsymbol{g}^{\prime}(t_{n,3}) \\
\end{matrix} \right\}
= \left\{ \begin{matrix}
\boldsymbol{g}^{\prime}(t_{n-1} \! + \! \tfrac{1}{6} \mathrm{d}t) \\
\boldsymbol{g}^{\prime}(t_{n-1} \! + \! \tfrac{1}{2} \mathrm{d}t) \\
\boldsymbol{g}^{\prime}(t_{n-1} \! + \! \tfrac{5}{6} \mathrm{d}t) \\
\end{matrix} \right\}
$$
with weights of quadrature being described by
$$
\boldsymbol{W}_n = \boldsymbol{X}^{-1} \boldsymbol{\mu}_n
\quad \text{wherein} \quad
\boldsymbol{X} = 
\begin{bmatrix}
1 & 1 & 1 \\
-1 & 0 & 1 \\
1 & 0 & 1
\end{bmatrix}
\quad \text{so} \quad
\boldsymbol{X}^{-1} = \frac{1}{2} \begin{bmatrix}
0 & -1 & 1 \\
2 & 0 & -2 \\
0 & 1 & 1
\end{bmatrix}
$$
where $\boldsymbol{X}$ is a Vandermonde matrix (Young [16] calls it the alternant matrix) and $\boldsymbol{\mu}_n$ is the $n^{\text{th}}$ moment matrix, which are consequences of expressing $\boldsymbol{f}^{\prime}$ in a Taylor series expanded about the midpoint to its local span of integration.

The global nodes of integration associate with times $t_n$ where $t_n = t_{n-1} + \mathrm{d}t \; \forall \; n$ given $t_0 = 0$, with $\mathrm{d}t$ being a distance separating neighboring global nodes. The local nodes of integration associate with times $t_{n,j}$ where $j = 1,2,3$ such that $t_{n,1} = t_{n-1} + \tfrac{1}{6} \mathrm{d}t$, $t_{n,2} = t_{n-1} + \tfrac{1}{2} \mathrm{d}t$ and $t_{n,3} = t_{n-1} + \tfrac{5}{6} \mathrm{d}t$, with a distance of $\tfrac{1}{3} \mathrm{d}t$ separating neighboring local nodes. These local nodes do not contain any global nodes, i.e., it uses an open quadrature method. Consequently, this solution strategy can, in principle, be applied to kernels that are singular at the upper limit of integration (like memory kernels: CCM, FLS and KWW).

### Moment Matrices

The greatest expense in implementing this numerical method is often in the computing its moment matrices. Fortunately, once gotten they can be reused in future solutions. There is a more efficient algorithm for implementing an exponential kernel that is based upon its recursive property [6], but that algorithm is restricted to that kernel alone. The algorithm presented below is applicable to all types of kernel functions, and is therefore versatile.

The moment matrices are solutions to integral equations. For $n=1$, moment $\boldsymbol{\mu}_1$ describes a $3 \times 3$ matrix whose elements are solutions to the Volterra integral
$$
\mu_{ij,1} = \frac{1}{h^{i-1}} \int_0^{(j - 1/2)h} \bigl( \tau - \tfrac{1}{2} (j - 1/2) h \bigr)^{i-1} \, K \bigl( (j - 1/2) h - \tau \bigr) \, \mathrm{d} \tau
$$
where $h = \tfrac{1}{3} \mathrm{d}t$ and $i,j=1,2,3$. The remaining moment matrices $\boldsymbol{\mu}_n$, where $n=2, 3, \ldots, N$, describe $3 \times 3$ matrices whose elements are solutions to the Fredholm integral
$$
\mu_{ij,n} = \frac{1}{h^{i-1}} \int_0^{\mathrm{d}t} \bigl( \tau - \tfrac{1}{2} \mathrm{d}t \bigr)^{i-1} \, K \bigl( (n-1) \mathrm{d}t + (j - 1/2) h - \tau \bigr) \, \mathrm{d} \tau
$$
These are integrals of the kernel function $K$ scaled by moment arms measured from the midpoints of their spans of integration. These moment arms are raised to powers of $i \! - \! 1$. They arise because of the Taylor series expansion imposed.

#### Gauss Quadrature

Employing Gauss' quadrature rule
$$
\int_0^t \mathfrak{f}(\tau) \, \mathrm{d}\tau = \frac{t}{2} \int_{-1}^1 \mathfrak{f} \bigl( \tfrac{t}{2} ( 1 + x ) \bigr) \, \mathrm{d} x \approx \frac{t}{2} \sum_{s=1}^S w_s \, \mathfrak{f} \bigl( \tfrac{t}{2} ( 1 + x_s ) \bigr)
$$
to the above moment integrals, selecting $S=3$, produces weights $w_s$ and nodes $x_s$ of quadrature
$$
w_s = \left\{ 5/9, \, 8/9, \, 5/9 \right\}^{\mathsf{T}}
\quad \text{and} \quad
x_s = \left\{ -\sqrt{3/5} , \, 0 , \, \sqrt{3/5} \right\}^{\mathsf{T}}
$$
Consequently, the first moment $\boldsymbol{\mu}_1$ is approximated by a quadrature rule of
$$
\mu_{ij,1} = \frac{(j \! - \! 1/2) \, \mathrm{d}t}{6} \sum_{s=1}^3 w_s m_{ij,s} K \left( \tfrac{1}{6} (j \! - \! 1/2) (1 \! - \! x_s)^{\vphantom{|}} \mathrm{d}t \right) , \quad i,j = 1,2,3
$$
wherein, for $s=1$
$$
m_{ij,1} = \begin{bmatrix}
	1 & 1 & 1 \\
	-\sqrt{3/80} & -\sqrt{27/80} & -\sqrt{15/16} \\
	3/80 & 27/80 & 15/16
\end{bmatrix}
$$
while, for $s=2$
$$
m_{ij,2} = \begin{bmatrix}
	1 & 1 & 1 \\
	0 & 0 & 0 \\
	0 & 0 & 0
\end{bmatrix} 
$$
where the zeros follow because $x_2$ resides at the midpoint, i.e., it has no moment arm, and for $s=3$
$$
m_{ij,3} = \begin{bmatrix}
	1 & 1 & 1 \\
	\sqrt{3/80} & \sqrt{27/80} & \sqrt{15/16} \\
	3/80 & 27/80 & 15/16
\end{bmatrix}
$$
Similarly, the remaining moments $\boldsymbol{\mu}_n$, $n = 2,3,\ldots,N$, are approximated by a quadrature rule of
$$
\mu_{ij,n} = \frac{\mathrm{d}t}{2} \sum_{s=1}^3 w_s v_{i,s} K \left( \bigl( n - \tfrac{1}{3} (5 - j) - \tfrac{1}{2} x_s \bigr)^{\vphantom{|}} \mathrm{d}t \right) , \quad i,j = 1,2,3
$$
wherein, for $s=1$
$$
v_{i,1} = \left\{ 1 , \, -\sqrt{27/20} , \, 27/20 \right\}^{\mathsf{T}}
$$
while, for $s=2$
$$
v_{i,2} = \left\{ 1, \, 0, \, 0 \right\}^{\mathsf{T}}
$$
and, for $s=3$
$$
v_{i,3} = \left\{ 1 , \, \sqrt{27/20} , \, 27/20 \right\}^{\mathsf{T}} 
$$
that collectively weigh an effect caused by the moment arms within a Taylor expansion.

### History Truncation

For those kernels $K$ that are monotonic-decreasing functions, a number $N_{\max}$ exists beyond which point memory of the past effectively fades away. How one assigns $N_{\max}$ will depend upon the kernel $K$ (specifically, its characteristic time), the global step size $\mathrm{d}t$, and the accuracy sought in a solution. Considering that such an $N_{\max}$ exists, then, for integration steps where $n \le N_{\max}$, a solution $\boldsymbol{f}^{\prime}_n$ advances along its path according to the linear equation
$$
\boldsymbol{f}^{\prime}_n = \Bigl( \boldsymbol{I} + c_n \boldsymbol{W}^{\mathsf{T}}_1 \Bigr)^{-1} \left( \boldsymbol{g}^{\prime}_n - c_n \sum_{m=1}^{n-1} \boldsymbol{W}^{\mathsf{T}}_{n-m+1} \boldsymbol{f}^{\prime}_{m} \right)
$$
that whenever $n > N_{\max}$ will continue to advance along its path, but now according to the linear equation
$$
\boldsymbol{f}^{\prime}_n = \Bigl( \boldsymbol{I} + c_n \boldsymbol{W}^{\mathsf{T}}_1 \Bigr)^{-1} \left( \boldsymbol{g}^{\prime}_n - c_n \sum_{m=1}^{N_{\max}-1} \boldsymbol{W}^{\mathsf{T}}_{N_{\max}-m+1} \boldsymbol{f}^{\prime}_{m+n-N_{\max}} \right)
$$
which removes forcing functions from $\boldsymbol{f}^{\prime}_1$ through $\boldsymbol{f}^{\prime}_{n-N_{\max}}$ from its summation history. These memories are so distant that they have effectively been forgotten.

## Solving the Resulting Differential Equations

Here a Volterra integral equation of the second kind, whose method of solution has just been laid out, describes a differential equation whose solution is sought. Specifically, the control $\boldsymbol{g}^{\prime}$ and response $\boldsymbol{f}^{\prime}$ rate functions can both be integrated using a Newton-Cotes formula of the form
$$
\boldsymbol{g}_{n} = \boldsymbol{g}_{n-1} + \frac{\mathrm{d}t}{8} \left( 3 \boldsymbol{g}^{\prime} (t_{n,1}) + 2 \boldsymbol{g}^{\prime} (t_{n,2})+ 3 \boldsymbol{g}^{\prime} (t_{n,3}) \right)
$$
and
$$
\boldsymbol{f}_{n} = \boldsymbol{f}_{n-1} + \frac{\mathrm{d}t}{8} \left( 3 \boldsymbol{f}^{\prime} (t_{n,1}) + 2 \boldsymbol{f}^{\prime} (t_{n,2})+ 3 \boldsymbol{f}^{\prime} (t_{n,3}) \right)
$$
with the $\boldsymbol{g}^{\prime}$ being known and the $\boldsymbol{f}^{\prime}$ being solutions to a Volterra integral equation. Time $t_{n}$ associates with the next global node of integration, while times $t_{n,j} = t_{n-1} + (\tfrac{1}{6} + \tfrac{1}{3}(j-1)) \mathrm{d}t$, $j=1,2,3$, associate with the three local nodes of integration spanning $[t_{n-1} , t_{n}]$.

There are applications where the control function may depend upon either or both the control and/or response functions, e.g., the biologic fiber model in the viscoelastic example described above, in which case one will also need solutions to be gotten at the three local nodes of integration, too, viz.,
$$
\boldsymbol{f}_{n,1} = \boldsymbol{f}_{n-1} + \frac{\mathrm{d}t}{72} \left(
17 \boldsymbol{f}^{\prime} (t_{n,1}) -7 \boldsymbol{f}^{\prime} (t_{n,2})+ 2 \boldsymbol{f}^{\prime} (t_{n,3}) \right)
$$
$$
\boldsymbol{f}_{n,2} = \boldsymbol{f}_{n-1} + \frac{\mathrm{d}t}{8} \left(
3 \boldsymbol{f}^{\prime} (t_{n,1}) + \boldsymbol{f}^{\prime} (t_{n,2}) \right)
$$
$$
\boldsymbol{f}_{n,3} = \boldsymbol{f}_{n-1} + \frac{5 \, \mathrm{d}t}{72} \left( 5 \boldsymbol{f}^{\prime} (t_{n,1}) + 5 \boldsymbol{f}^{\prime} (t_{n,2})+ 2 \boldsymbol{f}^{\prime} (t_{n,3}) \right)
$$
with like expressions establishing $\boldsymbol{g}_{n,j}$. The method of Newton-Cotes provides these formulæ, too.

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

1) BOX: parameters = $(τ_1, τ_2)$
2) CCM: parameters = $(α, τ)$
3) FLS: parameters = $(α, τ)$
4) KWW: parameters = $(α, τ)$
5) MCM: parameters = $(c_1, c_2, …, c_L, τ_1, τ_2, …, τ_L)$
6) MPL: parameters = $(α, τ)$
7) RFS: parameters = $(α, δ, τ)$
8) SLS: parameters = $(τ,)$

All memory functions return a tuple. Specifically, they return tuple $(k, \tau)$ wherein $k$ is the value of $K(t)$ and $\tau$ is its rate controlling characteristic time, both of which are instances of type `PhysicalFields.PhysicalScalar.`

### Constructing Weights of Quadrature for Volterra Integral Equations

The weights of quadrature for this, a solver of Young's [16] method, can be created with a call to the following function
```
function normalizedQuadratureWeights(systemOfUnits::String,
                                     dTime::PhysicalScalar,
                                     parameters::Tuple,
                                     kernel::Function,
                                     Nₘₐₓ::Integer,
                                     significantFigures::Integer=5)::ArrayOfPhysicalTensors

```
where the `kernel` is any of the eight memory functions addressed above, or one of your own design. This function is called internally, whose arguments include `systemOfUnits` and `parameters,` with its argument for `time` being an integer multiple of step size `dTime,` which is the global time-step size of the solver. The returned quadrature weights are contained within an instance of type `PhysicalFields.ArrayOfPhysicalTensors` whose length is the lesser value of `Nₘₐₓ` and that length determined from the characteristic time $\tau$ supplied by the `kernel` at the accuracy sought. This is established via argument `significantFigues,` which is bound to the integer interval [3, 9] with 5 being its default value. For example, at this default, the array of quadrature weights will have a length of $\min (N_{\text{max}}, N_t)$ where $N_t$ is that value whereat $K(N_t \, \mathrm{d}t) < 10^{-5}$.

The returned quadrature weights are of type `PhysicalFields.ArrayOfPhysicalTensors,` and as such, they can be stored to a JSON file for a later retrieval. 

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
    N::Integer                  # number of integration nodes in a solution path
    Nₘₐₓ::Integer               # maximum number of nodes whose history is kept
    n::MInteger                 # current node along a solution path
    # Arrays of length N+1 holding integrated variable rates at the global nodes
    f::ArrayOfPhysicalScalars   # array of integrated response function values
    g::ArrayOfPhysicalScalars   # array of integrated control function values
    t::ArrayOfPhysicalScalars   # array of times, the independent variable
    # Array of length 3N holding response function rates at over local intervalS
    f′::ArrayOfPhysicalScalars  # history array of response function rates
    # Array of Nₘₐₓ normalized weights of quadrature for a product integral
    W::ArrayOfPhysicalTensors   # array of matrices holding quadrature weights
end
```

### Internal constructors:

There are two such constructors, i.e., the first constructor is
```
function VolterraIntegralScalarEquation(systemOfUnits::String,
                                        N::Integer,
                                        dt::PhysicalScalar,
                                        f₀::PhysicalScalar,
                                        g₀::PhysicalScalar,
                                        W::ArrayOfPhysicalTensors)
```
wherein `f₀` and `g₀` are initial values for the response variable and the control variable, respectively, while the second constructor is
```
function VolterraIntegralScalarEquation(dt::PhysicalScalar,
                                        N::Integer,
                                        Nₘₐₓ::Integer,
                                        n::MInteger, 
                                        f::ArrayOfPhysicalScalars,
                                        g::ArrayOfPhysicalScalars,
                                        t::ArrayOfPhysicalScalars,
                                        f′::ArrayOfPhysicalScalars,
                                        W::ArrayOfPhysicalTensors)
```
The former is used to initially construct such a data structure, while the latter is used by `JSON3.jl` whenever a data structure of this type is to be recreated from a JSON file.

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

After a Volterra integral's data structure has been created, thereby initializing a problem, one can advance a solution along its path from node *n-1* to node *n* by sequentially calling the following function:
```
function advance!(vie::VolterraIntegralScalarEquation, g′ₙ::ArrayOfPhysicalScalars, cₙ::PhysicalScalar)
```
wherein `vie` is an instance of type `VolterraIntegralScalarEquation.` Argument `g′ₙ` is an array of length 3 that contains rates for the control function at times $t_{n,1} = t_{n-1} + \mathrm{d}t/6$, $t_{n,2} = t_{n-1} + \mathrm{d}t/2$ and $t_{n,3} = t_{n-1} + 5\mathrm{d}t/6$, while argument `cₙ` is a scalar field.

There are applications where a solution is to be secured iteratively, e.g., during an optimization problem. In such cases one can call the following method as many times as needed before advancing a solution to its next node along its path via a call to `advance!.` Such iterative refinements result from a call to
```
function update!(vie::VolterraIntegralScalarEquation, g′ₙ::ArrayOfPhysicalScalars, cₙ::PhysicalScalar)
```
where array `g′ₙ` is of length 3 containing rates for the control function refined via an iterative step in some external global solver, e.g., a finite element solver. One need not call method `update!` whenever the coefficient `cₙ`  and control rate `g′ₙ` are known explicitly.

## For Vector Equations:

```
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
end
```

### Internal constructors:

There are two such constructors, i.e., the first constructor is
```
function VolterraIntegralVectorEquation(systemOfUnits::String,
                                        N::Integer,
                                        dt::PhysicalScalar,
                                        f₀::PhysicalVector,
                                        g₀::PhysicalVector,
                                        W::ArrayOfPhysicalTensors)
```
wherein `f₀` and `g₀` are initial values for the response variable and the control variable, respectively, while the second constructor is
```
function VolterraIntegralVectorEquation(dt::PhysicalScalar,
                                        N::Integer,
                                        Nₘₐₓ::Integer,
                                        n::MInteger, 
                                        f::ArrayOfPhysicalVectors,
                                        g::ArrayOfPhysicalVectors,
                                        t::ArrayOfPhysicalScalars,
                                        f′::ArrayOfPhysicalVectors,
                                        W::ArrayOfPhysicalTensors)
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
function advance!(vie::VolterraIntegralVectorEquation, g′ₙ::ArrayOfPhysicalVectors, cₙ::PhysicalScalar)
```
wherein `vie` is an instance of type `VolterraIntegralVectorEquation.` Argument `g′ₙ` is an array of length 3 that contains rates for the control function at times $t_{n,1} = t_{n-1} + \mathrm{d}t/6$, $t_{n,2} = t_{n-1} + \mathrm{d}t/2$ and $t_{n,3} = t_{n-1} + 5\mathrm{d}t/6$, while argument `cₙ` is a scalar field.

There are applications where a solution is to be secured iteratively, e.g., during an optimization problem. In such cases one can call the following method as many times as needed before advancing a solution to its next node along its path via a call to `advance!.` Such iterative refinements result from a call to
```
function update!(vie::VolterraIntegralVectorEquation, g′ₙ::ArrayOfPhysicalVectors, cₙ::PhysicalScalar)
```
where array `g′ₙ` is of length 3 containing rates for the control function refined via an iterative step in some external global solver, e.g., a finite element solver. One need not call method `update!` whenever the coefficient `cₙ`  and control rate `g′ₙ` are known explicitly.

## For Tensor Equations:

```
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
end
```

### Internal constructors:

There are two such constructors, i.e., the first constructor is
```
function VolterraIntegralTensorEquation(systemOfUnits::String,
                                        N::Integer,
                                        dt::PhysicalScalar,
                                        f₀::PhysicalTensor,
                                        g₀::PhysicalTensor,
                                        W::ArrayOfPhysicalTensors)
```
wherein `f₀` and `g₀` are initial values for the response variable and the control variable, respectively, while the second constructor is
```
function VolterraIntegralTensorEquation(dt::PhysicalScalar,
                                        N::Integer,
                                        Nₘₐₓ::Integer,
                                        n::MInteger, 
                                        f::ArrayOfPhysicalTensors,
                                        g::ArrayOfPhysicalTensors,
                                        t::ArrayOfPhysicalScalars,
                                        f′::ArrayOfPhysicalTensors,
                                        W::ArrayOfPhysicalTensors)
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
function advance!(vie::VolterraIntegralTensorEquation, g′ₙ::ArrayOfPhysicalTensors, cₙ::PhysicalScalar)
```
wherein `vie` is an instance of type `VolterraIntegralTensorEquation.`  Argument `g′ₙ` is an array of length 3 that contains rates for the control function at times $t_{n,1} = t_{n-1} + \mathrm{d}t/6$, $t_{n,2} = t_{n-1} + \mathrm{d}t/2$ and $t_{n,3} = t_{n-1} + 5\mathrm{d}t/6$, while argument `cₙ` is a scalar field.

There are applications where a solution is to be secured iteratively, e.g., during an optimization problem. In such cases one can call the following method as many times as needed before advancing a solution to its next node along its path via a call to `advance!.` Such iterative refinements result from a call to
```
function update!(vie::VolterraIntegralTensorEquation, g′ₙ::ArrayOfPhysicalTensors, cₙ::PhysicalScalar)
```
where array `g′ₙ` is of length 3 containing rates for the control function refined via an iterative step in some external global solver, e.g., a finite element solver. One need not call method `update!` whenever the coefficient `cₙ`  and control rate `g′ₙ` are known explicitly.


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

### Version 0.1.3

Bug fixes.

### Version 0.1.2

The coefficient c within a Volterra integral equation is now taken to be a function, assigned during construction. It is introduced via the solver, viz., as an argument in methods advance! and update!, instead of as a field in the data structure.

### Version 0.1.1

Rearranged order of fields in the three viscoelastic data structures so that they conform with the order of other data structures used for solvers in this family of software.

### Version 0.1.0

Initial release.