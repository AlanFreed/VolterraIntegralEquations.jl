# VolterraIntegralEquations

This module provides a solver for Volterra integral equations of the second kind [11]. Specifically, this solver solves an integro-differential equation of the form
$$
\boldsymbol{f}^{\prime}(t) = \boldsymbol{g}^{\prime}(t) - \int_0^t k(t-τ) \, \boldsymbol{f}^{\prime}(τ) \, \mathrm{d}τ
$$
where the lower limit of integration is an initial time, taken to be 0, while the upper limit of integration is some current time $t$, with $\tau$ denoting a dummy variable of integration. For illustrative purposes, the independent variable in this document is taken to be time, but this is not necessary. 

This equation becomes a Volterra integral equation of the first kind whenever $\boldsymbol{g}^{\prime}$ is $\boldsymbol{0}$.

This Volterra integral equation describes an ordinary differential equation that, in turn, must be solved, i.e., it is an integro-differential equation. Here derivative $\boldsymbol{g}^{\prime} = \mathrm{d}\boldsymbol{g}/\mathrm{d}t$ is the known rate of some control function $\boldsymbol{g}$, while derivative $\boldsymbol{f}^{\prime} = \mathrm{d}\boldsymbol{f}/\mathrm{d}t$ is the unknown rate of a response function $\boldsymbol{f}$ to be solved for, whose evolution is characterized by a Volterra integral equation. Kernel $k(t-τ)$ is also a known quantity.

Functions $\boldsymbol{f}$ and $\boldsymbol{g}$ have the same physical units. They may be scalars, vectors or tensors in stature. Kernel $k$ is a memory function, which is the derivative of a relaxation function. It is scalar valued with physical units of $t^{-1}$, e.g., reciprocal time whenever $t$ has units of time.

To use this module, you will need to add the following private repositories to your project:

```
using Pkg
Pkg.add(url = "https://github.com/AlanFreed/PhysicalFields.jl")
Pkg.add(url = "https://github.com/AlanFreed/VolterraIntegralEquations.jl")
```

This software was written to accompany a book the authors are writing [6].

## Example: A Linear Viscoelastic Fiber

A viscoelastic fiber can be modeled as a Volterra integral equation of the second kind, viz., [5]
$$
\frac{\mathrm{d} \sigma}{\mathrm{d} t} =
	E \, \frac{\mathrm{d} \epsilon}{\mathrm{d} t} -
	\frac{E_0 - E_{\infty}}{E_{\infty}}
	\int_0^t K (t - \tau) \, \frac{\mathrm{d} \sigma}{\mathrm{d}\tau} \,
	\mathrm{d}\tau
$$
where $\sigma$ is stress and $\epsilon$ is strain, and $E_{\infty} \; (> 0)$ and $E_0 \; (> E_{\infty})$ are its rubbery and glassy moduli, respectively, while $K$ is a positive, monotonic-decreasing kernel that represents a memory function whose units are reciprocal time.  The first term on the right-hand side provides a glassy elastic change in stress that is attenuated by the second term, which introduces a viscous loss to this change in stress.

The elastic tangent modulus $E$ for a linear Hookean fiber comes from the compliance
$$
\frac{1}{E} = \frac{1}{E_0}
$$
wherein $E_0$ denotes a glassy modulus. While, for a nonlinear biologic fiber, its glassy-like elastic modulus associates with a compliance of [7]
$$
\frac{1}{E} = \frac{1}{E_0} + \frac{\beta + (\sigma - \sigma_r) / E_0 - \epsilon}{\beta E_r + 2(\sigma - \sigma_r)}
$$
wherein $\sigma_r$ denotes a residual stress, and $\beta$ designates a limiting state for an internal strain (where strain is caused by molecular reconfiguration).  The elastic tangent modulus associated with a fiber's strain-free reference configuration is designated as $E_r$ ($> 0$). Within the biologic fiber's linear region of its overall nonlinear response (where strain is caused by molecular stretching), moduli $E_\infty$ ($> E_r$) and $E_0$ ($> E_{\infty}$) denote its rubbery and glassy moduli, respectively.

### Volterra Functions

For a viscoelastic fiber, the unknown forcing function to be ascertained is a stress rate, viz.,
$$
f^{\prime}(t) = \frac{\mathrm{d} \sigma}{\mathrm{d} t}
$$
that is to be solved for in terms of:
1) a known kernel function
$$
k(t) = c \, K(t) \quad \text{wherein} \quad c = \frac{E_0 - E_{\infty}}{E_{\infty}}
$$
which is a modulus-scaled memory function, and 
2) a known control function
$$
g^{\prime}(t) = E \, \frac{\mathrm{d}\epsilon}{\mathrm{d}t}
$$
which is a modulus-scaled strain rate.

## Memory Functions

A selection of positive, monotonic-decreasing, viscoelastic, memory functions $K(t)$, with units of reciprocal time, have been preprogrammed into the software. Memory functions are the derivative of relaxation functions, the latter of which are more commonly found in the literature. See Freed [5,7] for a discussion of these functions.

1) **BOX**: the *𝑏𝑜𝑥* energy dissipation model of Neuber [10], a.k.a. Fung's **Q**uasi-**L**inear **V**iscoelastic (QLV) kernel, has a memory function of
$$
K(t) = \frac{\exp(-t/\tau_2) - \exp(-t/\tau_1)}{t \, \ln(\tau_2/\tau_1)}
$$
wherein $\tau_1$ and $\tau_2$ are characteristic times, ordered so that $0 < \tau_1 < \tau_2$. It is between these two times where dissipation is considered to occur. The BOX memory kernel is weakly singular.

2) **CCM**: **C**ole and **C**ole's [3,4] power-law **M**odel has a memory function of
$$
K(t) = \frac{\alpha}{t} \, \left( \frac{t}{\tau} \right)^{\alpha} 
\frac{1}{\bigl( 1 + (t / \tau)^{\alpha} \bigr)^2}
$$
wherein $\tau$ is a characteristic time and $\alpha \in (0,1]$ is the exponent of a power law. The CCM memory kernel is weakly singular whenever $\alpha \in (0,1)$. The CCM kernel equates with the SLS kernel at time 0 whenever $\alpha = 1$.

3) **FLS**: Caputo and Mainardi's [1,2] **F**ractional **L**inear **S**olid has a memory function of
$$
K(t) = -\frac{E_{\alpha,0} \left( - \left( t / \tau \right)^{\alpha} \right)}{t}
$$
wherein $\tau$ is a characteristic time, and $\alpha \in (0, 1]$ is a fractional order of evolution, with $E_{\alpha,\beta}(t)$ being the two-parameter Mittag-Leffler function. The FLS model contains the SLS model as a special case; specifically, they are equivalent whenever $\alpha = 1$. Mainardi's memory kernel becomes infinite at $K(0)$, provided $\alpha \ne 1$, i.e., the FLS memory kernel is weakly singular.

4) **KWW**: **K**ohlrausch's [8] and **W**illiams & **W**atts' [12] stretched exponential has a memory function of
$$
K(t) = \alpha \, \left( \frac{t}{\tau} \right)^{\alpha} \;
\frac{\exp \bigl( -(t/\tau)^{\alpha} \bigr)}{t}
$$
wherein $\tau$ is a characteristic time and $\alpha \in (0,1]$ is an exponent for the power of the argument in the exponential. The KWW model contains the SLS model as a special case; specifically, they are equivalent whenever $\alpha = 1$. The KWW memory kernel is weakly singular whenever $\alpha \in (0,1)$.

5) **MCM**: **M**axwell's **C**hain **M**odel, a.k.a. the Prony series model, has a memory function of
$$
K(t) = \sum_{i=1}^n (c_i / \tau_i) \, \exp(-t/\tau_i)
$$
whose coefficients $c_i$ are positive and sum as $\sum_{i=1}^n c_i = 1$, and whose characteristic times $\tau_i$ are ordered such that $0 < \tau_1 < \tau_2 < \cdots < \tau_n$. This is the popular Prony series that pervades the viscoelastic literature.

6) **MPL**: Williams' [13] **M**odified **P**ower-**L**aw model has a memory function of
$$
K(t) = \frac{\alpha}{\tau} \frac{1}{(1 + t / \tau)^{1+\alpha}}
$$
wherein $\tau$ is a characteristic time and $\alpha$ ($>0$) is in the exponent of a power law. The MPL kernel equates with the SLS kernel at time 0 whenever $\alpha = 1$.

7) **RFS**: Freed and Rajagopal's [7] **R**egularized **F**ractional linear **S**olid has a memory function of
$$
K(t) = \frac{-1}{E_{\alpha,1} \bigl( - ( \delta / \tau )^{\alpha} \bigr)} 
\frac{E_{\alpha,0} \bigl( - \bigl( ( \delta + t ) / \tau \bigr)^{\alpha} \bigr)}{\delta + t}
$$
wherein $\tau$ is a characteristic time, $\alpha \in (0, 1)$ is a fractional order of evolution, and $\delta$ is a regularization parameter satisfying the implicit equation
$$
\frac{\delta}{\tau} = -E_{\alpha,0} \left( - \left( 
\frac{\delta}{\tau} \right)^{\alpha}\right)
$$
which effectively shifts the singularity of an FLS kernel a short distance into negative time so the RFS model equates with the SLS model at time 0. As a consequence, the RFS memory kernel is finite valued throughout its domain of application.

8) **SLS**: Zener's [15] **S**tandard **L**inear **S**olid, a.k.a. the Maxwell-Debye kernel, has a memory function of
$$
K(t) = \frac{\mathrm{exp}(-t / \tau)}{\tau}
$$
wherein $\tau$ is a characteristic time.

The weakly singular memory kernels BOX, CCM and KWW could likely be regularized via a similar strategy laid out in [7] for regularizing the FLS kernel to become the RFS kernel, i.e., a time shift of $\delta$ could be introduced so that the weakly singular kernel evaluated at time $\delta$ equals the SLS kernel evaluated at time 0, i.e., $1/\tau$.

## Numerical Method

This solver implements a numerical method developed in an appendix of a book that is currently being written by the authors of this software [6]. It builds upon another numerical method developed by Young [14] which comes from the actuarial sciences literature.

Solutions to Volterra integrals of the second kind can be advanced through a block-by-block algorithmic strategy of
$$
\boldsymbol{f}^{\prime}(t_1) = \boldsymbol{g}^{\prime}(t_1) - \int_{0}^{t_1} k(t_1 - \tau) \, \boldsymbol{f}^{\prime}(\tau) \, \mathrm{d}\tau
$$
$$
\boldsymbol{f}^{\prime}(t_2) = \boldsymbol{g}^{\prime}(t_2) - \int_{0}^{t_2} k(t_2 - \tau) \, \boldsymbol{f}^{\prime}(\tau) \, \mathrm{d}\tau 
$$
$$
\vdots
$$
$$
\boldsymbol{f}^{\prime}(t_N) = \boldsymbol{g}^{\prime}(t_N) - \int_{0}^{t_N} k(t_N - \tau) \, \boldsymbol{f}^{\prime}(\tau) \, \mathrm{d}\tau
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
\boldsymbol{f}^{\prime}(t_2) & = \boldsymbol{g}^{\prime}(t_2) - \int_{0}^{t_1} k(t_2 - \tau) \, \boldsymbol{f}^{\prime}(\tau) \, \mathrm{d}\tau \\
& \phantom{= \boldsymbol{g}^{\prime}(t_3)} \; -
\int_{t_1}^{t_2} k(t_2 - \tau) \, \boldsymbol{f}^{\prime}(\tau) \, \mathrm{d}\tau
\end{aligned}
$$
$$
\begin{aligned}
\boldsymbol{f}^{\prime}(t_3) & = \boldsymbol{g}^{\prime}(t_3) - \int_{0}^{t_1} k(t_3 - \tau) \, \boldsymbol{f}^{\prime}(\tau) \, \mathrm{d}\tau \\
& \phantom{= \boldsymbol{g}^{\prime}(t_3)} \; -
		\int_{t_1}^{t_2} k(t_3 - \tau) \, \boldsymbol{f}^{\prime}(\tau) \, \mathrm{d}\tau \\
& \phantom{= \boldsymbol{g}^{\prime}(t_3)} \; -
		\int_{t_2}^{t_3} k(t_3 - \tau) \, \boldsymbol{f}^{\prime}(\tau) \, \mathrm{d}\tau
\end{aligned}
$$
$$
\vdots
$$
wherein, e.g., the product integral $\int_{0}^{t_1} k(t_2 - \tau) \, \boldsymbol{f}^{\prime}(\tau) \, \mathrm{d}\tau$ has a kernel $k$ whose fixed time $t_2$ lies outside the interval of integration, in this case $[0, t_1]$, over which time $\tau$ spans, while the product integral $\int_{t_1}^{t_2} k(t_2 - \tau) \, \boldsymbol{f}^{\prime}(\tau) \, \mathrm{d}\tau$ has a kernel $k$ whose fixed time, in this case $t_2$, is now the upper limit of integration. These two integrals will have different rules of quadratures representing them.

### Quadratures

There are product integrals of the Fredholm type, viz.,
$$
\int_{t_l}^{t_{l+1}} k(t_{n+1} - \tau) \, \boldsymbol{f}^{\prime}(\tau) \, \mathrm{d}\tau =
	\sum_{j=1}^J W_j (k, \mathrm{d}t) \, \boldsymbol{f}^{\prime}(t_l + t_j)  + \boldsymbol{\varepsilon}
$$
where $l = 0,1,2, \ldots, n \! - \! 1,$ with $t_j \in [0, \mathrm{d}t]$ representing the $J$ local quadrature nodes for a selected integration rule, recalling that $\mathrm{d}t = t_{l+1} - t_l$. Here $W_j$ denotes a weight of quadrature for which $\boldsymbol{\varepsilon}$ is its local truncation error.

There are also product integrals of the Volterra type, viz.,
$$
\int_{t_n}^{t_{n+1}} k(t_{n+1} - \tau) \, \boldsymbol{f}^{\prime}(\tau) \, \mathrm{d} \tau = \sum_{j=1}^J W_j (k, \mathrm{d}t) \, \boldsymbol{f}^{\prime}(t_n + t_j) + \boldsymbol{\varepsilon}
$$
where the upper limit of integration $t_{n+1}$ now appears as an argument in the kernel function $k$.

These two quadrature rules result in similar, but different, weights of quadrature between them.

In practice, it is considered that the global interval of integration $\mathrm{d}t$ is sufficiently small so that its local truncation error $\boldsymbol{\varepsilon}$ can be safely neglected. Even so, the user ought to be cognizant of this consideration whenever one is assessing the validity of one's solution.

### Block-by-Block Solution Strategy

Solutions to Volterra integral equations of the second kind take on the form of a sequentially solved linear equation
$$
\boldsymbol{f}^{\prime}_{n} = \Bigl( \boldsymbol{I} + \boldsymbol{W}^{\mathsf{T}}_{\!1} \Bigr)^{-1} \left( \boldsymbol{g}^{\prime}_n - \sum_{m=1}^{n-1} \boldsymbol{W}^{\mathsf{T}}_{\!n-m+1} \boldsymbol{f}^{\prime}_{m} \right) ,
\qquad n = 1, 2, \ldots, N
$$
where $N$ specifies the number of global integration steps, called nodes, that are required to traverse a solution's path, whose weights of quadrature are $\boldsymbol{W}_n$, with weight $\boldsymbol{W}_1$ being distinct in form from all the others. In order for this algorithm to work, it is necessary that the matrix $\boldsymbol{I} + \boldsymbol{W}_1^{\mathsf{T}}$ not be singular so that its inverse exists. 

For the method implemented here, there are three local nodes of integration per single global node, i.e., $J = 3$, and as such, the control function $\boldsymbol{g}^{\prime}_n$ and the forcing function $\boldsymbol{f}^{\prime}_n$ are both vectors of length 3, while the weights of quadrature $\boldsymbol{W}_n$ are $3 \times 3$ matrices. Specifically, a mid-point quadrature rule is selected wherein
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
\boldsymbol{f}^{\prime}(t_n \! + \! \tfrac{1}{6} \mathrm{d}t) \\
\boldsymbol{f}^{\prime}(t_n \! + \! \tfrac{1}{2} \mathrm{d}t) \\
\boldsymbol{f}^{\prime}(t_n \! + \! \tfrac{5}{6} \mathrm{d}t) \\
\end{matrix} \right\}
\quad \text{and} \quad
\boldsymbol{g}^{\prime}_n = \left\{ \begin{matrix}
\boldsymbol{g}^{\prime}(t_{n,1}) \\
\boldsymbol{g}^{\prime}(t_{n,2}) \\
\boldsymbol{g}^{\prime}(t_{n,3}) \\
\end{matrix} \right\}
= \left\{ \begin{matrix}
\boldsymbol{g}^{\prime}(t_n \! + \! \tfrac{1}{6} \mathrm{d}t) \\
\boldsymbol{g}^{\prime}(t_n \! + \! \tfrac{1}{2} \mathrm{d}t) \\
\boldsymbol{g}^{\prime}(t_n \! + \! \tfrac{5}{6} \mathrm{d}t) \\
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
where Young [3] calls $\boldsymbol{X}$ the alternant matrix and $\boldsymbol{\mu}_n$ the $n^{\text{th}}$ moment matrix, which are consequences of expressing $\boldsymbol{f}^{\prime}$ in a Taylor series expanded about the midpoint to its local span of integration.

The global nodes of integration associate with times $t_n$ where $t_n = t_{n-1} + \mathrm{d}t \; \forall \; n$ given $t_0 = 0$, with $\mathrm{d}t$ being a distance separating neighboring global nodes. The local nodes of integration associate with times $t_{n,j}$ where $j = 1,2,3$ and where $t_{n,1} = t_n + \tfrac{1}{6} \mathrm{d}t$, $t_{n,2} = t_n + \tfrac{1}{2} \mathrm{d}t$ and $t_{n,3} = t_n + \tfrac{5}{6} \mathrm{d}t$, with a distance of $\tfrac{1}{3} \mathrm{d}t$ separating neighboring local nodes. These local nodes do not contain any global nodes. Consequently, this solution strategy can, in principle, be applied to kernels that are singular at the upper limit of integration (like memory kernels: BOX, CCM, FLS and KWW).

### Moment Matrices

The greatest expense in implementing this numerical method is in computing its moment matrices. Fortunately, once gotten they can be reused in future solutions. There is a more efficient algorithm for implementing an exponential kernel that is based upon its recursive property [5], but that algorithm is restricted to that kernel alone. The algorithm presented below is applicable to all types of kernel functions, and is therefore versatile.

The moment matrices are solutions to integral equations. For $n=1$, moment $\boldsymbol{\mu}_1$ describes a $3 \times 3$ matrix whose elements are solutions to the integral equation
$$
\mu_{ij,1} = \frac{1}{h^{i-1}} \int_0^{(j - 1/2)h} \bigl( \tau - \tfrac{1}{2} (j - 1/2) h \bigr)^{i-1} \, k \bigl( (j - 1/2) h - \tau \bigr) \, \mathrm{d} \tau
$$
where $h = \tfrac{1}{3} \mathrm{d}t$ and $i,j=1,2,3$. The remaining moment matrices $\boldsymbol{\mu}_n$, where $n=2, 3, \ldots, N$, describe $3 \times 3$ matrices whose elements are solutions to the integral equation
$$
\mu_{ij,n} = \frac{1}{h^{i-1}} \int_0^{\mathrm{d}t} \bigl( \tau - \tfrac{1}{2} \mathrm{d}t \bigr)^{i-1} \, k \bigl( (n-1) \mathrm{d}t + (j - 1/2) h - \tau \bigr) \, \mathrm{d} \tau
$$
These are integrals of the kernel function $k$ scaled by moment arms measured from the midpoints of their spans of integration. These moment arms are raised to powers of $i \! - \! 1$. They arise because of the Taylor series expansion imposed.

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
\mu_{ij,1} = \frac{(j \! - \! 1/2) \, \mathrm{d}t}{6} \sum_{s=1}^3 w_s m_{ij,s} k \left( \tfrac{1}{6} (j \! - \! 1/2) (1 \! - \! x_s)^{\vphantom{|}} \mathrm{d}t \right) , \quad i,j = 1,2,3
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
\mu_{ij,n} = \frac{\mathrm{d}t}{2} \sum_{s=1}^3 w_s v_{i,s} k \left( \bigl( n - \tfrac{1}{3} (5 - j) - \tfrac{1}{2} x_s \bigr)^{\vphantom{|}} \mathrm{d}t \right) , \quad i,j = 1,2,3
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

For those kernels $k$ that are monotonic-decreasing functions, a number $N_{\max}$ exists beyond which point memory of the past effectively fades away. How one assigns $N_{\max}$ will depend upon the kernel $k$ (specifically, its characteristic time), the global step size $\mathrm{d}t$, and the accuracy sought in a solution. Considering that such an $N_{\max}$ exists, then, for integration steps where $n \le N_{\max}$, a solution $\boldsymbol{f}^{\prime}_n$ advances along its path according to the linear equation
$$
\boldsymbol{f}^{\prime}_n = \Bigl( \boldsymbol{I} + \boldsymbol{W}^{\mathsf{T}}_1 \Bigr)^{-1} \left( \boldsymbol{g}^{\prime}_n - \sum_{m=1}^{n-1} \boldsymbol{W}^{\mathsf{T}}_{n-m+1} \boldsymbol{f}^{\prime}_{m} \right)
$$
that when $n > N_{\max}$ it continues to advance along its path, but now according to the linear equation
$$
\boldsymbol{f}^{\prime}_n = \Bigl( \boldsymbol{I} + \boldsymbol{W}^{\mathsf{T}}_1 \Bigr)^{-1} \left( \boldsymbol{g}^{\prime}_n - \sum_{m=1}^{N_{\max}-1} \boldsymbol{W}^{\mathsf{T}}_{N_{\max}-m+1} \boldsymbol{f}^{\prime}_{m+n-N_{\max}} \right)
$$
which removes forcing functions from $\boldsymbol{f}^{\prime}_1$ through $\boldsymbol{f}^{\prime}_{n-N_{\max}}$ from its summation history. Their memories are so distant that they have effectively been forgotten.

## Solving the Resulting Differential Equations

Here a Volterra integral equation of the second kind, whose method of solution has just been laid out, describes a differential equation whose solution is sought. Specifically, the control $\boldsymbol{g}^{\prime}$ and response $\boldsymbol{f}^{\prime}$ function rates can both be integrated using a Newton-Cotes formula of the form
$$
\boldsymbol{g}_{n+1} = \boldsymbol{g}_n + \frac{\mathrm{d}t}{8} \left( 3 \boldsymbol{g}^{\prime} (t_{n,1}) + 2 \boldsymbol{g}^{\prime} (t_{n,2})+ 3 \boldsymbol{g}^{\prime} (t_{n,3}) \right)
$$
and
$$
\boldsymbol{f}_{n+1} = \boldsymbol{f}_n + \frac{\mathrm{d}t}{8} \left( 3 \boldsymbol{f}^{\prime} (t_{n,1}) + 2 \boldsymbol{f}^{\prime} (t_{n,2})+ 3 \boldsymbol{f}^{\prime} (t_{n,3}) \right)
$$
with the $\boldsymbol{g}^{\prime}$ being known and the $\boldsymbol{f}^{\prime}$ being solutions to a Volterra integral equation. Time $t_{n+1}$ associates with the next global node of integration, while times $t_{n,j}$, $j=1,2,3$ associate with the three local nodes of integration spanning $[t_n , t_{n+1}]$ whereat the $\boldsymbol{f}^{\prime}_{n,j} = \boldsymbol{f}^{\prime}(t_{n,j})$ have been determined.

There are applications where the control function may depend upon either or both the control and/or response functions, e.g., the biologic fiber model in the viscoelastic example described above, in which case one will also need solutions to be gotten at the three local nodes of integration, too, viz.,
$$
\boldsymbol{f}_{n,1} = \boldsymbol{f}_n + \frac{\mathrm{d}t}{72} \left(
17 \boldsymbol{f}^{\prime} (t_{n,1}) -7 \boldsymbol{f}^{\prime} (t_{n,2})+ 2 \boldsymbol{f}^{\prime} (t_{n,3}) \right)
$$
$$
\boldsymbol{f}_{n,2} = \boldsymbol{f}_n + \frac{\mathrm{d}t}{8} \left(
3 \boldsymbol{f}^{\prime} (t_{n,1}) + \boldsymbol{f}^{\prime} (t_{n,2}) \right)
$$
$$
\boldsymbol{f}_{n,3} = \boldsymbol{f}_n + \frac{5 \, \mathrm{d}t}{72} \left( 5 \boldsymbol{f}^{\prime} (t_{n,1}) + 5 \boldsymbol{f}^{\prime} (t_{n,2})+ 2 \boldsymbol{f}^{\prime} (t_{n,3}) \right)
$$
with like expressions establishing $\boldsymbol{g}_{n,j}$. The method of Newton-Cotes provides these formulæ, too.

## Software

This software uses the `PhysicalFields.jl` package.

### Memory Functions

Normalized memory functions, i.e., the $K(t)$ in $k(t) = (E_0 - E_{\infty}) \, K(t) / E_{\infty}$, like those introduced above are callable in software through the following functions. For those kernels that are not weakly singular
```
function MCM(systemOfUnits::String, time::PhysicalScalar, parameters::Tuple)::PhysicalScalar
function MPL(systemOfUnits::String, time::PhysicalScalar, parameters::Tuple)::PhysicalScalar
function RFS(systemOfUnits::String, time::PhysicalScalar, parameters::Tuple)::PhysicalScalar
function SLS(systemOfUnits::String, time::PhysicalScalar, parameters::Tuple)::PhysicalScalar
```
and for those kernels that are weakly singular
```
function BOX(systemOfUnits::String, time::PhysicalScalar, parameters::Tuple)::PhysicalScalar
function CCM(systemOfUnits::String, time::PhysicalScalar, parameters::Tuple)::PhysicalScalar
function FLS(systemOfUnits::String, time::PhysicalScalar, parameters::Tuple)::PhysicalScalar
function KWW(systemOfUnits::String, time::PhysicalScalar, parameters::Tuple)::PhysicalScalar
```
while other memory functions of one's own design can be introduced, too, provided they have an interface of:
```
function <myMemoryFunction>(systemOfUnits::String, time::PhysicalScalar, parameters::Tuple)::PhysicalScalar
```
Herein string `systemOfUnits` is (at present) either "SI" or "CGS", with `time` specifying the argument of `K(t)` sent to the memory function for evaluation, and whose `parameters,` or material constants, are supplied via a Tuple. All memory functions return an instance of type `PhysicalFields.PhysicalScalar.`

The tuples that are to be supplied for the incorporated memory functions, i.e. the model's material constants, are:

1) BOX: parameters = $(τ_1, τ_2)$
2) CCM: parameters = $(α, τ)$
3) FLS: parameters = $(α, τ)$
4) KWW: parameters = $(α, τ)$
5) MCM: parameters = $(c_1, c_2, …, c_n, τ_1, τ_2, …, τ_n)$
6) MPL: parameters = $(α, τ)$
7) RFS: parameters = $(α, δ, τ)$
8) SLS: parameters = $(τ,)$

### Constructing Weights of Quadrature for Volterra Integral Equations

The weights of quadrature for this, a solver of Young's [14] method, can be created with a call to the following function
```
function normalizedQuadratureWeights(K::Function,
                                     systemOfUnits::String,
                                     N::Integer,
                                     dTime::PhysicalScalar, 
                                     parameters::Tuple)::ArrayOfPhysicalTensors

```
where `K` is any of the eight memory functions addressed above, or one of your own design. This function is called internally whose arguments include `systemOfUnits` and `parameters,` with its argument for time being constructed from `dTime` and `N.` Argument `dTime` denotes the global time step of the solver, with there being `N` solution nodes uniformly spaced over the solution path. These quadrature weights are stored in an instance of type `PhysicalFields.ArrayOfPhysicalTensors,` and as such, can be stored to a JSON file for a later retrieval.

### Solver for Volterra Integral Equations of the Second Kind

The exported solvers for Volterra integral equations are implementations of the abstract class
```
abstract type VolterraIntegralEquation end
```
and they come in three flavors: for scalar equations, for vector equations, and for tensor equations.

#### For scalar equations:
```
struct VolterraIntegralScalarEquation <: VolterraIntegralEquation
    # Dimensioning fields
    dt::PhysicalScalar          # distance separating neighboring global integration nodes
    N::Integer                  # number of integration nodes in a solution path
    Nₘₐₓ::Integer               # maximum number of nodes whose history is kept
    n::MInteger                 # current node along a solution path
    # Arrays of length N+1 holding integrated variable rates at the global nodes
    f::ArrayOfPhysicalScalars   # array of integrated response function values
    g::ArrayOfPhysicalScalars   # array of integrated control function values
    t::ArrayOfPhysicalScalars   # array of times, the independent variable
    # Array of length 3N holding response function rates at the local nodes
    f′::ArrayOfPhysicalScalars  # array of response function rates
    # Coefficient that scales the normalized weights of quadrature
    c::PhysicalScalar           # e.g., c = (E₀ - E∞)/E∞ in viscoelaticity
    # Array of Nₘₐₓ normalized weights of quadrature for a product integral
    W::ArrayOfPhysicalTensors   # array of matrices holding the quadrature weights
end
```
which has two constructors, i.e., the first constructor is
```
function VolterraIntegralScalarEquation(systemOfUnits::String,
                                        N::Integer,
                                        dt::PhysicalScalar,
                                        f₀::PhysicalScalar,
                                        g₀::PhysicalScalar,
                                        c::PhysicalScalar,
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
                                        c::PhysicalScalar,
                                        W::ArrayOfPhysicalTensors)
```
The former is used to initially construct such a data structure, while the latter is used by `JSON3.jl` whenever a data structure of this type is to be recreated from a JSON file.

#### For vector equations:
```
struct VolterraIntegralVectorEquation <: VolterraIntegralEquation
    # Dimensioning fields
    dt::PhysicalScalar          # distance separating neighboring global integration nodes
    N::Integer                  # number of integration nodes in a solution path
    Nₘₐₓ::Integer               # maximum number of nodes whose history is kept
    n::MInteger                 # current node along a solution path
    # Arrays of length N+1 holding integrated variable rates at the global nodes
    f::ArrayOfPhysicalVectors   # array of integrated response function values
    g::ArrayOfPhysicalVectors   # array of integrated control function values
    t::ArrayOfPhysicalScalars   # array of times, the independent variable
    # Array of length 3N holding response function rates at the local nodes
    f′::ArrayOfPhysicalVectors  # array of response function rates
    # Coefficient that scales the normalized weights of quadrature
    c::PhysicalScalar           # e.g., c = (E₀ - E∞)/E∞ in viscoelaticity
    # Array of Nₘₐₓ normalized weights of quadrature for a product integral
    w::ArrayOfPhysicalTensors   # array of matrices holding the quadrature weights
end
```
which has two constructors, i.e., the first constructor is
```
function VolterraIntegralVectorEquation(systemOfUnits::String, 
                                        N::Integer,
                                        dt::PhysicalScalar,
                                        f₀::PhysicalVector,
                                        g₀::PhysicalVector,
                                        c::PhysicalScalar,
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
                                        c::PhysicalScalar,
                                        W::ArrayOfPhysicalTensors)
```
The former is used to initially construct such a data structure, while the latter is used by `JSON3.jl` whenever a data structure of this type is to be recreated from a JSON file.

#### For tensor equations:
```
struct VolterraIntegralTensorEquation <: VolterraIntegralEquation
    # Dimensioning fields
    dt::PhysicalScalar          # distance separating neighboring global integration nodes
    N::Integer                  # number of integration nodes in a solution path
    Nₘₐₓ::Integer               # maximum number of nodes whose history is kept
    n::MInteger                 # current node along a solution path
    # Arrays of length N+1 holding integrated variable rates at the global nodes
    f::ArrayOfPhysicalTensors   # array of integrated response function values
    g::ArrayOfPhysicalTensors   # array of integrated control function values
    t::ArrayOfPhysicalScalars   # array of times, the independent variable
    # Array of length 3N holding response function rates at the local nodes
    f′::ArrayOfPhysicalTensors  # array of response function rates
    # Coefficient that scales the normalized weights of quadrature
    c::PhysicalScalar           # e.g., c = (E₀ - E∞)/E∞ in viscoelaticity
    # Array of Nₘₐₓ normalized weights of quadrature for a product integral
    w::ArrayOfPhysicalTensors   # array of matrices holding the quadrature weights
end
```
which has two constructors, i.e., the first constructor is
```
function VolterraIntegralTensorEquation(systemOfUnits::String, 
                                        N::Integer,
                                        dt::PhysicalScalar,
                                        f₀::PhysicalTensor,
                                        g₀::PhysicalTensor,
                                        c::PhysicalScalar,
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
                                        c::PhysicalScalar,
                                        W::ArrayOfPhysicalTensors)
```
The former is used to initially construct such a data structure, while the latter is used by `JSON3.jl` whenever a data structure of this type is to be recreated from a JSON file.

### Methods

#### Working with JSON files.

These objects are persistent in the sense that they can be stored-to and retrieved-from a file; specifically, one can call the following method to write an object to a JSON file.
```
function toFile(y::VolterraIntegralScalarEquation, json_stream::IOStream)
function toFile(y::VolterraIntegralVectorEquation, json_stream::IOStream)
function toFile(y::VolterraIntegralTensorEquation, json_stream::IOStream)
```
Or one can call the following method to read an object from a JSON file.
```
function fromFile(::Type{VolterraIntegralScalarEquation}, json_stream::IOStream)::VolterraIntegralScalarEquation
function fromFile(::Type{VolterraIntegralVectorEquation}, json_stream::IOStream)::VolterraIntegralVectorEquation
function fromFile(::Type{VolterraIntegralTensorEquation}, json_stream::IOStream)::VolterraIntegralTensorEquation
```
where an appropriate `json_stream` comes from functions
```
function PhysicalFields.openJSONReader(my_dir_path::String, my_file_name::String)::IOStream
function PhysicalFields.openJSONWriter(my_dir_path::String, my_file_name::String)::IOStream
function PhysicalFields.closeJSONStream(json_stream::IOStream)
```

#### Calling the solver to advance a solution.

After a Volterra integral's data structure has been created, thereby initializing a problem, one can advance a solution along its path from node *n* to node *n+1* by sequentially calling the following function:
```
function advance!(VIE::VolterraIntegralEquation, g′ₙ::Tuple)
```
wherein `VIE` is an instance of either `VolterraIntegralScalarEquation`, `VolterraIntegralVectorEquation` or `VolterraIntegralTensorEquation,` while tuple `g′ₙ` contains rates of the control function at times $t_{n,1} = t_n + \mathrm{d}t/6$, $t_{n,2} = t_n + \mathrm{d}t/2$ and $t_{n,3} = t_n + 5\mathrm{d}/6$, viz., `g′ₙ` = $\bigl( g^{\prime} (t_{n,1}), g^{\prime} (t_{n,2}), g^{\prime} (t_{n,3}) \bigr)$, which are scalar, vector or tensor valued, as appropriate for the `VIE` object supplied.

There are applications where a solution is to be secured iteratively, e.g., during an optimization problem. In such cases one can call the following method as many times as needed before advancing a solution to its next node along its path via a call to `advance!.` Such iterative refinements result from a call to
```
function update!(VIE::VolterraIntegralEquation, g′ₙ::Tuple)
```
where the tuple `g′ₙ` of control rates is to be refined via an iterative step in a more global solver than this one. One need not call method `update!` whenever the control rates `g′ₙ` are known explicitly.  

## References

1) Caputo, M. and Mainardi, F., "Linear models of dissipation in anelastic solids," *Rivista del Nuoro Cimento*, **1** (1971), 161-198.

2) Caputo, M. and Mainardi, F., "A new dissipation model based on memory mechanism," *Pure and Applied Geophysics*, **91** (1971), 134-147.

3) Cole, K.S. and Cole, R.H., "Dispersion and absorption in dielectrics I. Alternating current characteristics," *Journal of Chemical Physics*, **9** (1941), 342-351.

4) Cole, K.S. and Cole, R.H., "Dispersion and absorption in dielectrics II. Direct current characteristics," *Journal of Chemical Physics*, **10** (1942), 98-105.

5) Freed, A.D., *Soft Solids: A primer to the theoretical mechanics of materials*, Modeling and Simulation in Science, Engineering and Technology. Basel: Birkhäuser, 2014.

6) Freed, A.D. and Clayton, J.D., *Application of Laplace Stretch in Alveolar Mechanics: A Case Study of Blast and Blunt Trauma*, in preparation.

7) Freed, A.D. and Rajagopal, K.R. "A viscoelastic model for describing the response of biological fibers," *ACTA Mechanica*, **227** (2016), 3367-3380.

8) Kohlrausch, R., "Ueber das Dellmann'sche Elektrometer," *Annalen der Physik und Chemie*, **72** (1847), 353-405.

9) Maxwell, J.C., "On the dynamical theory of gases," *Philosophical Transactions of the Royal Society, London*, **157** (1867), 49-88.

10) Neubert, H.K., "A simple model representing internal damping in solid materials," *The Aeronautical Quarterly*, **14** (1963), 187-210.

11) Volterra, V. *Theory of functionals and of integral and integro-differential equations*. Glasgow: Blackie and Son, 1930.

12) Williams, G. and Watts, D.C., "Non-symmetrical dielectric relaxation behaviour arising from a simple empirical decay function," *Transactions of the Faraday Society*, **66** (1970), 80-85.

13) Williams, M.L., "Structural analysis of viscoelastic materials," *AIAA Journal*, **2** (1964), 785-808.

14) Young, A., "Approximate product integration," *Proceedings of the Royal Society, London*, **A-224** (1954), 552-561.

15) Zener, C., *Elasticity and Anelasticity of Metals*. Chicago: University of Chicago Press, 1948.

## Version History

### Version 0.1.1

Rearranged order of fields in the three viscoelastic data structures so that they conform with the order of other data structures used for solvers in this family of software.

### Version 0.1.0

Initial release.