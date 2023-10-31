# Getting Started

## General Viewpoint

`Overall Task`: Solve *linear system of equations*

$$
\boldsymbol{Ax}=\boldsymbol{b}
$$

In the discussion afterwards, we assume that $\boldsymbol{A}$ is invertible, $\boldsymbol{A}\in \mathbb{R} ^{n\times n}$ but $n$ is very large.

In scientific computing, there are three important issues:

- ***Efficiency***
- ***Accuracy***
- ***Stability***: a concept measuring how sensitive the algorithm to perturbation

## Introductory Examples

Here are two examples for linear systems in scientific computing.

### ODE

Consider:

$$
\begin{cases}
	-u''\left( x \right) =f\left( x \right)\\
	u\left( 0 \right) =1, u\left( 1 \right) =2\\
\end{cases}, x\in \left( 0,1 \right) 
$$

where $f(x)$ is given. We want to know $u$ using the computer.

The important step is to do the discretization.

We split the interval $[0,1]$ into: $0=x_0<x_1<\cdots <x_{i-1}<x_i<x_{i+1}<\cdots <x_N=1$, where $x_{i+1}-x_i\triangleq h=\frac{1}{N}, \left( i=0,1,\cdots ,N-1 \right)$. We just need to know the solutions on these points for $u(x_i)$.

Let $u\left( x_i \right) \triangleq u_i$.

Using the finite difference method, we know that:

$$
u\prime\left( x \right) =\lim_{\varepsilon \rightarrow 0} \frac{u\left( x+\varepsilon \right) -u\left( x \right)}{\varepsilon}\approx \frac{u\left( x+\varepsilon \right) -u\left( x \right)}{\varepsilon}
$$ 

if $\varepsilon$ is small. We take $\varepsilon=h$, then $u\prime\left( x \right) \approx \frac{u\left( x+h \right) -u\left( x \right)}{h}$. We can get $\frac{u\left( x_i+h \right) -u\left( x_i \right)}{h}=\frac{u_{i+1}-u_i}{h}$.

Similarly, we can get $u''\left( x \right) \approx \frac{u_{i+1}-2u_i+u_{i-1}}{h^2}$. Also $u''\left( x \right) \propto O\left( h^2 \right)$.

Therefore, $-\frac{u_{i+1}-2u_i+u_{i-1}}{h^2}=f\left( x_i \right)$

To sum up, $\boldsymbol{u}=\left[ u_1,u_2,\cdots ,u_{N-1} \right]^T$ is a linear system that satisfies:

$$
\boldsymbol{A\cdot u}=\boldsymbol{f}
$$

We can get:

$$
\boldsymbol{u}=\left[ \begin{array}{c}
	u_1\\
	\vdots\\
	u_{N-1}\\
\end{array} \right] , \boldsymbol{A}=\left[ \begin{matrix}
	2&		-1&		0&		\cdots&		0\\
	-1&		2&		-1&		\cdots&		0\\
	0&		-1&		\ddots&		\ddots&		\vdots\\
	\vdots&		\vdots&		\ddots&		2&		-1\\
	0&		0&		\cdots&		-1&		2\\
\end{matrix} \right] , \boldsymbol{f}=\left[ \begin{array}{c}
	h^2f\left( x_1 \right) +?\\
	h^2f\left( x_2 \right)\\
	\vdots\\
	h^2f\left( x_N \right) +?\\
\end{array} \right] 
$$

### De-blurring Problem in Image Processing

The degrading image model is:

$$
\underset{\mathrm{observed} \ \mathrm{image}}{\underbrace{g\left( x \right) }}=\underset{\mathrm{true} \ \mathrm{image}}{\underbrace{u\left( x \right) }}+\int_{\varOmega}{\underset{\mathrm{blurring} \ \mathrm{kernel}}{\underbrace{k\left( x-y \right) }}u\left( y \right) \mathrm{d}y}+\underset{\mathrm{noise}}{\underbrace{\eta \left( x \right) }}
$$

Sometimes $k\left( x \right) =\frac{1}{c}e^{-\frac{x^2}{2}}$ which is Gaussian blur. Usually:

$$
k\left( x-y \right) =\begin{cases}
	c, if\,\,\left| x-y \right|<\varepsilon\\
	0, \mathrm{otherwise}\\
\end{cases}
$$

To simplify, we ignore the noise: let $\eta \left( x \right) =0$.

Given $g(x)$, we want to find $u(x)$. Assume that $k(x)$ is given and the image/signal is 1D. We can model the linear system as:

$$
\boldsymbol{u}=\left[ \begin{array}{c}
	u_1\\
	u_2\\
	\vdots\\
	u_n\\
\end{array} \right] , \boldsymbol{g}=\left[ \begin{array}{c}
	g_1\\
	g_2\\
	\vdots\\
	g_n\\
\end{array} \right] , \boldsymbol{A}=\left[ \begin{matrix}
	1&		0&		\cdots&		0\\
	0&		1&		\cdots&		0\\
	\vdots&		\vdots&		\ddots&		\vdots\\
	0&		0&		\cdots&		1\\
\end{matrix} \right] +\left[ \begin{matrix}
	k_0&		k_1&		k_2&		\cdots&		k_{n-1}\\
	k_1&		k_0&		k_1&		\cdots&		k_{n-2}\\
	k_2&		k_1&		k_0&		\cdots&		k_{n-3}\\
	\vdots&		\ddots&		\ddots&		\ddots&		\vdots\\
	k_{n-1}&		k_{n-2}&		k_{n-3}&		\cdots&		k_0\\
\end{matrix} \right] ,
\\
\boldsymbol{Au}=\boldsymbol{g}
$$

$\boldsymbol{A}$ is called [Toeplitz matrix](https://en.wikipedia.org/wiki/Toeplitz_matrix).
