# Conditioning

## Introduction

Let's start with an example:

Given $x\in \mathbb{R} ^+$ , compute $\sqrt{x}$ . We are considering the problem not the algorithm here.

We denote the problem as:

$$
\delta f=f\left( x+\delta x \right) -f\left( x \right) :\begin{cases}
	x\longmapsto \sqrt{x}=f\left( x \right)\\
	\Downarrow\\
	x+\delta x\longmapsto \sqrt{x+\delta x}=f\left( x+\delta x \right)\\
\end{cases}
$$

Rather than considering $\frac{\left\| \delta f \right\|}{\left\| fx \right\|}$ , we can consider the relative ratio:

$$
\frac{\frac{\left\| \delta f \right\|}{\left\| f \right\|}}{\frac{\left\| \delta x \right\|}{\left\| x \right\|}}
$$

For $f\left( x \right) =\sqrt{x}$ :

$$
\frac{\frac{\left\| \delta f \right\|}{\left\| f \right\|}}{\frac{\left\| \delta x \right\|}{\left\| x \right\|}}=\frac{\frac{\left| \sqrt{x+\delta x}-\sqrt{x} \right|}{\sqrt{x}}}{\frac{\left| \delta x \right|}{\left| x \right|}}=\frac{\frac{\left| x+\delta x-x \right|}{\sqrt{x}\left( \sqrt{x+\delta x}+\sqrt{x} \right)}}{\frac{\left| \delta x \right|}{\left| x \right|}}
\\
=\frac{\left| x \right|}{\sqrt{x}\left( \sqrt{x+\delta x}+\sqrt{x} \right)}\approx \frac{\left| x \right|}{\sqrt{x}\left( 2\sqrt{x} \right)}=\frac{1}{2}
$$

It means that the output does not amplify the permutation/mistake of the input.

## Definition

Assume $\delta x$ is small, we define the **condition number** as:

$$
\lim_{\alpha \rightarrow 0} \underset{\left\| \delta x \right\| \leqslant \alpha}{\mathrm{sup}}\frac{\frac{\left\| \delta f \right\|}{\left\| f \right\|}}{\frac{\left\| \delta x \right\|}{\left\| x \right\|}}
$$

Here $f$ is any problem.

## Matrix Conditioning Number

### Formula of Matrix Conditioning Number

>If $\boldsymbol{A}\in \mathbb{R} ^{n\times n}$ , then **matrix condition number** can be defined as:
>
>$$
\kappa \left( \boldsymbol{A} \right) =\left\| \boldsymbol{A} \right\| \cdot \left\| \boldsymbol{A} \right\| ^{-1}
$$

Why is this definition important?

### Proof of the Formula

Consider $\boldsymbol{Ax}=\boldsymbol{b}$, perturbations may occur to $\boldsymbol{A}, \boldsymbol{x}, \boldsymbol{b}$ .

#### Case 1

$\boldsymbol{x}\rightarrow \boldsymbol{x}+\delta \boldsymbol{x}$ , $\boldsymbol{A}$ is given.

$$
\boldsymbol{Ax}=\boldsymbol{b},
\\
\boldsymbol{A}\left( \boldsymbol{x}+\delta \boldsymbol{x} \right) =\boldsymbol{b}+\delta \boldsymbol{b}, \delta \boldsymbol{b}=\boldsymbol{A}\delta \boldsymbol{x},
$$

$$
\underset{\left\| \delta \boldsymbol{x} \right\| \leqslant \alpha}{\mathrm{sup}}\frac{\frac{\left\| \boldsymbol{A}\left( \boldsymbol{x}+\delta \boldsymbol{x} \right) -\boldsymbol{Ax} \right\|}{\left\| \boldsymbol{Ax} \right\|}}{\frac{\left\| \delta \boldsymbol{x} \right\|}{\left\| \boldsymbol{x} \right\|}}=\underset{\left\| \delta \boldsymbol{x} \right\| \leqslant \alpha}{\mathrm{sup}}\frac{\frac{\left\| \boldsymbol{A}\delta \boldsymbol{x} \right\|}{\left\| \boldsymbol{Ax} \right\|}}{\frac{\left\| \delta \boldsymbol{x} \right\|}{\left\| \boldsymbol{x} \right\|}}=\underset{\left\| \delta \boldsymbol{x} \right\| \leqslant \alpha}{\mathrm{sup}}\frac{\left\| \boldsymbol{A}\delta \boldsymbol{x} \right\|}{\left\| \delta \boldsymbol{x} \right\|}\cdot \frac{\left\| \boldsymbol{x} \right\|}{\left\| \boldsymbol{Ax} \right\|}
\\
=\frac{\left\| \boldsymbol{x} \right\|}{\left\| \boldsymbol{Ax} \right\|}\cdot \underset{\left\| \delta \boldsymbol{x} \right\| \leqslant \alpha}{\mathrm{sup}}\frac{\left\| \boldsymbol{A}\delta \boldsymbol{x} \right\|}{\left\| \delta \boldsymbol{x} \right\|}=\frac{\left\| \boldsymbol{x} \right\|}{\left\| \boldsymbol{Ax} \right\|}\cdot \left\| \boldsymbol{A} \right\| 
$$

Introduce $\boldsymbol{y}=\boldsymbol{Ax}$ , $\boldsymbol{A}$ is invertible, then:

$$
\boldsymbol{x}=\boldsymbol{A}^{-1}\boldsymbol{y},
$$

$$
\frac{\left\| \boldsymbol{x} \right\|}{\left\| \boldsymbol{Ax} \right\|}\left\| \boldsymbol{A} \right\| =\frac{\left\| \boldsymbol{A}^{-1}\boldsymbol{y} \right\|}{\left\| \boldsymbol{y} \right\|}\left\| \boldsymbol{A} \right\| \underset{\mathrm{achievable}}{\underbrace{\leqslant }}\left\| \boldsymbol{A}^{-1} \right\| \cdot \left\| \boldsymbol{A} \right\| 
$$

#### Case 2

Perturbation happens to $\boldsymbol{b}$ :

$$
\boldsymbol{Ax}=\boldsymbol{b}, \boldsymbol{x}=\boldsymbol{A}^{-1}\boldsymbol{b},
\\
\boldsymbol{b}\rightarrow \boldsymbol{b}+\delta \boldsymbol{b}
$$

Using case 1, we can easily have the condition number:

$$
\left\| \left( \boldsymbol{A}^{-1} \right) ^{-1} \right\| \cdot \left\| \boldsymbol{A}^{-1} \right\| =\left\| \boldsymbol{A} \right\| \cdot \left\| \boldsymbol{A}^{-1} \right\| 
$$

#### Case 3

Given $\boldsymbol{A}, \boldsymbol{b}$ , and let $\boldsymbol{A}\rightarrow \boldsymbol{A}+\delta \boldsymbol{A}$ . Then it incurs $\boldsymbol{x}\rightarrow \boldsymbol{x}+\delta \boldsymbol{x}$ . We can get:

$$
\left( \boldsymbol{A}+\delta \boldsymbol{A} \right) \left( \boldsymbol{x}+\delta \boldsymbol{x} \right) =\boldsymbol{b},
\\
\Longrightarrow \boldsymbol{Ax}+\boldsymbol{A}\cdot \delta \boldsymbol{x}+\delta \boldsymbol{A}\cdot \boldsymbol{x}+\delta \boldsymbol{A}\cdot \delta \boldsymbol{x}=\boldsymbol{b}
$$

Usually $\delta \boldsymbol{A}\cdot \delta \boldsymbol{x}$ is smaller than $\delta \boldsymbol{x}$ or $\delta \boldsymbol{A}$, then we can also get:

$$
\boldsymbol{A}\cdot \delta \boldsymbol{x}+\delta \boldsymbol{A}\cdot \boldsymbol{x}=0,
\\
\boldsymbol{A}\cdot \delta \boldsymbol{x}=-\delta \boldsymbol{A}\cdot \boldsymbol{x},
\\
\delta \boldsymbol{x}=\boldsymbol{A}^{-1}\left( -\delta \boldsymbol{A} \right) \boldsymbol{x},
\\
\left\| \delta \boldsymbol{x} \right\| \leqslant \left\| \boldsymbol{A}^{-1} \right\| \left\| \delta \boldsymbol{A} \right\| \left\| \boldsymbol{x} \right\| ,
$$

$$
\Longrightarrow \underset{\left\| \delta \boldsymbol{A} \right\| \le \alpha}{\mathrm{sup}}\frac{\frac{\left\| \delta \boldsymbol{x} \right\|}{\left\| \boldsymbol{x} \right\|}}{\frac{\left\| \delta \boldsymbol{A} \right\|}{\left\| \boldsymbol{A} \right\|}}\le \underset{\left\| \delta \boldsymbol{A} \right\| \le \alpha}{\mathrm{sup}}\frac{\left\| \boldsymbol{A} \right\| \left\| \boldsymbol{A}^{-1} \right\| \left\| \delta \boldsymbol{A} \right\| \left\| \boldsymbol{x} \right\|}{\left\| \delta \boldsymbol{A} \right\| \left\| \boldsymbol{x} \right\|}
\\
=\underset{\left\| \delta \boldsymbol{A} \right\| \le \alpha}{\mathrm{sup}}\left\| \boldsymbol{A} \right\| \left\| \boldsymbol{A}^{-1} \right\| =\left\| \boldsymbol{A} \right\| \left\| \boldsymbol{A}^{-1} \right\| 
$$

Based on the three conditions above, we can derive the conclusion:

$$
\left\| \boldsymbol{A} \right\| \left\| \boldsymbol{A}^{-1} \right\| =\kappa \left( \boldsymbol{A} \right) 
$$

### Extensions

#### More Examples

Here are a few more examples:

If:

$$
\boldsymbol{A}=\left[ \begin{matrix}
	\lambda _1&		&		0\\
	&		\ddots&		\\
	0&		&		\lambda _n\\
\end{matrix} \right] , \lambda _i\ne 0,\lambda _1\geqslant \lambda _2\geqslant \cdots \geqslant \lambda _n
$$

then $\kappa \left( \boldsymbol{A} \right) =\frac{\lambda _1}{\lambda _n}$

Let's look at another example. If:

$$
\boldsymbol{A}=\boldsymbol{U}\mathbf{\Sigma }\boldsymbol{V}^*, \mathbf{\Sigma }=\left[ \begin{matrix}
	\sigma _1&		&		0\\
	&		\ddots&		\\
	0&		&		\sigma _n\\
\end{matrix} \right] , \sigma _1\geqslant \sigma _2\geqslant \cdots \geqslant \sigma _n>0
$$

then:

$$
\left\| \boldsymbol{A} \right\| _2=\sigma _1,
\\
\left\| \boldsymbol{A}^{-1} \right\| _2=\frac{1}{\sigma _n}, \left( \boldsymbol{A}^{-1}=\boldsymbol{V}\mathbf{\Sigma }^{-1}\boldsymbol{U}^* \right) 
\\
\Longrightarrow \kappa \left( \boldsymbol{A} \right) =\left\| \boldsymbol{A} \right\| _2\cdot \left\| \boldsymbol{A}^{-1} \right\| _2=\frac{\sigma _1}{\sigma _n}
$$

#### Useful Properties

The $\kappa \left( \boldsymbol{A} \right)$ satisfies $\kappa \left( \boldsymbol{A} \right) \geqslant 1$ .

If $\kappa \left( \boldsymbol{A} \right)$ is large, the question is called **ill-conditioned**. If $\kappa \left( \boldsymbol{A} \right)$ is small, it is called **well-conditioned**.

