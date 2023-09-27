# Stability

Stability is for the algorithms.

Two concepts:

- [***Backward Stability***](#backward-stability)
- [***General Stability***](#general-stability)

Problem: $f: X\mapsto Y$

Algorithm: $\tilde{f}: X\mapsto Y$

The algorithm $\tilde{f}$ is used to solve the problem $f$

## Backward Stability

### Definition of Backward Stability

`Definition` (*Backward Stability*): An algorithm $\tilde{f}$ is called **backward stable** if for each $x\in X$ :

$$
\tilde{f}\left( x \right) =f\left( \tilde{x} \right) 
$$

for some $\tilde{x}$ with $\frac{\left\| x-\tilde{x} \right\|}{\left\| x \right\|}\leqslant O\left( \varepsilon _{\mathrm{machine}} \right)$

> *In words, a backward stable algorithm gives exactly the right solution to a nearly right problem.*

`Theorem`: If $\tilde{f}$ is a backward stable algorithm for a problem $f: X\mapsto Y$ with condition number $\kappa \left( f \right)$ , then the relative error satisfies:

$$
\frac{\left\| \tilde{f}\left( x \right) -f\left( x \right) \right\|}{\left\| f\left( x \right) \right\|}=O\left( \kappa \left( f \right) \cdot \varepsilon _{\mathrm{machine}} \right) 
$$

`Proof`: $\tilde{f}$ is backward stable: $\forall x\in X$ , $\exists \tilde{x}\in X$ such that $\tilde{f}\left( x \right) =f\left( \tilde{x} \right) , \frac{\left\| x-\tilde{x} \right\|}{\left\| x \right\|}\leqslant O\left( \varepsilon _{\mathrm{machine}} \right)$ . Then we can get:

$$
\frac{\left\| \tilde{f}\left( x \right) -f\left( x \right) \right\|}{\left\| f\left( x \right) \right\|}=\frac{\left\| f\left( \tilde{x} \right) -f\left( x \right) \right\|}{\left\| f\left( x \right) \right\|}\leqslant \kappa \left( f \right) \cdot \frac{\left\| \tilde{x}-x \right\|}{\left\| x \right\|}\leqslant O\left( \kappa \left( f \right) \cdot \varepsilon _{\mathrm{machine}} \right) 
$$

`Remark`: If $\boldsymbol{Ax}=\boldsymbol{b}$ is ill-conditioned, one must expect to lose $\log _{10}\kappa \left( \boldsymbol{A} \right)$ digits in computing the solution with a backward stable method except under very special cases.

### Example of Backward Stability

`Question`: How can we know an algorithm is backward stable or not? (tough and tedious)

`Example`: *Forward substitution* is backward stable.

Consider $\boldsymbol{Rx}=\boldsymbol{b}$ , where $\boldsymbol{R}\in \mathbb{R} ^{m\times m}$ is a lower triangular matrix. The forward substitution steps are:

$$
x_j=\small{\frac{b_j-\sum_{k=1}^{j-1}{x_kr_{kj}}}{r_{jj}}}
$$

where $r_{ij}$ is the element of $\boldsymbol{R}$ .

Denote $+,-,\times ,\div$ as exact operations and $\oplus ,\ominus ,\otimes ,\oslash$ as operations with error.

#### Case 1

Case $m=1$ : Exact version:

$$
r_{11}x_1=b_1\Leftrightarrow x_1=\frac{b_1}{r_{11}}=b_1\div r_{11}
$$

Computed version:

$$
\tilde{x}_1=b_1\oslash r_{11}=\frac{b_1}{r_{11}}\left( 1+\varepsilon _1 \right) ,
\\
\varepsilon _1\sim O\left( \varepsilon _{\mathrm{machine}} \right) 
$$

Introduce $\left( 1+\varepsilon _1 \right) \left( 1+\varepsilon _1\prime \right) =1$ . We can rewrite as:

$$
\tilde{x}_1=\frac{b_1}{r_{11}\left( 1+\varepsilon _1\prime \right)}\Longleftrightarrow r_{11}\left( 1+\varepsilon _1\prime \right) \tilde{x}_1=b_1
$$

The outcome can be viewed as an exact solution for a slightly perturbed problem.

#### Case 2

Case $m=2$ : 

Step 1: 

$$
\tilde{x}_1=b_1\oslash r_{11}=\frac{b_1}{r_{11}\left( 1+\varepsilon _1\prime \right)}
$$

Step 2:

$$
\tilde{x}_2=\left[ b_2\ominus \left( \tilde{x}_1\otimes r_{21} \right) \right] \oslash r_{22}=\left[ b_2\ominus \left( \tilde{x}_1\cdot r_{21} \right) \cdot \left( 1+\varepsilon _2 \right) \right] \oslash r_{22}
\\
=\left[ b_2-\left( \tilde{x}_1\cdot r_{21} \right) \left( 1+\varepsilon _2 \right) \right] \left( 1+\varepsilon _3 \right) \oslash r_{22}
\\
=\left[ b_2-\left( \tilde{x}_1\cdot r_{21} \right) \left( 1+\varepsilon _2 \right) \right] \left( 1+\varepsilon _3 \right) \div r_{22}\left( 1+\varepsilon _4 \right) 
\\
=\frac{b_2-\left( \tilde{x}_1\cdot r_{21} \right) \left( 1+\varepsilon _2 \right)}{r_{22}}\left( 1+\varepsilon _3 \right) \left( 1+\varepsilon _4 \right) 
$$

Introduce $\left( 1+\varepsilon _i \right) \left( 1+\varepsilon _i\prime \right) =1$ , then:

$$
\tilde{x}_2=\frac{b_2-\left( \tilde{x}_1\cdot r_{21} \right) \left( 1+\varepsilon _2 \right)}{r_{22}\left( 1+\varepsilon _3\prime \right) \left( 1+\varepsilon _4\prime \right)}
\\
=\frac{b_2-\tilde{x}_1\cdot \left( r_{21}\left( 1+\varepsilon _2 \right) \right)}{r_{22}\left( 1+2\varepsilon _5\prime \right)}
$$

Then:

$$
\left( r_{21}\left( 1+\varepsilon _2 \right) \right) \cdot \tilde{x}_1+r_{22}\left( 1+2\varepsilon _5\prime \right) \cdot \tilde{x}_2=b_2
$$

We can claim that solving $\boldsymbol{Rx}=\boldsymbol{b}$ by a computer with a numerical solution $\hat{\boldsymbol{x}}$ can be viewed as $\left( \boldsymbol{R}+\delta \boldsymbol{R} \right) \hat{\boldsymbol{x}}=\boldsymbol{b}$ exactly. Also:

$$
\delta \boldsymbol{R}=\left[ \begin{matrix}
	\varepsilon _1\prime r_{11}&		0\\
	\varepsilon _2r_{21}&		2\varepsilon _5\prime r_{22}\\
\end{matrix} \right]
$$

#### General Case

For general $m$ , $\hat{\boldsymbol{x}}$ satisfies:

$$
\frac{\left\| \delta \boldsymbol{R} \right\|}{\left\| \boldsymbol{R} \right\|}\leqslant \left\| \left[ \begin{matrix}
	1&		&		&		O\\
	\vdots&		2&		&		\\
	\vdots&		\vdots&		\ddots&		\\
	m-1&		m-1&		\cdots&		m\\
\end{matrix} \right] \right\| \cdot \varepsilon _{\mathrm{machine}}
$$

Finally, we can claim that forward substitution algorithm computes the solution $\hat{\boldsymbol{x}}$ of $\boldsymbol{Rx}=\boldsymbol{b}$ satisfying 

$$
\left( \boldsymbol{R}+\delta \boldsymbol{R} \right) \hat{\boldsymbol{x}}=\boldsymbol{b}
$$

for some $\delta \boldsymbol{R}$ satisfying

$$
\frac{\left\| \delta \boldsymbol{R} \right\|}{\left\| \boldsymbol{R} \right\|}=O\left( \varepsilon _{\mathrm{machine}} \right) 
$$

Also, forward and backward substitutions are all backward stable.

## General Stability

### Definition of General Stability

Problem: $f: X\mapsto Y$

Algorithm: $\tilde{f}: X\mapsto Y$

`Definition`: $\tilde{f}$ is called **stable** if for $\forall x\in X$ ,

$$
\frac{\left\| \tilde{f}\left( x \right) -f\left( \tilde{x} \right) \right\|}{\left\| f\left( \tilde{x} \right) \right\|}=O\left( \kappa \left( f \right) \cdot \varepsilon _{\mathrm{machine}} \right) 
$$

for some $\tilde{x}$ satisfying $\frac{\left\| \tilde{x}-x \right\|}{\left\| x \right\|}=O\left( \varepsilon _{\mathrm{machine}} \right)$

`Conclusion`: If $\tilde{f}$ is backward stable, then it is also stable.

> *In words, a stable algorithm gives nearly the right solution to nearly the right algorithm.*

### Conclusions on Gaussian Elimination

`Theorem`: Let $\boldsymbol{A}=\boldsymbol{LU}$ of a nonsingular matrix $\boldsymbol{A}\in \mathbb{R} ^{m\times m}$ be computed by *Gaussian Elimination without pivoting*, then:

$$
\tilde{\boldsymbol{L}}\tilde{\boldsymbol{U}}=\boldsymbol{A}+\delta \boldsymbol{A}
$$

for some $\delta \boldsymbol{A}\in \mathbb{R} ^{m\times m}$ satisfying $\frac{\left\| \delta \boldsymbol{A} \right\|}{\left\| \boldsymbol{L} \right\| \left\| \boldsymbol{U} \right\|}=O\left( \varepsilon _{\mathrm{machine}} \right)$

If $\left\| \boldsymbol{L} \right\| \left\| \boldsymbol{U} \right\| =O\left( \left\| \boldsymbol{A} \right\| \right)$ , then Gaussian Elimination without pivoting is stable. Otherwise, the method is not stable.

`Theorem`: *Gaussian Elimination with partial pivoting* is backward stable:

$$
\tilde{\boldsymbol{L}}\tilde{\boldsymbol{U}}=\tilde{\boldsymbol{P}}\left( \boldsymbol{A}+\delta \boldsymbol{A} \right) 
$$

for some $\delta \boldsymbol{A}\in \mathbb{R} ^{m\times m}$ satisfying $\frac{\left\| \delta \boldsymbol{A} \right\|}{\left\| \boldsymbol{A} \right\|}=O\left( \rho \cdot \varepsilon _{\mathrm{machine}} \right)$ where

$$
\rho =\frac{\max_{i,j} \left| \tilde{u}_{ij} \right|}{\max_{i,j} \left| a_{ij} \right|}
$$

is called **growth factor**.

`Theorem`: If $\boldsymbol{A}$ is symmetric positive definite (SPD), then *Cholesky Factorization* is backward stable:

$$
\tilde{\boldsymbol{R}}^T\tilde{\boldsymbol{R}}=\boldsymbol{A}+\delta \boldsymbol{A}
$$

for some $\delta \boldsymbol{A}$ satisfying $\frac{\left\| \delta \boldsymbol{A} \right\|}{\left\| \boldsymbol{A} \right\|}=O\left( \varepsilon _{\mathrm{machine}} \right)$
