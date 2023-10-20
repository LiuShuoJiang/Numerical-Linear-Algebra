# Classical Eigenvalue Algorithms

## Preprocessing: Phase One

### General Steps for Phase One

Eigenvalue problem phase 1: $\boldsymbol{A}\rightarrow \boldsymbol{H}\,\,\left( \boldsymbol{A}\sim \boldsymbol{H} \right)$, $\boldsymbol{H}$ is a upper Hessenberg matrix. We get $\boldsymbol{A}=\boldsymbol{XHX}^{-1}$.

For example, let $\boldsymbol{A}\in \mathbb{R} ^{n\times n}$, the steps for phase 1 can be shown as:

$$
\boldsymbol{A}=\left[ \begin{matrix}
	\times&		\times&		\times&		\times\\
	\times&		\times&		\times&		\times\\
	\times&		\times&		\times&		\times\\
	\times&		\times&		\times&		\times\\
\end{matrix} \right] \xrightarrow{{\boldsymbol{Q}_1}^*\left( \mathrm{left} \right)}\left[ \begin{matrix}
	\times&		\times&		\times&		\times\\
	\times&		\times&		\times&		\times\\
	0&		\times&		\times&		\times\\
	0&		\times&		\times&		\times\\
\end{matrix} \right] \xrightarrow[\boldsymbol{Q}_1\left( \mathrm{right} \right)]{}\left[ \begin{matrix}
	\times&		\times&		\times&		\times\\
	\times&		\times&		\times&		\times\\
	0&		\times&		\times&		\times\\
	0&		\times&		\times&		\times\\
\end{matrix} \right] 
$$

$$
\xrightarrow{{\boldsymbol{Q}_2}^*\left( \mathrm{left} \right)}\left[ \begin{matrix}
	\times&		\times&		\times&		\times\\
	\times&		\times&		\times&		\times\\
	0&		\times&		\times&		\times\\
	0&		0&		\times&		\times\\
\end{matrix} \right] \xrightarrow[\boldsymbol{Q}_2\left( \mathrm{right} \right)]{}\left[ \begin{matrix}
	\times&		\times&		\times&		\times\\
	\times&		\times&		\times&		\times\\
	0&		\times&		\times&		\times\\
	0&		0&		\times&		\times\\
\end{matrix} \right] =\boldsymbol{H}
$$

### Householder Reduction

Algorithm: **Householder Reduction** to upper Hessenberg form:

For $k=1:m-2$:

- $\boldsymbol{x}=\boldsymbol{A}_{k+1:m,k}$
- $\boldsymbol{v}_k=\mathrm{sign}\left( x_1 \right) \left\| \boldsymbol{x} \right\| \boldsymbol{e}_1+\boldsymbol{x}$
- $\boldsymbol{v}_k=\frac{\boldsymbol{v}_k}{\left\| \boldsymbol{v}_k \right\|}$
- $\boldsymbol{A}_{k+1:m,k:m}=\boldsymbol{A}_{k+1:m,k:m}-2\boldsymbol{v}_k\left( {\boldsymbol{v}_k}^*\cdot \boldsymbol{A}_{k+1:m,k:m} \right)$
- $\boldsymbol{A}_{1:m,k+1:m}=\boldsymbol{A}_{1:m,k+1:m}-2\left( \boldsymbol{A}_{1:m,k+1:m}\cdot \boldsymbol{v}_k \right) {\boldsymbol{v}_k}^*$

End

The cost of the algorithm is $O\left( \frac{10}{3}m^3 \right)$.

Why does symmetric positive definite (***SPD***) matrix is easier to compute? Most stable situation in physics.

If $\boldsymbol{A}$ is symmetric, the cost is reduced to $O\left( \frac{4}{3}m^3 \right)$, and $\boldsymbol{H}$ is tridiagonal and also symmetric.

`Theorem` (***backward stability***): Let the Hessenberg reduction $\boldsymbol{A}=\boldsymbol{QHQ}^*$ be computed by the householder algorithm, then we have $\tilde{\boldsymbol{Q}}\tilde{\boldsymbol{H}}\tilde{\boldsymbol{Q}}^*=\boldsymbol{A}+\delta \boldsymbol{A}$ where for some $\delta \boldsymbol{A}$ satisfying $\frac{\left\| \delta \boldsymbol{A} \right\|}{\left\| \boldsymbol{A} \right\|}=O\left( \varepsilon _{\mathrm{machine}} \right)$.

We discuss some classical eigenvalue solvers below. Assume $\boldsymbol{A}$ is SPD, and eigenvalues $\lambda _1,\lambda _2,\cdots \lambda$ satisfy $\lambda _1\geqslant \lambda _2\geqslant \cdots \geqslant \lambda _m$. Their corresponding eigenvectors are $\boldsymbol{q}_1,\boldsymbol{q}_2,\cdots \boldsymbol{q}_m$.

## Power Iteration

**Power Iteration** algorithm:

Pick $\boldsymbol{v}^{\left( 0 \right)}$ as some vector with $\left\| \boldsymbol{v}^{\left( 0 \right)} \right\| =1$.

For $k=1,2,\cdots$:

- $\boldsymbol{w}=\boldsymbol{Av}^{\left( k-1 \right)}$
- $\boldsymbol{v}^{\left( k \right)}=\frac{\boldsymbol{w}}{\left\| \boldsymbol{w} \right\|}$
- $\lambda ^{\left( k \right)}=\left( \boldsymbol{v}^{\left( k \right)} \right) ^T\boldsymbol{Av}^{\left( k \right)}$

End

Claim of convergence: $\lambda ^{\left( k \right)}\rightarrow \lambda _1,\boldsymbol{v}^{\left( k \right)}\rightarrow \boldsymbol{q}_1$ as $k\rightarrow +\infty$.

Proof: We can write $\boldsymbol{v}^{\left( k \right)}=c_k\boldsymbol{A}^k\boldsymbol{v}^{\left( 0 \right)}$. we can also get $\boldsymbol{v}^{\left( 0 \right)}=a_1\boldsymbol{q}_1+a_2\boldsymbol{q}_2+\cdots +a_m\boldsymbol{q}_m$. Therefore, 

$$
\boldsymbol{v}^{\left( k \right)}=c_k\boldsymbol{A}^k\boldsymbol{v}^{\left( 0 \right)}=c_k\boldsymbol{A}^k\left( a_1\boldsymbol{q}_1+a_2\boldsymbol{q}_2+\cdots +a_m\boldsymbol{q}_m \right) 
$$

$$
=c_k\left( a_1\boldsymbol{A}^k\boldsymbol{q}_1+a_2\boldsymbol{A}^k\boldsymbol{q}_2+\cdots +a_m\boldsymbol{A}^k\boldsymbol{q}_m \right) 
$$

$$
=c_k\left( a_1{\lambda _1}^k\boldsymbol{q}_1+a_2{\lambda _2}^k\boldsymbol{q}_2+\cdots +a_m{\lambda _m}^k\boldsymbol{q}_m \right) 
$$

$$
=c_k{\lambda _1}^k\left( a_1\boldsymbol{q}_1+a_2\frac{{\lambda _2}^k}{{\lambda _1}^k}\boldsymbol{q}_2+\cdots +a_m\frac{{\lambda _m}^k}{{\lambda _1}^k}\boldsymbol{q}_m \right) 
$$

We know that $\left| \frac{\lambda _i}{\lambda _1} \right|<1, i\geqslant 2$. Therefore, $\left( \frac{{\lambda _i}^k}{{\lambda _1}^k} \right) \rightarrow 0$ as $k\rightarrow +\infty$. We can get the result.




