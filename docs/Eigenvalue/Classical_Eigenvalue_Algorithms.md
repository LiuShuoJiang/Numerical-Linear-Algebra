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

Algorithm (**Householder Reduction** to upper Hessenberg form):

For $k=1:m-2$:

- $\boldsymbol{x}=\boldsymbol{A}_{k+1:m,k}$;
- $\boldsymbol{v}_k=\mathrm{sign}\left( x_1 \right) \left\| \boldsymbol{x} \right\| \boldsymbol{e}_1+\boldsymbol{x}$;
- $\boldsymbol{v}_k=\frac{\boldsymbol{v}_k}{\left\| \boldsymbol{v}_k \right\|}$;
- $\boldsymbol{A}_{k+1:m,k:m}=\boldsymbol{A}_{k+1:m,k:m}-2\boldsymbol{v}_k\left( {\boldsymbol{v}_k}^*\cdot \boldsymbol{A}_{k+1:m,k:m} \right)$;
- $\boldsymbol{A}_{1:m,k+1:m}=\boldsymbol{A}_{1:m,k+1:m}-2\left( \boldsymbol{A}_{1:m,k+1:m}\cdot \boldsymbol{v}_k \right) {\boldsymbol{v}_k}^*$;

End

The cost of the algorithm is $O\left( \frac{10}{3}m^3 \right)$.

Why does symmetric positive definite (***SPD***) matrix is easier to compute? Most stable situation in physics.

If $\boldsymbol{A}$ is symmetric, the cost is reduced to $O\left( \frac{4}{3}m^3 \right)$, and $\boldsymbol{H}$ is tridiagonal and also symmetric.

`Theorem` (***backward stability***): Let the Hessenberg reduction $\boldsymbol{A}=\boldsymbol{QHQ}^*$ be computed by the householder algorithm, then we have $\tilde{\boldsymbol{Q}}\tilde{\boldsymbol{H}}\tilde{\boldsymbol{Q}}^*=\boldsymbol{A}+\delta \boldsymbol{A}$ where for some $\delta \boldsymbol{A}$ satisfying $\frac{\left\| \delta \boldsymbol{A} \right\|}{\left\| \boldsymbol{A} \right\|}=O\left( \varepsilon _{\mathrm{machine}} \right)$.

We discuss some classical eigenvalue solvers below. Assume $\boldsymbol{A}$ is SPD, and eigenvalues $\lambda _1,\lambda _2,\cdots \lambda$ satisfy $\lambda _1\geqslant \lambda _2\geqslant \cdots \geqslant \lambda _m$. Their corresponding eigenvectors are $\boldsymbol{q}_1,\boldsymbol{q}_2,\cdots \boldsymbol{q}_m$.

## Power Iteration

Algorithm (**[Power Iteration](https://en.wikiversity.org/wiki/Numerical_Analysis/Power_iteration_examples)**):

Pick $\boldsymbol{v}^{\left( 0 \right)}$ as some vector with $\left\| \boldsymbol{v}^{\left( 0 \right)} \right\| =1$;

For $k=1,2,\cdots$:

- $\boldsymbol{w}=\boldsymbol{Av}^{\left( k-1 \right)}$;
- $\boldsymbol{v}^{\left( k \right)}=\frac{\boldsymbol{w}}{\left\| \boldsymbol{w} \right\|}$;
- $\lambda ^{\left( k \right)}=\left( \boldsymbol{v}^{\left( k \right)} \right) ^T\boldsymbol{Av}^{\left( k \right)}$;

End

`Claim of Convergence`: $\lambda ^{\left( k \right)}\rightarrow \lambda _1,\boldsymbol{v}^{\left( k \right)}\rightarrow \boldsymbol{q}_1$ as $k\rightarrow +\infty$.

`Proof`: We can write $\boldsymbol{v}^{\left( k \right)}=c_k\boldsymbol{A}^k\boldsymbol{v}^{\left( 0 \right)}$. We can also get $\boldsymbol{v}^{\left( 0 \right)}=a_1\boldsymbol{q}_1+a_2\boldsymbol{q}_2+\cdots +a_m\boldsymbol{q}_m$. Therefore, 

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

If ${\lambda}_1 > {\lambda}_2$, we know that $\left| \frac{\lambda _i}{\lambda _1} \right|<1, i\geqslant 2$. Therefore, $\left( \frac{{\lambda _i}^k}{{\lambda _1}^k} \right) \rightarrow 0$ as $k\rightarrow +\infty$. We can get the result.

In conclusion, $\boldsymbol{v}^{\left( k \right)}\rightarrow \pm \boldsymbol{q}_1; \lambda ^{\left( k \right)}=\left( \boldsymbol{v}^{\left( k \right)} \right) ^T\boldsymbol{Av}^{\left( k \right)}\rightarrow {\boldsymbol{q}_1}^T\boldsymbol{Aq}_1=\lambda _1$.

`Theorem`: Suppose $\left| \lambda _1 \right|>\left| \lambda _2 \right|\geqslant \cdots \geqslant \left| \lambda _m \right|>0$, and ${\boldsymbol{q}_1}^T\boldsymbol{v}^{\left( 0 \right)}\ne 0 \left( a_1\ne 0 \right)$, then the Power Iteration has:

$$
\left\| \boldsymbol{v}^{\left( k \right)}-\left( \pm \boldsymbol{q}_1 \right) \right\| =O\left( \left| \frac{\lambda _2}{\lambda _1} \right|^k \right) ;
$$

$$
\left\| \lambda ^{\left( k \right)}-\lambda _1 \right\| =O\left( \left| \frac{\lambda _2}{\lambda _1} \right|^{2k} \right) 
$$

as $k\rightarrow +\infty$. $\pm 1$ sign means one of them is correct.

**Comments**:

1. ${\boldsymbol{q}_1}^T\boldsymbol{v}^{\left( 0 \right)}\ne 0$ is only required for the theory purpose. In practice, this is rarely needed because the computer may perturb.
2. The eigenvalue convergence is ***quadratic***, while the eigenvector convergence is ***linear***.
3. Power Iteration only gives $\lambda _1$ and $\boldsymbol{q}_1$.
4. The rate of convergence is determined by $\left| \frac{\lambda _2}{\lambda _1} \right|$ ("Spectral Gap"). If $\left| \frac{\lambda _2}{\lambda _1} \right|\sim 1$, the convergence is slow.

To capture the smallest (in absolute value sense) eigenvalue, we can use Power Iteration to $\boldsymbol{A}^{-1}$: $\boldsymbol{A}^{-1}\boldsymbol{v}^{\left( 0 \right)}=\boldsymbol{w}, \boldsymbol{v}^{\left( 1 \right)}=\frac{\boldsymbol{w}}{\left\| \boldsymbol{w} \right\|}$. This gives us $\frac{1}{\lambda _m}, \boldsymbol{q}_m$. In practice, we should solve linear systems instead of getting the inverse directly: $\boldsymbol{A}^{-1}\boldsymbol{v}^{\left( 0 \right)}=\boldsymbol{w}\Longleftrightarrow \boldsymbol{Aw}=\boldsymbol{v}^{\left( 0 \right)}$. Do $\boldsymbol{A}=\boldsymbol{LU}$ first, then solve $\boldsymbol{LUw}=\boldsymbol{v}^{\left( 0 \right)}$ or generally $\boldsymbol{LUw}=\boldsymbol{v}^{\left( k \right)}$.

## Inverse Iteration with Shift

`Question`: How can we use Power Iteration to find the eigenvalue/eigenvector close to a number $\mu$ ?

**Inverse Iteration with Shift** can solve this problem.

`Claim`: If $\mu \in \mathbb{R}$ is not an eigenvalue of $\boldsymbol{A}$, then $\left( \boldsymbol{A}-\mu \mathbf{I} \right) ^{-1}$ and $\boldsymbol{A}$ has the same eigenvectors, and their corresponding eigenvalues are $\left( \lambda _i-\mu \right) ^{-1}$ and $\lambda _i$ respectively.

`Proof`: $\left( \boldsymbol{A}-\mu \mathbf{I} \right) ^{-1}$ exists. If $\boldsymbol{q}_i$ is an eigenvector of $\boldsymbol{A}$, then $\boldsymbol{Aq}_i=\lambda _i\boldsymbol{q}_i$. We can get $\boldsymbol{Aq}_i-\mu \boldsymbol{q}_i=\left( \lambda _i-\mu \right) \boldsymbol{q}_i\Longleftrightarrow \left( \boldsymbol{A}-\mu \mathbf{I} \right) \boldsymbol{q}_i=\left( \lambda _i-\mu \right) \boldsymbol{q}_i$, and $\frac{1}{\lambda _i-\mu}\boldsymbol{q}_i=\left( \boldsymbol{A}-\mu \mathbf{I} \right) ^{-1}\boldsymbol{q}_i$. Therefore, $\boldsymbol{q}_i$ is an eigenvector of $\left( \boldsymbol{A}-\mu \mathbf{I} \right) ^{-1}$.

If $\mu$ is the closest to $\lambda _J$, then $\left( \lambda _J-\mu \right) ^{-1}$ is the largest eigenvalue for $\left( \boldsymbol{A}-\mu \mathbf{I} \right) ^{-1}$.

Apply Power Iteration to $\left( \boldsymbol{A}-\mu \mathbf{I} \right) ^{-1}$, we can obtain $\boldsymbol{q}_J$ and $\left( \lambda _J-\mu \right) ^{-1}$.

Algorithm (**Inverse Iteration with Shift**):

$\boldsymbol{v}^{\left( 0 \right)}$ is selected with $\left\| \boldsymbol{v}^{\left( 0 \right)} \right\| =1$; Also solve $\boldsymbol{LU}=\boldsymbol{A}-\mu \mathbf{I}$;

For $k=1,2,\cdots$:

- Solve $\boldsymbol{LUw}=\boldsymbol{v}^{\left( k-1 \right)}$ for $\boldsymbol{w}$;
- $\boldsymbol{v}^{\left( k \right)}=\frac{\boldsymbol{w}}{\left\| \boldsymbol{w} \right\|}$;
- $\lambda ^{\left( k \right)}=\left( \boldsymbol{v}^{\left( k \right)} \right) ^T\boldsymbol{Av}^{\left( k \right)}$;

End

`Theorem`: Assume $\lambda _J$ is closest to $\mu$, and $\left| \mu -\lambda _J \right|<\left| \mu -\lambda _K \right|\leqslant \left| \mu -\lambda _i \right|, \left( i\ne J \right)$. We can get:

$$
\left\| \boldsymbol{v}^{\left( k \right)}-\left| \pm \boldsymbol{q}_j \right| \right\| =O\left( \left| \frac{\mu -\lambda _J}{\mu -\lambda _K} \right|^k \right) ;
$$

$$
\left\| \lambda _J-\lambda ^{\left( k \right)} \right\| =O\left( \left| \frac{\mu -\lambda _J}{\mu -\lambda _K} \right|^{2k} \right) 
$$

as $k\rightarrow +\infty$.

## Rayleigh Quotient Iteration

### Rayleigh Quotient

Given matrix $\boldsymbol{A}$, the **Rayleigh Quotient** of a vector $\boldsymbol{x}\ne 0$ is the ratio:

$$
r\left( \boldsymbol{x} \right) =\frac{\boldsymbol{x}^T\boldsymbol{Ax}}{\boldsymbol{x}^T\boldsymbol{x}}
$$

If $\boldsymbol{x}$ is an eigenvector, then $r\left( \boldsymbol{x} \right)$ is its corresponding eigenvalue: $r\left( \boldsymbol{x} \right) =\frac{\boldsymbol{x}^T\boldsymbol{Ax}}{\boldsymbol{x}^T\boldsymbol{x}}=\frac{\boldsymbol{x}^T\lambda \boldsymbol{x}}{\boldsymbol{x}^T\boldsymbol{x}}=\lambda$.

`Claim`: If $\boldsymbol{x}$ is close to an eigenvector, then $r\left( \boldsymbol{x} \right)$ is also close to the corresponding eigenvalue (due to continuity).

Consider $\frac{\partial r\left( \boldsymbol{x} \right)}{\partial x_j}=\frac{2}{\boldsymbol{x}^T\boldsymbol{x}}\cdot \left( \boldsymbol{Ax}-r\left( \boldsymbol{x} \right) \cdot \boldsymbol{x} \right) _j, j=1,2,\cdots ,m$ (How to get this result?). We can get:

$$
\nabla r\left( \boldsymbol{x} \right) =\left[ \begin{array}{c}
	\frac{\partial r\left( \boldsymbol{x} \right)}{\partial x_1}\\
	\vdots\\
	\frac{\partial r\left( \boldsymbol{x} \right)}{\partial x_m}\\
\end{array} \right] =\frac{2}{\boldsymbol{x}^T\boldsymbol{x}}\cdot \left( \boldsymbol{Ax}-r\left( \boldsymbol{x} \right) \cdot \boldsymbol{x} \right) 
$$

If $\boldsymbol{x}$ is an eigenvector of $\boldsymbol{A}$, then $\nabla r\left( \boldsymbol{x} \right) =\mathbf{0}$. If $\boldsymbol{x}\ne 0$ satisfies $\nabla r\left( \boldsymbol{x} \right) =\mathbf{0}$, then $\boldsymbol{x}$ is an eigenvector of $\boldsymbol{A}$.

Take the *Taylor Expansion* of $r\left( \boldsymbol{x} \right)$ at an eigenvector $\boldsymbol{q}_J$. We know:

$$
r\left( \boldsymbol{x} \right) =r\left( \boldsymbol{q}_J \right) +\left( \boldsymbol{x}-\boldsymbol{q}_J \right) \nabla r\left( \boldsymbol{q}_J \right) +\frac{1}{2}\left\| \boldsymbol{x}-\boldsymbol{q}_J \right\| ^2\nabla ^2r\left( \boldsymbol{q}_J \right) +\cdots 
$$

If $\boldsymbol{x}$ is in a small neighborhood of $\boldsymbol{q}_J$, we can get:

$$
r\left( \boldsymbol{x} \right) =r\left( \boldsymbol{q}_J \right) +\left( \boldsymbol{x}-\boldsymbol{q}_J \right) \underset{0}{\underbrace{\nabla r\left( \boldsymbol{q}_J \right) }}+O\left( \left\| \boldsymbol{x}-\boldsymbol{q}_J \right\| ^2 \right) 
$$

$$
=r\left( \boldsymbol{q}_J \right) +O\left( \left\| \boldsymbol{x}-\boldsymbol{q}_J \right\| ^2 \right) =\lambda _J+O\left( \left\| \boldsymbol{x}-\boldsymbol{q}_J \right\| ^2 \right) 
$$

### Rayleigh Quotient Iteration Algorithm

Algorithm (**Rayleigh Quotient Iteration**):

Pick vector $\boldsymbol{v}^{\left( 0 \right)}$ as $\left\| \boldsymbol{v}^{\left( 0 \right)} \right\| =1$; Also $\lambda ^{\left( 0 \right)}=\boldsymbol{v}^{\left( 0 \right)}\boldsymbol{Av}^{\left( 0 \right)}$;

For $k=1,2,\cdots$:

- Solve $\left( \boldsymbol{A}-\lambda ^{\left( k-1 \right)}\mathbf{I} \right) \boldsymbol{w}=\boldsymbol{v}^{\left( k-1 \right)}$ for $\boldsymbol{w}$;
- $\boldsymbol{v}^{\left( k \right)}=\frac{\boldsymbol{w}}{\left\| \boldsymbol{w} \right\|}$;
- $\lambda ^{\left( k \right)}=\left( \boldsymbol{v}^{\left( k \right)} \right) ^T\boldsymbol{Av}^{\left( k \right)}$;

End

`Theorem`: Rayleigh Quotient Iteration converges to an eigenvalue/eigenvector pair for all except a set of *measure zero* starting vectors, i.e., $\left( \boldsymbol{v}^{\left( 0 \right)} \right) ^T\boldsymbol{q}_J\ne 0$. When it converges, the ultimate rate of convergence is ***cubic***:

$$
\left\| \boldsymbol{v}^{\left( k \right)}-\left( \pm \boldsymbol{q}_J \right) \right\| =O\left( \left\| \boldsymbol{v}^{\left( k-1 \right)}-\left( \pm \boldsymbol{q}_J \right) \right\| ^3 \right) ;
$$

$$
\left\| \lambda ^{\left( k \right)}-\lambda _J \right\| =O\left( \left| \lambda ^{\left( k-1 \right)}-\lambda _J \right|^3 \right) 
$$

**Remark**:

1. The cubic convergence is super fast. E.g., if $\left\| \boldsymbol{v}^{\left( 0 \right)}-\left( \pm \boldsymbol{q}_J \right) \right\| =0.1$, then $\left\| \boldsymbol{v}^{\left( 1 \right)}-\left( \pm \boldsymbol{q}_J \right) \right\| =O\left( 0.1^3 \right) =10^{-3}$, $\left\| \boldsymbol{v}^{\left( 2 \right)}-\left( \pm \boldsymbol{q}_J \right) \right\| =10^{-9}$, $\left\| \boldsymbol{v}^{\left( 3 \right)}-\left( \pm \boldsymbol{q}_J \right) \right\| =10^{-27}$.
2. It reduces the number of iterations, but since we have to solve the linear system for each iteration, each iteration has much higher cost. Rayleigh Quotient Iteration is less frequently used by Inverse Iteration in practice.

