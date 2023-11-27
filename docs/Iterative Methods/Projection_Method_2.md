# Projection Method (Continued)

## Krylov Subspace

A different way to introduce CG method:

**Krylov Subspace**: Given $\boldsymbol{r}_0$, let

$$
\mathbf{K}_n=\mathrm{span}\left\{ \boldsymbol{r}_0,\boldsymbol{Ar}_0,\boldsymbol{A}^2\boldsymbol{r}_0,\cdots ,\boldsymbol{A}^{n-1}\boldsymbol{r}_0 \right\} \triangleq \mathbf{K}_n\left( \boldsymbol{A},\boldsymbol{r}_0 \right) 
$$

`Claim`: If we take $\boldsymbol{L}_n=\mathbf{K}_n$, when $\boldsymbol{A}$ is SPD, by applying projection method, we obtain the CG method.

We need to have basis for $\boldsymbol{L}_n=\mathbf{K}_n$.

`Goal`: Construct an orthonormal basis for $\mathbf{K}_n$, denoted by $\boldsymbol{V}_n$.

We can apply Modified Gram-Schmidt algorithm to find the orthonormal basis for the Krylov Subspace $\mathbf{K}_n$. This resulting algorithm is called **Lancozs algorithm**.

Algorithm ( **Lancozs** ):

Choose an initial vector $\boldsymbol{v}_1$ with $\left\| \boldsymbol{v}_1 \right\| =1$;

Set $\beta _1=0$;

For $j=1,2,\cdots ,n$, do:

- $\boldsymbol{w}_j=\boldsymbol{Av}_j-\beta _j\boldsymbol{v}_{j-1}$;
- $\alpha _j=\left< \boldsymbol{w}_j,\boldsymbol{v}_j \right>$;
- $\boldsymbol{w}_j=\boldsymbol{w}_j-\alpha _j\boldsymbol{v}_j$;
- $\beta _{j+1}=\left\| \boldsymbol{w}_j \right\|$; If $\beta _{j+1}=0$, stop;
- $\boldsymbol{v}_{j+1}=\frac{\boldsymbol{w}_j}{\beta _{j+1}}$;

End

$\left\{ \boldsymbol{v}_1,\cdots ,\boldsymbol{v}_n \right\}$ form an orthonormal basis for $\mathbf{K}_n$.

$$
\beta _{j+1}\boldsymbol{v}_{j+1}=\boldsymbol{w}_j=\boldsymbol{Av}_j-\alpha _j\boldsymbol{v}_j-\beta _j\boldsymbol{v}_{j-1};
$$

$$
\boldsymbol{Av}_j=\beta _{j+1}\boldsymbol{v}_{j+1}+\alpha _j\boldsymbol{v}_j+\beta _j\boldsymbol{v}_{j-1}
$$

Note that:

$$
\boldsymbol{V}_n=\left[ \begin{matrix}
	\boldsymbol{v}_1&		\cdots&		\boldsymbol{v}_n\\
\end{matrix} \right] , \boldsymbol{AV}_n=\boldsymbol{V}_n\boldsymbol{T}_n;
$$

$$
\boldsymbol{T}_n=\left[ \begin{matrix}
	\alpha _1&		\beta _2&		&		&		&		\boldsymbol{O}\\
	\beta _2&		\alpha _2&		\beta _3&		&		&		\\
	&		\beta _3&		\alpha _3&		\beta _4&		&		\\
	&		&		\ddots&		\ddots&		\ddots&		\\
	&		&		&		\ddots&		\ddots&		\beta _n\\
	\boldsymbol{O}&		&		&		&		\beta _n&		\alpha _n\\
\end{matrix} \right] 
$$

We know that:

$$
{\boldsymbol{V}_n}^T\boldsymbol{AV}_n=\boldsymbol{T}_n
$$

Apply projection method:

$$
\begin{cases}
	\boldsymbol{x}^{\left( n \right)}=\boldsymbol{x}^{\left( 0 \right)}+\boldsymbol{V}_n\boldsymbol{y}^{\left( n \right)}\\
	\boldsymbol{y}^{\left( n \right)}=\left( {\boldsymbol{V}_n}^T\boldsymbol{AV}_n \right) ^{-1}{\boldsymbol{V}_n}^T\boldsymbol{r}^{\left( 0 \right)}\\
\end{cases};
$$

$$
\boldsymbol{y}^{\left( n \right)}=\left( \boldsymbol{T}_n \right) ^{-1}\beta _1\boldsymbol{e}_1, \beta _1=\left\| \boldsymbol{r}^{\left( 0 \right)} \right\| 
$$

Apply LU factorization to find $\left( \boldsymbol{T}_n \right) ^{-1}$:

$$
\boldsymbol{T}_n=\boldsymbol{L}_n\boldsymbol{U}_n, 
$$

$$
\boldsymbol{L}_n=\left[ \begin{matrix}
	1&		&		&		\boldsymbol{O}\\
	\lambda _2&		\ddots&		&		\\
	&		\ddots&		\ddots&		\\
	\boldsymbol{O}&		&		\lambda _n&		1\\
\end{matrix} \right] , \boldsymbol{U}_n=\left[ \begin{matrix}
	\eta _1&		\beta _2&		&		\boldsymbol{O}\\
	&		\ddots&		\ddots&		\\
	&		&		\ddots&		\beta _n\\
	\boldsymbol{O}&		&		&		\eta _n\\
\end{matrix} \right] ;
$$

$$
\boldsymbol{x}^{\left( n \right)}=\boldsymbol{x}^{\left( 0 \right)}+\boldsymbol{V}_n\boldsymbol{y}^{\left( n \right)}=\boldsymbol{x}^{\left( 0 \right)}+\boldsymbol{V}_n\left( \boldsymbol{T}_n \right) ^{-1}\beta _1\boldsymbol{e}_1
$$

$$
=\boldsymbol{x}^{\left( 0 \right)}+\underset{\boldsymbol{P}_n}{\underbrace{\boldsymbol{V}_n\left( \boldsymbol{U}_n \right) ^{-1}}}\cdot \underset{\boldsymbol{z}_n}{\underbrace{\left( \boldsymbol{L}_n \right) ^{-1}\beta _1\boldsymbol{e}_1}}
$$

$$
=\boldsymbol{x}^{\left( 0 \right)}+\boldsymbol{P}_n\boldsymbol{z}_n=\boldsymbol{x}^{\left( 0 \right)}+\boldsymbol{P}_{n-1}\boldsymbol{z}_{n-1}+\xi _n\boldsymbol{p}_n
$$

$$
=\boldsymbol{x}^{\left( n-1 \right)}+\xi _n\boldsymbol{p}_n
$$

We need to determine $\xi _n\boldsymbol{p}_n$ such that $\boldsymbol{x}^{\left( 0 \right)}$ is the solution of the projection method.

Summarize what we have got so far:

$$
\boldsymbol{x}^{\left( n \right)}=\boldsymbol{x}^{\left( 0 \right)}+\boldsymbol{P}_n\boldsymbol{z}_n,
$$

$$
\boldsymbol{P}_n=\boldsymbol{V}_n{\boldsymbol{U}_n}^{-1}, \boldsymbol{z}_n={\boldsymbol{L}_n}^{-1}\beta _1\boldsymbol{e}_1;
$$

$$
\boldsymbol{P}_n=\left[ \begin{matrix}
	\boldsymbol{p}_1&		\boldsymbol{p}_2&		\cdots&		\boldsymbol{p}_n\\
\end{matrix} \right] ,
$$

$$
\boldsymbol{x}_n=\boldsymbol{x}^{\left( 0 \right)}+\boldsymbol{P}_{n-1}\boldsymbol{z}_{n-1}+\xi _n\boldsymbol{p}_n=\boldsymbol{x}^{\left( n-1 \right)}+\xi _n\boldsymbol{p}_n
$$

Note that:

$$
\boldsymbol{P}_n=\boldsymbol{V}_n{\boldsymbol{U}_n}^{-1}\Leftrightarrow \boldsymbol{P}_n\boldsymbol{U}_n=\boldsymbol{V}_n,
$$

$$
\left[ \begin{matrix}
	\boldsymbol{p}_1&		\boldsymbol{p}_2&		\cdots&		\boldsymbol{p}_n\\
\end{matrix} \right] \cdot \left[ \begin{matrix}
	\eta _1&		\beta _2&		0&		\cdots&		\boldsymbol{O}\\
	0&		\eta _2&		\beta _3&		\ddots&		\vdots\\
	&		\ddots&		\ddots&		\ddots&		0\\
	&		&		\ddots&		\ddots&		\beta _n\\
	\boldsymbol{O}&		&		&		0&		\eta _n\\
\end{matrix} \right] =\left[ \begin{matrix}
	\boldsymbol{v}_1&		\boldsymbol{v}_2&		\cdots&		\boldsymbol{v}_n\\
\end{matrix} \right] ,
$$

$$
\boldsymbol{v}_n=\beta _n\boldsymbol{p}_{n-1}+\eta _n\boldsymbol{p}_n
$$

$$
\Rightarrow \boldsymbol{p}_n=\frac{1}{\eta _n}\left( \boldsymbol{v}_n-\beta _n\boldsymbol{p}_{n-1} \right) 
$$

`Claim`: $\boldsymbol{v}_n\parallel \boldsymbol{r}_n$. (Why?  Prove by induction)

Then:

$$
\boldsymbol{p}_n=\frac{1}{\eta _n}\left( \gamma _n\boldsymbol{r}_n-\beta _n\boldsymbol{p}_{n-1} \right) 
$$

Now, this projection method becomes CG method.

## Properties of Conjugate Gradient Method

`Theorem`: In CG method,

- $\left< \boldsymbol{r}^{\left( j \right)},\boldsymbol{r}^{\left( i \right)} \right> =0, i\ne j$;
- $\left< \boldsymbol{Ap}^{\left( j \right)},\boldsymbol{p}^{\left( i \right)} \right> =0,i\ne j$.

`Proof`: Note that:

$$
\boldsymbol{P}_n=\boldsymbol{V}_n{\boldsymbol{U}_n}^{-1},
$$

$$
{\boldsymbol{P}_n}^T\boldsymbol{AP}_n=\left( {\boldsymbol{U}_n}^{-1} \right) ^T{\boldsymbol{V}_n}^T\boldsymbol{AV}_n{\boldsymbol{U}_n}^{-1}
$$

$$
=\left( {\boldsymbol{U}_n}^{-1} \right) ^T\boldsymbol{T}_n{\boldsymbol{U}_n}^{-1}=\left( {\boldsymbol{U}_n}^{-1} \right) ^T\boldsymbol{L}_n\boldsymbol{U}_n{\boldsymbol{U}_n}^{-1}=\left( {\boldsymbol{U}_n}^{-1} \right) ^T\boldsymbol{L}_n
$$

Since $\left( {\boldsymbol{U}_n}^{-1} \right) ^T\boldsymbol{L}_n$ is lower triangular, and ${\boldsymbol{P}_n}^T\boldsymbol{AP}_n$ is symmetric, then we know that ${\boldsymbol{P}_n}^T\boldsymbol{AP}_n$ is diagonal. Therefore,

$$
{\boldsymbol{p}_i}^T\boldsymbol{Ap}_j=0, i\ne j\Longleftrightarrow \left< \boldsymbol{Ap}^{\left( j \right)},\boldsymbol{p}^{\left( i \right)} \right> =0,i\ne j
$$

We call $\boldsymbol{p}^{\left( j \right)},\boldsymbol{p}^{\left( i \right)}$ are ***conjugate*** to each other.

`Theorem`: Assume $\boldsymbol{A}$ is SPD, $\boldsymbol{L}=\boldsymbol{K}$, then a vector $\boldsymbol{x}$ is the result of a projection method onto $\boldsymbol{K}$, orthogonal to $\boldsymbol{L}$ with initial guess $\boldsymbol{x}^{(0)}$ if and only if $\boldsymbol{x}$ minimizes the $\boldsymbol{A}$-norm of the error over $\boldsymbol{x}^{\left( 0 \right)}+\boldsymbol{K}$, i.e.,

$$
E\left( \boldsymbol{y} \right) =\left\| \boldsymbol{x}^*-\boldsymbol{y} \right\| _{\boldsymbol{A}}^{2}=\left( \boldsymbol{x}^*-\boldsymbol{y} \right) ^T\boldsymbol{A}\left( \boldsymbol{x}^*-\boldsymbol{y} \right) 
$$

where $\boldsymbol{x}^*$ is the true solution of $\boldsymbol{Ax}=\boldsymbol{b}$, and:

$$
\boldsymbol{x}=\mathrm{arg}\min E\left( \boldsymbol{y} \right) , \boldsymbol{y}\in \boldsymbol{x}^{\left( 0 \right)}+\boldsymbol{K}
$$

`Proof` (idea): To show that

$$
\frac{\partial E}{\partial \boldsymbol{y}}\mid_{\boldsymbol{y}=\boldsymbol{x}}^{}=\boldsymbol{v}^T\boldsymbol{r}=0,\forall \boldsymbol{v}\in \boldsymbol{K}=\boldsymbol{L}, \boldsymbol{r}=\boldsymbol{b}-\boldsymbol{Ax}
$$

`Theorem`: In exact mathematics, CG method converges in at most $m$ steps if $\boldsymbol{A}$ is SPD in $\mathbb{R} ^{m\times m}$.

`Theorem`: In general, CG method converges according to the following estimate ( $\boldsymbol{A}$ is SPD):

$$
\left\| \boldsymbol{x}^{\left( n \right)}-\boldsymbol{x}^{\left( * \right)} \right\| _2\leqslant \left( \frac{\sqrt{\kappa \left( \boldsymbol{A} \right)}-1}{\sqrt{\kappa \left( \boldsymbol{A} \right)}+1} \right) ^n\cdot \left\| \boldsymbol{x}^{\left( 0 \right)}-\boldsymbol{x}^{\left( * \right)} \right\| _2
$$

`Theorem`: The Steepest Descent Method converges with the following estimate ( $\boldsymbol{A}$ is SPD):

$$
\left\| \boldsymbol{x}^{\left( n \right)}-\boldsymbol{x}^{\left( * \right)} \right\| _{\boldsymbol{A}}\leqslant \left( \frac{\kappa \left( \boldsymbol{A} \right) -1}{\kappa \left( \boldsymbol{A} \right) +1} \right) ^n\cdot \left\| \boldsymbol{x}^{\left( 0 \right)}-\boldsymbol{x}^{\left( * \right)} \right\| _{\boldsymbol{A}}
$$

The difference between CG and Steepest Descent Method is significant.

`Example`: If $\kappa \left( \boldsymbol{A} \right) =999$, then:

$$
\frac{\kappa \left( \boldsymbol{A} \right) -1}{\kappa \left( \boldsymbol{A} \right) +1}=0.998,
$$

$$
n=100: \left( 0.998 \right) ^{100}\approx 1-100\times 0.002=0.8;
$$

$$
\frac{\sqrt{\kappa \left( \boldsymbol{A} \right)}-1}{\sqrt{\kappa \left( \boldsymbol{A} \right)}+1}=0.93,
$$

$$
n=0.8100: \left( 0.93 \right) ^{100}\approx 7\times 10^{-4}
$$

**Preconditioning Techniques**: Reduce the condition number of a linear system, so that the convergence is faster.
