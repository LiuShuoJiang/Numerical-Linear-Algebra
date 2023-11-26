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

