# Generalized Minimal Residual Method (GMRES)

The **Generalized Minimal Residual Method (GMRES)** is an iterative method for the numerical solution of an indefinite non-symmetric system of linear equations.

## Arnoldi Method

We apply Krylov Subspace based Projection Method to non-symmetric matrices. Pick:

$$
\boldsymbol{A}\mathbf{K}=\mathbf{L}
$$

where $\mathbf{K}=\mathbf{K}_n$ is Krylov Subspace:

$$
\mathbf{K}_n=\mathrm{span}\left\{ \boldsymbol{r}^{\left( 0 \right)},\boldsymbol{Ar}^{\left( 0 \right)},\boldsymbol{A}^2\boldsymbol{r}^{\left( 0 \right)},\cdots ,\boldsymbol{A}^{n-1}\boldsymbol{r}^{\left( 0 \right)} \right\} 
$$

and:

$$
\mathbf{L}=\boldsymbol{A}\mathbf{K}_n=\mathrm{span}\left\{ \boldsymbol{Ar}^{\left( 0 \right)},\boldsymbol{A}^2\boldsymbol{r}^{\left( 0 \right)},\boldsymbol{A}^3\boldsymbol{r}^{\left( 0 \right)},\cdots ,\boldsymbol{A}^n\boldsymbol{r}^{\left( 0 \right)} \right\} 
$$

Now apply ***modified Gram-Schmidt*** to Krylov Subspace. If $\boldsymbol{A}$ is non-symmetric, this process is called **Arnoldi Procedure**. If $\boldsymbol{A}$ is SPD, this process is called **Lancozs Method**.

`Algorithm` ( **Arnoldi Procedure** ):

$\boldsymbol{r}^{\left( 0 \right)}\ne 0$ is an arbitrary vector;

$\boldsymbol{q}_1=\frac{\boldsymbol{r}^{\left( 0 \right)}}{\left\| \boldsymbol{r}^{\left( 0 \right)} \right\| _2}$;

For $n=1,2,3,\cdots$:

- $\boldsymbol{v}=\boldsymbol{Aq}_n$;
- For $j=1,2,\cdots ,n$:
    - $h_{jn}=\left< \boldsymbol{q}_j,\boldsymbol{v} \right>$;
    - $\boldsymbol{v}=\boldsymbol{v}-h_{jn}\boldsymbol{q}_j$;
- End;
- $h_{n+1,n}=\left\| \boldsymbol{v} \right\|$;
- $\boldsymbol{q}_{n+1}=\frac{\boldsymbol{v}}{h_{n+1,n}}$;

End

Arnoldi Procedure creates:

$$
\boldsymbol{Q}_n=\left[ \begin{matrix}
	\boldsymbol{q}_1&		\boldsymbol{q}_2&		\cdots&		\boldsymbol{q}_n\\
\end{matrix} \right] , \boldsymbol{H}_n=\left[ h_{ij} \right] ;
$$

$$
\Longrightarrow {\boldsymbol{Q}_n}^T\boldsymbol{AQ}_n=\boldsymbol{H}_n
$$

Also denote:

$$
\overline{\boldsymbol{H}}_n=\left[ \begin{matrix}
	\boldsymbol{H}_n&		\boldsymbol{O}\\
	\boldsymbol{O}&		h_{n+1,n}\\
\end{matrix} \right] 
$$

Then we can verify:

$$
\boldsymbol{AQ}_n=\boldsymbol{Q}_{n+1}\overline{\boldsymbol{H}}_n
$$

## Introduction to GMRES

### Derivation

`Theorem`: Let $\boldsymbol{A}$ be a square nonsingular matrix. $\mathbf{L}=\boldsymbol{A}\mathbf{K}$. The vector $\boldsymbol{x}\prime$ is the result of the Projection Method onto $\mathbf{K}$, orthogonal to $\mathbf{L}$ with the starting point $\boldsymbol{x}^{\left( 0 \right)}$ *if and only if* $\boldsymbol{x}\prime$ minimizes the 2-norm of the residual vector over $\boldsymbol{x}\in \boldsymbol{x}^{\left( 0 \right)}+\mathbf{K}$. i.e., define:

$$
R\left( \boldsymbol{x} \right) =\left\| \boldsymbol{b}-\boldsymbol{Ax} \right\| _2
$$

$$
\Rightarrow \boldsymbol{x}\prime=\mathrm{arg}\min R\left( \boldsymbol{x} \right) , \boldsymbol{x}\in \boldsymbol{x}^{\left( 0 \right)}+\mathbf{K}
$$

Based on the argument above, we can do a quick derivation of GMRES. Now select:

$$
\mathbf{K}=\mathbf{K}_n,\mathbf{L}=\boldsymbol{A}\mathbf{K}_n
$$

And $\left\{ \boldsymbol{q}_1,\cdots ,\boldsymbol{q}_n \right\}$ is an orthonormal basis for $\mathbf{K}$. Also $\boldsymbol{r}^{\left( n \right)}=\boldsymbol{b}-\boldsymbol{Ax}^{\left( n \right)}$.

The solution of Projection Method:

$$
\boldsymbol{x}^{\left( n \right)}=\mathrm{arg}\min_{\boldsymbol{x}\in \boldsymbol{x}^{\left( 0 \right)}+\mathbf{K}} \left\| \boldsymbol{b}-\boldsymbol{Ax} \right\| _2
$$

$$
\boldsymbol{r}^{\left( n \right)}\bot \mathbf{L},
$$

$$
\boldsymbol{x}^{\left( n \right)}\in \boldsymbol{x}^{\left( 0 \right)}+\mathbf{K}\Longrightarrow \boldsymbol{x}^{\left( n \right)}=\boldsymbol{x}^{\left( 0 \right)}+\boldsymbol{Q}_n\boldsymbol{y},\ \boldsymbol{y}\in \mathbb{R} ^n
$$

Then:

$$
\left\| \boldsymbol{b}-\boldsymbol{Ax} \right\| _2=\left\| \boldsymbol{b}-\boldsymbol{A}\left( \boldsymbol{x}^{\left( 0 \right)}+\boldsymbol{Q}_n\boldsymbol{y} \right) \right\| _2
$$

$$
=\left\| \left( \boldsymbol{b}-\boldsymbol{Ax}^{\left( 0 \right)} \right) -\boldsymbol{AQ}_n\boldsymbol{y} \right\| _2=\left\| \boldsymbol{r}^{\left( 0 \right)}-\boldsymbol{AQ}_n\boldsymbol{y} \right\| _2
$$

$$
=\left\| \beta \boldsymbol{q}_1-\boldsymbol{Q}_{n+1}\overline{\boldsymbol{H}}_n\boldsymbol{y} \right\| _2
$$

$$
=\left\| \boldsymbol{Q}_{n+1}\left( \beta \boldsymbol{e}_1-\overline{\boldsymbol{H}}_n\boldsymbol{y} \right) \right\| _2=\left\| \beta \boldsymbol{e}_1-\overline{\boldsymbol{H}}_n\boldsymbol{y} \right\| _2
$$

We can find:

$$
\min_{\boldsymbol{x}\in \boldsymbol{x}^{\left( 0 \right)}+\mathbf{K}} \left\| \boldsymbol{b}-\boldsymbol{Ax} \right\| _2=\min_{\boldsymbol{y}\in \mathbb{R} ^n} \left\| \beta \boldsymbol{e}_1-\overline{\boldsymbol{H}}_n\boldsymbol{y} \right\| _2
$$

Therefore, the idea of GMRES algorithm is to find the least square solution to:

$$
\min_{\boldsymbol{y}\in \mathbb{R} ^n} \left\| \beta \boldsymbol{e}_1-\overline{\boldsymbol{H}}_n\boldsymbol{y} \right\| _2
$$

for $\boldsymbol{y}$, and update:

$$
\boldsymbol{x}^{\left( n \right)}=\boldsymbol{x}^{\left( 0 \right)}+\boldsymbol{Q}_n\boldsymbol{y}
$$

### Algorithm Steps

`Algorithm` ( **GMRES** ):

Compute $\boldsymbol{r}^{\left( 0 \right)}=\boldsymbol{b}-\boldsymbol{Ax}^{\left( 0 \right)}, \beta =\left\| \boldsymbol{r}^{\left( 0 \right)} \right\| , \boldsymbol{q}_1=\frac{\boldsymbol{r}^{\left( 0 \right)}}{\beta}$;

Define a $(n+1)\times n$ matrix $\overline{\boldsymbol{H}}_n=\left\{ h_{ij} \right\} _{\left( n+1 \right) \times n}$, and set $\overline{\boldsymbol{H}}_n=\boldsymbol{O}$;

For $j=1,2,\cdots ,n$:

- $\boldsymbol{w}_j=\boldsymbol{Aq}_j$;
- For $i=1,2,\cdots ,j$:
    - $h_{ij}=\left< \boldsymbol{w}_j,\boldsymbol{q}_i \right>$;
    - $\boldsymbol{w}_j=\boldsymbol{w}_j-h_{ij}\boldsymbol{q}_i$;
- End;
- $h_{j+1,j}=\left\| \boldsymbol{w}_j \right\| _2$;
- $\boldsymbol{q}_{j+1}=\frac{\boldsymbol{w}_j}{h_{j+1,j}}$;

End;

Compute $\boldsymbol{y}_n$ to minimize $\left\| \beta \boldsymbol{e}_1-\overline{\boldsymbol{H}}_n\boldsymbol{y} \right\| _2$;

Set $\boldsymbol{x}^{\left( n \right)}=\boldsymbol{x}^{\left( 0 \right)}+\boldsymbol{Q}_n\boldsymbol{y}_n$;

End of procedure.

`Note`: If $h_{j+1,j} =0$, break! (This is a lucky breakdown.) In this case, the solution is in $\boldsymbol{x}^{\left( 0 \right)}+\mathbf{K}$.

## Further Discussion

In practice, ***restarting*** GMRES can avoid large $n$.

Theoretically, GMRES converges quickly if $\kappa \left( \boldsymbol{V} \right)$ is not loo large, where $\boldsymbol{V}$ is the eigenmatrix of $\boldsymbol{A}$.

*Bi-CG*, *CGN* are other similar iterative methods.
