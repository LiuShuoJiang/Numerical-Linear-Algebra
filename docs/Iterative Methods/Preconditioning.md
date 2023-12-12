# Preconditioning

## General Introduction

Idea: Perform Preconditioning on Conjugate Gradient method.

**Preconditioning** is to reduce the condition number of the system so that it can converge faster. We want to solve $\boldsymbol{Ax}=\boldsymbol{b}$ where $\kappa \left( \boldsymbol{A} \right) \gg 1$.

If $\boldsymbol{M}$ is invertible, then $\boldsymbol{M}^{-1}\boldsymbol{Ax}=\boldsymbol{M}^{-1}\boldsymbol{b}$. We need to select a good $\boldsymbol{M}$ so that $\kappa \left( \boldsymbol{M}^{-1}\boldsymbol{A} \right) \ll \kappa \left( \boldsymbol{A} \right)$. $\boldsymbol{M}$ is called **preconditioner**.

Here are two criteria that we choose for a good $\boldsymbol{M}$:

1. $\kappa \left( \boldsymbol{M}^{-1}\boldsymbol{A} \right) \ll \kappa \left( \boldsymbol{A} \right)$;
2. $\boldsymbol{M}^{-1}$ is easy to find ( $\boldsymbol{Mz}=\boldsymbol{r}$ is easy to solve).

`Question`: How to apply CG to the preconditioned system $\boldsymbol{M}^{-1}\boldsymbol{Ax}=\boldsymbol{M}^{-1}\boldsymbol{b}$ ? Note that $\boldsymbol{M}^{-1}\boldsymbol{A}$ is no longer guaranteed to be a symmetric matrix.

Assume $\boldsymbol{A}$ is SPD.

## Option 1: Split Preconditioner

Pick a SPD matrix $\boldsymbol{M}$. Using Cholesky factorization: $\boldsymbol{M}=\boldsymbol{LL}^T$. Then we have:

$$
\underbrace{\boldsymbol{L}^{-1}\boldsymbol{A}\left( \boldsymbol{L}^{-1} \right) ^T}\cdot \underset{\boldsymbol{u}}{\underbrace{\boldsymbol{L}^T\boldsymbol{x}}}=\boldsymbol{L}^{-1}\boldsymbol{b}
$$

$$
\Longrightarrow \begin{cases}
	\boldsymbol{L}^T\boldsymbol{x}=\boldsymbol{u}\\
	\underset{\mathrm{SPD}}{\underbrace{\boldsymbol{L}^{-1}\boldsymbol{A}\left( \boldsymbol{L}^{-1} \right) ^T}}\cdot \boldsymbol{u}=\boldsymbol{L}^{-1}\boldsymbol{b}\\
\end{cases}
$$

How to pick $\boldsymbol{L}$ such that $\kappa \left( \boldsymbol{L}^{-1}\boldsymbol{A}\left( \boldsymbol{L}^{-1} \right) ^T \right) \ll \kappa \left( \boldsymbol{A} \right)$ ?

`Strategy`: $\boldsymbol{L}$ is constructed directly from $\boldsymbol{A}$.

`Task` (HW question): Write an algorithm to use CG for $\boldsymbol{L}^{-1}\boldsymbol{A}\left( \boldsymbol{L}^{-1} \right) ^T\boldsymbol{u}=\boldsymbol{L}^{-1}\boldsymbol{b}$.

## Option 2

### Some Concepts

Assume $\boldsymbol{A}$ is SPD. $\boldsymbol{A}^T=\boldsymbol{A}$.

We claim that $\forall \boldsymbol{x},\boldsymbol{y}\in \mathbb{R} ^n$, $\boldsymbol{A}$ is symmetric if $\left< \boldsymbol{x},\boldsymbol{Ay} \right> =\left< \boldsymbol{Ax},\boldsymbol{y} \right>$.

Note that if $\boldsymbol{A}$ is symmetric,

$$
\left< \boldsymbol{x},\boldsymbol{Ay} \right> =\boldsymbol{x}^T\left( \boldsymbol{Ay} \right) =\boldsymbol{x}^T\boldsymbol{A}^T\boldsymbol{y}=\left( \boldsymbol{Ax} \right) ^T\boldsymbol{y}=\left< \boldsymbol{Ax},\boldsymbol{y} \right> 
$$

If $\boldsymbol{M}$ is SPD, we introduce $\boldsymbol{M}$ - **inner product**: $\forall \boldsymbol{x},\boldsymbol{y}\in \mathbb{R} ^n$, define:

$$
\left< \boldsymbol{x},\boldsymbol{y} \right> _{\boldsymbol{M}}\triangleq \left< \boldsymbol{Mx},\boldsymbol{y} \right> =\left< \boldsymbol{x},\boldsymbol{My} \right> 
$$

`Claim`: $\boldsymbol{M}^{-1}\boldsymbol{A}$ is symmetric with respect to $\boldsymbol{M}$ -inner product.

`Verification`: Beacuse $\boldsymbol{A}$ is symmetric,

$$
\left< \left( \boldsymbol{M}^{-1}\boldsymbol{A} \right) \boldsymbol{x}, \boldsymbol{y} \right> _{\boldsymbol{M}}=\left< \boldsymbol{M}\left( \boldsymbol{M}^{-1}\boldsymbol{A} \right) \boldsymbol{x}, \boldsymbol{y} \right> =\left< \boldsymbol{Ax},\boldsymbol{y} \right> 
$$

$$
=\left< \boldsymbol{x},\boldsymbol{Ay} \right> =\left< \boldsymbol{x},\boldsymbol{M}\left( \boldsymbol{M}^{-1}\boldsymbol{A} \right) \boldsymbol{y} \right> =\left< \boldsymbol{x},\left( \boldsymbol{M}^{-1}\boldsymbol{A} \right) \boldsymbol{y} \right> _{\boldsymbol{M}}
$$

### Preconditioned Conjugate Gradients (PCG)

The idea is to replace the standard inner product in CG by the $\boldsymbol{M}$ -inner product for the preconditioned system to **PCG**.

| CG      | PCG |
| :---------------------------:      |    :-------------------------------------:   |
| $\boldsymbol{Ax}=\boldsymbol{b}$                                                                                             | $\boldsymbol{M}^{-1}\boldsymbol{Ax}=\boldsymbol{M}^{-1}\boldsymbol{b}$       |
| $\boldsymbol{r}_j=\boldsymbol{b}-\boldsymbol{Ax}_j$                                                                          |   $\boldsymbol{z}_j=\boldsymbol{M}^{-1}\left( \boldsymbol{b}-\boldsymbol{Ax}_j \right) =\boldsymbol{M}^{-1}\boldsymbol{r}_j,\left( \boldsymbol{Mz}_j=\boldsymbol{r}_j \right)$      |
| $\alpha _j=\frac{\left< \boldsymbol{r}_j,\boldsymbol{r}_j \right>}{\left< \boldsymbol{p}_j,\boldsymbol{Ap}_j \right>}$       | $\alpha _j=\frac{\left< \boldsymbol{z}_j,\boldsymbol{z}_j \right> _{\boldsymbol{M}}}{\left< \boldsymbol{p}_j,\boldsymbol{M}^{-1}\boldsymbol{Ap}_j \right> _{\boldsymbol{M}}}$ |
|  $\tau _j=\frac{\left< \boldsymbol{r}_{j+1},\boldsymbol{r}_{j+1} \right>}{\left< \boldsymbol{r}_j,\boldsymbol{r}_j \right>}$ | $\tau _j=\frac{\left< \boldsymbol{z}_{j+1},\boldsymbol{z}_{j+1} \right> _{\boldsymbol{M}}}{\left< \boldsymbol{z}_j,\boldsymbol{z}_j \right> _{\boldsymbol{M}}}$  |

Also note:

$$
\tau _j=\frac{\left< \boldsymbol{z}_{j+1},\boldsymbol{z}_{j+1} \right> _{\boldsymbol{M}}}{\left< \boldsymbol{z}_j,\boldsymbol{z}_j \right> _{\boldsymbol{M}}}=\frac{\left< \boldsymbol{Mz}_{j+1},\boldsymbol{z}_{j+1} \right>}{\left< \boldsymbol{Mz}_j,\boldsymbol{z}_j \right>}=\frac{\left< \boldsymbol{r}_{j+1},\boldsymbol{z}_{j+1} \right>}{\left< \boldsymbol{r}_j,\boldsymbol{z}_j \right>},
$$

$$
\alpha _j=\frac{\left< \boldsymbol{z}_j,\boldsymbol{z}_j \right> _{\boldsymbol{M}}}{\left< \boldsymbol{p}_j,\boldsymbol{M}^{-1}\boldsymbol{Ap}_j \right> _{\boldsymbol{M}}}=\frac{\left< \boldsymbol{Mz}_j,\boldsymbol{z}_j \right>}{\left< \boldsymbol{p}_j,\boldsymbol{M}\left( \boldsymbol{M}^{-1}\boldsymbol{A} \right) \boldsymbol{p}_j \right>}=\frac{\left< \boldsymbol{r}_j,\boldsymbol{z}_j \right>}{\left< \boldsymbol{p}_j,\boldsymbol{Ap}_j \right>}
$$

`Algorithm` ( **PCG** ):

Compute $\boldsymbol{r}_0=\boldsymbol{b}-\boldsymbol{Ax}^{\left( 0 \right)}$, and set $\boldsymbol{p}_0=\boldsymbol{r}_0$;

Solve $\boldsymbol{Mz}_0=\boldsymbol{r}_0$;

For $j=0,1,2,\cdots$ until convergence, do:

- $\alpha _j=\frac{\left< \boldsymbol{r}_j,\boldsymbol{z}_j \right>}{\left< \boldsymbol{p}_j,\boldsymbol{Ap}_j \right>}$;
- $\boldsymbol{x}^{\left( j+1 \right)}=\boldsymbol{x}^{\left( j \right)}+\alpha _j\boldsymbol{p}_j$;
- $\boldsymbol{r}_{j+1}=\boldsymbol{r}_j-\alpha _j\boldsymbol{Ap}_j$;
- ( *Critical change* ) Solve $\boldsymbol{Mz}_{j+1}=\boldsymbol{r}_{j+1}$;
- $\tau _j=\frac{\left< \boldsymbol{r}_{j+1},\boldsymbol{z}_{j+1} \right>}{\left< \boldsymbol{r}_j,\boldsymbol{z}_j \right>}$;
- $\boldsymbol{p}_{j+1}=\boldsymbol{z}_{j+1}+\tau _j\boldsymbol{p}_j$;

End

`Question`: How to pick $\boldsymbol{M}$ ?

- ***Multigrid***;
- ***Domain Decomposition***.
