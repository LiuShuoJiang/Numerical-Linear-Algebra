# Conjugate Gradient Method (CG)

## Introduction to Conjugate Gradient Algorithm

`Algorithm` ( **Conjugate Gradient Method** ):

Compute $\boldsymbol{r}^{\left( 0 \right)}=\boldsymbol{b}-\boldsymbol{Ax}^{\left( 0 \right)}$, and set $\boldsymbol{p}^{\left( 0 \right)}=\boldsymbol{r}^{\left( 0 \right)}$;

For $j=0,1,\cdots$ until convergence, do:

- Step length: $\alpha _j=\frac{\left< \boldsymbol{r}^{\left( j \right)},\boldsymbol{r}^{\left( j \right)} \right>}{\left< \boldsymbol{p}^{\left( j \right)},\boldsymbol{Ap}^{\left( j \right)} \right>}$;
- Approximate solution: $\boldsymbol{x}^{\left( j+1 \right)}=\boldsymbol{x}^{\left( j \right)}+\alpha _j\boldsymbol{p}^{\left( j \right)}$;
- Residual: $\boldsymbol{r}^{\left( j+1 \right)}=\boldsymbol{r}^{\left( j \right)}-\alpha _j\boldsymbol{Ap}^{\left( j \right)}$;
- Improvement this step: $\tau _j=\frac{\left< \boldsymbol{r}^{\left( j+1 \right)},\boldsymbol{r}^{\left( j+1 \right)} \right>}{\left< \boldsymbol{r}^{\left( j \right)},\boldsymbol{r}^{\left( j \right)} \right>}$;
- Search direction: $\boldsymbol{p}^{\left( j+1 \right)}=\boldsymbol{r}^{\left( j+1 \right)}+\tau _j\boldsymbol{p}^{\left( j \right)}$;

End

The most expensive part for computation of this algorithm is computing $\boldsymbol{Ap}$, whose cost is $O(m^2)$. However, It is needed for only *once* per iteration!

`Question`: How can we derive the CG method?

## Derivation from Projection Method

Perform the projection method with $\mathbf{K}=\mathbf{L}=\mathrm{span}\left\{ \boldsymbol{r}^{\left( j \right)},\boldsymbol{p}^{\left( j-1 \right)} \right\}$:

$$
\boldsymbol{x}^{\left( j+1 \right)}=\boldsymbol{x}^{\left( j \right)}+\boldsymbol{\delta }^{\left( j \right)}
$$

where $\boldsymbol{\delta }^{\left( j \right)}\in \mathbf{K}$. Then:

$$
\boldsymbol{r}^{\left( j+1 \right)}=\boldsymbol{b}-\boldsymbol{Ax}^{\left( j+1 \right)}
$$

Based on the requirement of Projection Method:

$$
\boldsymbol{r}^{\left( j+1 \right)}\bot \mathbf{L}=\mathrm{span}\left\{ \boldsymbol{r}^{\left( j \right)},\boldsymbol{p}^{\left( j-1 \right)} \right\}
$$

$$
\Longrightarrow \left< \boldsymbol{r}^{\left( j+1 \right)},\boldsymbol{r}^{\left( j \right)} \right> =0, \left< \boldsymbol{r}^{\left( j+1 \right)},\boldsymbol{p}^{\left( j-1 \right)} \right> =0
$$

Since $\boldsymbol{\delta }^{\left( j \right)}\in \mathbf{K}$, we can write (in a strange way) $\boldsymbol{\delta }^{\left( j \right)}=\alpha _j\left( \boldsymbol{r}^{\left( j \right)}+\tau _{j-1}\boldsymbol{p}^{\left( j-1 \right)} \right)$, where $\boldsymbol{r}^{\left( j \right)}+\tau _{j-1}\boldsymbol{p}^{\left( j-1 \right)}$ is the *new search direction*.

$$
\boldsymbol{p}^{\left( j \right)}=\boldsymbol{r}^{\left( j \right)}+\tau _{j-1}\boldsymbol{p}^{\left( j-1 \right)}
\\
\Rightarrow \boldsymbol{r}^{\left( j \right)}=\boldsymbol{p}^{\left( j \right)}-\tau _{j-1}\boldsymbol{p}^{\left( j-1 \right)}
$$

where $\alpha _j,\tau _{j-1}$ are to be determined.

`Claim`: We can find (How to prove?):

$$
\boldsymbol{r}^{\left( j+1 \right)}=\boldsymbol{r}^{\left( j \right)}-\alpha _j\boldsymbol{Ap}^{\left( j \right)}
$$

`Proof`:

$$
\boldsymbol{p}^{\left( j \right)}=\boldsymbol{r}^{\left( j \right)}+\tau _{j-1}\boldsymbol{p}^{\left( j-1 \right)},\Leftrightarrow \boldsymbol{r}^{\left( j \right)}=\boldsymbol{p}^{\left( j \right)}-\tau _{j-1}\boldsymbol{p}^{\left( j-1 \right)};
$$

Then:

$$
\boldsymbol{r}^{\left( j+1 \right)}=\boldsymbol{b}-\boldsymbol{Ax}^{\left( j+1 \right)}=\boldsymbol{b}-\boldsymbol{A}\left( \boldsymbol{x}^{\left( j \right)}+\boldsymbol{\delta }^{\left( j \right)} \right) 
$$

$$
=\left( \boldsymbol{b}-\boldsymbol{Ax}^{\left( j \right)} \right) -\boldsymbol{A\delta }^{\left( j \right)}=\boldsymbol{r}^{\left( j \right)}-\boldsymbol{A}\delta ^{\left( j \right)}
$$

$$
=\boldsymbol{r}^{\left( j \right)}-\boldsymbol{A}\alpha _j\left( \boldsymbol{r}^{\left( j \right)}+\tau _{j-1}\boldsymbol{p}^{\left( j-1 \right)} \right) 
$$

$$
=\boldsymbol{r}^{\left( j \right)}-\alpha _j\boldsymbol{A}\left( \boldsymbol{p}^{\left( j \right)}-\tau _{j-1}\boldsymbol{p}^{\left( j-1 \right)} \right) -\alpha _j\tau _{j-1}\boldsymbol{Ap}^{\left( j-1 \right)}
$$

$$
=\boldsymbol{r}^{\left( j \right)}-\alpha _j\boldsymbol{Ap}^{\left( j \right)}
$$

End of proof.

Because:

$$
0=\left< \boldsymbol{r}^{\left( j+1 \right)},\boldsymbol{r}^{\left( j \right)} \right> =\left< \boldsymbol{r}^{\left( j \right)}-\alpha _j\boldsymbol{Ap}^{\left( j \right)},\boldsymbol{r}^{\left( j \right)} \right> =\left< \boldsymbol{r}^{\left( j \right)},\boldsymbol{r}^{\left( j \right)} \right> -\alpha _j\left< \boldsymbol{r}^{\left( j \right)},\boldsymbol{Ap}^{\left( j \right)} \right> 
$$

We have:

$$
\alpha _j=\frac{\left< \boldsymbol{r}^{\left( j \right)},\boldsymbol{r}^{\left( j \right)} \right>}{\left< \boldsymbol{r}^{\left( j \right)},\boldsymbol{Ap}^{\left( j \right)} \right>}
$$

Also:

$$
0=\left< \boldsymbol{r}^{\left( j+1 \right)},\boldsymbol{p}^{\left( j-1 \right)} \right> \Rightarrow \tau _{j-1}=\frac{\left< \boldsymbol{r}^{\left( j \right)},\boldsymbol{r}^{\left( j \right)} \right>}{\left< \boldsymbol{r}^{\left( j-1 \right)},\boldsymbol{r}^{\left( j-1 \right)} \right>}
$$

## Properties of Conjugate Gradient Method

It can be proved that CG satisfies the following properties:

- $\left< \boldsymbol{r}^{\left( j \right)},\boldsymbol{r}^{\left( i \right)} \right> =0, i\ne j;$
- $\left< \boldsymbol{Ap}^{\left( j \right)},\boldsymbol{p}^{\left( i \right)} \right> =0,i\ne j$.

Note that for $\alpha _j$, the formulas just derived here are slightly different from those described in the algorithm, but they are in fact equivalent. We can show:

$$
\alpha _j=\frac{\left< \boldsymbol{r}^{\left( j \right)},\boldsymbol{r}^{\left( j \right)} \right>}{\left< \boldsymbol{r}^{\left( j \right)},\boldsymbol{Ap}^{\left( j \right)} \right>}=\frac{\left< \boldsymbol{r}^{\left( j \right)},\boldsymbol{r}^{\left( j \right)} \right>}{\left< \boldsymbol{p}^{\left( j \right)},\boldsymbol{Ap}^{\left( j \right)} \right>},
$$

$$
\Longleftarrow \left< \boldsymbol{r}^{\left( j \right)},\boldsymbol{Ap}^{\left( j \right)} \right> =\left< \boldsymbol{p}^{\left( j \right)},\boldsymbol{Ap}^{\left( j \right)} \right> ;
$$

$$
\left< \boldsymbol{r}^{\left( j \right)},\boldsymbol{Ap}^{\left( j \right)} \right> =\left< \boldsymbol{p}^{\left( j \right)}-\tau _{j-1}\boldsymbol{p}^{\left( j-1 \right)},\boldsymbol{Ap}^{\left( j \right)} \right> 
$$

$$
=\left< \boldsymbol{p}^{\left( j \right)},\boldsymbol{Ap}^{\left( j \right)} \right> -\tau _{j-1}\left< \boldsymbol{p}^{\left( j-1 \right)},\boldsymbol{Ap}^{\left( j \right)} \right> 
$$

$$
=\left< \boldsymbol{p}^{\left( j \right)},\boldsymbol{Ap}^{\left( j \right)} \right> 
$$

In this way we prove the equivalence of formulas.
