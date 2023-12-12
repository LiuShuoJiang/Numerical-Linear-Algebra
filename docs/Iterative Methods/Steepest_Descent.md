# Steepest Descent Method

## Introduction to Steepest Descent

`Example`: Assume $\boldsymbol{A}$ is SPD, and $\boldsymbol{r}$ is the current residual. We pick $\mathbf{L}=\mathbf{K}=\mathrm{span}\left\{ \boldsymbol{r} \right\}$ if $\boldsymbol{r} \ne 0$. Then:

$$
\boldsymbol{\delta }=\alpha \boldsymbol{r}; \boldsymbol{\delta }\in \mathbf{L}, \alpha \in \mathbb{R} 
$$

Also:

$$
\boldsymbol{x}^{\mathrm{new}}=\boldsymbol{x}+\alpha \boldsymbol{r},
$$

$$
\boldsymbol{r}^{\mathrm{new}}=\boldsymbol{b}-\boldsymbol{Ax}^{\mathrm{new}}=\boldsymbol{r}-\alpha \boldsymbol{Ar}
$$

Note that:

$$
\boldsymbol{r}^{\mathrm{new}}\bot \boldsymbol{r}\Leftrightarrow \boldsymbol{r}^{\mathrm{new}}\bot \mathbf{L},
$$

$$
\left( \boldsymbol{r}^{\mathrm{new}} \right) ^T\boldsymbol{r}=0,\left( \boldsymbol{r}-\alpha \boldsymbol{Ar} \right) ^T\boldsymbol{r}=0,
$$

$$
\Longrightarrow \boldsymbol{r}^T\boldsymbol{r}-\alpha \boldsymbol{r}^T\boldsymbol{Ar}=0
$$

We get:

$$
\alpha =\frac{\boldsymbol{r}^T\boldsymbol{r}}{\boldsymbol{r}^T\boldsymbol{Ar}}, \boldsymbol{x}^{\mathrm{new}}=\boldsymbol{x}+\alpha \boldsymbol{r}
$$

`Algorithm`( **Steepest Descent Method** ):

For $k=1,2,\cdots$:

- $\boldsymbol{r}^{\left( k \right)}=\boldsymbol{b}^{\left( k \right)}-\boldsymbol{Ax}^{\left( k \right)}$;
- $\alpha ^{\left( k \right)}=\frac{\left( \boldsymbol{r}^{\left( k \right)} \right) ^T\boldsymbol{r}^{\left( k \right)}}{\left( \boldsymbol{r}^{\left( k \right)} \right) ^T\boldsymbol{Ar}^{\left( k \right)}}$;
- $\boldsymbol{x}^{\left( k+1 \right)}=\boldsymbol{x}^{\left( k \right)}+\alpha ^{\left( k \right)}\boldsymbol{r}^{\left( k \right)}$;

End

## Discussion on Steepest Descent

`Question`: Why is it called "Steepest Descent"? Why is $\alpha$ optimal?

### First Perspective (From Error Viewpoint)

Assume $\boldsymbol{x}^*$ is the true solution that we want to compute. We define the error as $\boldsymbol{e}=\boldsymbol{x}-\boldsymbol{x}^*$. Also define $f\left( \boldsymbol{x} \right) =\boldsymbol{e}^T\boldsymbol{Ae}\triangleq \left\| \boldsymbol{e} \right\| _{\boldsymbol{A}}^{2}$ as the $\boldsymbol{A}$-norm ( $\boldsymbol{A}$ is SPD here).

Note that since $f\left( \boldsymbol{x} \right)$ is *convex* (quadratic function), there is a unique minimizer, $\boldsymbol{x}^*$ satisfying:

$$
\mathrm{arg}\min_{\boldsymbol{x}} f\left( \boldsymbol{x} \right) =\boldsymbol{x}^*
$$

We use gradient descent to find the minimizing search along the direction of *negative gradient*:

$$
-\nabla f\left( \boldsymbol{x} \right) =-2\boldsymbol{A}\left( \boldsymbol{x}-\boldsymbol{x}^* \right) 
$$

$$
=-2\left( \boldsymbol{Ax}-\boldsymbol{Ax}^* \right) =-2\left( \boldsymbol{Ax}-\boldsymbol{b} \right) 
$$

$$
=2\left( \boldsymbol{b}-\boldsymbol{Ax} \right) =2\boldsymbol{r}
$$

This means that **the residual and negative gradient have the same direction**.

$$
\boldsymbol{x}^{\mathrm{new}}=\boldsymbol{x}+\underset{\mathrm{step}\ \mathrm{size}}{\underbrace{\alpha }}\cdot \boldsymbol{r}
$$

What is the best step size? Actually we would like:

$$
\min_{\alpha} f\left( \boldsymbol{x}+\alpha \boldsymbol{r} \right) 
$$

Therefore, we want to find $\alpha$ such that $\frac{\mathrm{d}}{\mathrm{d}\alpha}\left( f\left( \boldsymbol{x}+\alpha \boldsymbol{r} \right) \right) =0$. It turns out that (how to prove?):

$$
\alpha =\frac{\boldsymbol{r}^T\boldsymbol{r}}{\boldsymbol{r}^T\boldsymbol{Ar}}
$$

This is the optimal solution for $\alpha$.

### Second Perspective (Quadratic Optimization)

We define:

$$
g\left( \boldsymbol{x} \right) =\frac{1}{2}\boldsymbol{x}^T\boldsymbol{Ax}-\boldsymbol{b}^T\boldsymbol{x}
$$

Note that:

$$
\boldsymbol{x}^*=\mathrm{arg}\min_{\boldsymbol{x}} g\left( \boldsymbol{x} \right) 
$$

It turns out that:

$$
-\nabla g\left( \boldsymbol{x} \right) =-\left( \boldsymbol{Ax}-\boldsymbol{b} \right) =\boldsymbol{r}
$$

Similarly, we get:

$$
\alpha =\frac{\boldsymbol{r}^T\boldsymbol{r}}{\boldsymbol{r}^T\boldsymbol{Ar}}
$$

is the optimal choice along $\boldsymbol{r}$.

We can examine the [level set](https://en.wikipedia.org/wiki/Level_set) of $f\left( \boldsymbol{x} \right)$ or $g\left( \boldsymbol{x} \right)$ to get geometric properties.

## Further Remarks on Steepest Descent

Consider two cases for $\boldsymbol{A}$:

- `Case 1`: Eigenvalues $\frac{\lambda _{\max}}{\lambda _{\min}}\sim O\left( 1 \right)$ (nearly a circle, *fast* convergence);
- `Case 2`: Eigenvalues $\frac{\lambda _{\max}}{\lambda _{\min}}\gg 1$ (nearly a very flat oval, *slow* convergence).

`Conclusion`: The **condition number** of $\boldsymbol{A}$ determines the speed of convergence:

- If $\kappa \left( \boldsymbol{A} \right) \gg 1$, we call $\boldsymbol{A}$ ***ill-conditioned***;
- Otherwise, we call $\boldsymbol{A}$ ***well-conditioned***.
