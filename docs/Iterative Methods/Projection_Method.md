# Projection Methods

## General framework for Projection Methods

Given two subspaces $\mathbf{L},\mathbf{K}$, and an initial guess $\boldsymbol{x}^{\left( 0 \right)}$, we want to find $\boldsymbol{\delta }\in \mathbf{K}$ such that $\boldsymbol{x}=\boldsymbol{x}^{\left( 0 \right)}+\boldsymbol{\delta }$ produces a ***residual*** $\boldsymbol{r}=\boldsymbol{b}-\boldsymbol{Ax}$ that is orthogonal to $\mathbf{L}$.

$$
\boldsymbol{r}=\boldsymbol{b}-\boldsymbol{Ax}=\boldsymbol{b}-\boldsymbol{A}\left( \boldsymbol{x}^{\left( 0 \right)}+\boldsymbol{\delta } \right) 
$$

$$
=\left( \boldsymbol{b}-\boldsymbol{Ax}^{\left( 0 \right)} \right) -\boldsymbol{A\delta }=\boldsymbol{r}^{\left( 0 \right)}-\boldsymbol{A\delta };
$$

$$
\boldsymbol{r}\bot \mathbf{L}
$$

![Projection](./Projection.png)

Let $\left\{ \boldsymbol{v}_i \right\} _{i=1}^{n}$ be a basis of $\mathbf{K}$, and let $\left\{ \boldsymbol{w}_i \right\} _{i=1}^{n}$ be a basis of $\mathbf{L}$. Note that we assume $\mathrm{dim}\left( \mathbf{K} \right) =\mathrm{dim}\left( \mathbf{L} \right)$. Then:

$$
\boldsymbol{r}\bot \mathbf{L}\ \ \Leftrightarrow \ \ \boldsymbol{r}\bot \boldsymbol{w}_i,\ \left( i=1,2,\cdots ,n \right) ;
$$

$$
\left< \boldsymbol{w}_i,\ \boldsymbol{r} \right> =0,\ \left( i=1,2,\cdots ,n \right) 
$$

Denote:

$$
\boldsymbol{W}=\left[ \begin{matrix}
	\boldsymbol{w}_1&		\cdots&		\boldsymbol{w}_n\\
\end{matrix} \right];
$$

$$
\boldsymbol{V}=\left[ \begin{matrix}
	\boldsymbol{v}_1&		\cdots&		\boldsymbol{v}_n\\
\end{matrix} \right] 
$$


Then:

$$
\boldsymbol{W}^T\boldsymbol{r}=\mathbf{0}
$$

For $\boldsymbol{\delta }\in \mathbf{K}$, we have:

$$
\boldsymbol{\delta }=\sum_{i=1}^n{y_i\boldsymbol{v}_i}=\boldsymbol{Vy};\ \boldsymbol{y}=\left[ \begin{array}{c}
	y_1\\
	\vdots\\
	y_n\\
\end{array} \right] 
$$

We get:

$$
\boldsymbol{x}=\boldsymbol{x}^{\left( 0 \right)}+\boldsymbol{\delta }=\boldsymbol{x}^{\left( 0 \right)}+\boldsymbol{Vy};
$$

$$
\boldsymbol{r}=\boldsymbol{b}-\boldsymbol{Ax}=\boldsymbol{r}^{\left( 0 \right)}-\boldsymbol{A\delta }=\boldsymbol{r}^{\left( 0 \right)}-\boldsymbol{AVy};
$$

$$
\mathbf{0}=\boldsymbol{W}^T\boldsymbol{r}=\boldsymbol{W}^T\left( \boldsymbol{r}^{\left( 0 \right)}-\boldsymbol{AVy} \right) =\boldsymbol{W}^T\boldsymbol{r}^{\left( 0 \right)}-\boldsymbol{W}^T\boldsymbol{AVy};
$$

$$
\boldsymbol{W}^T\boldsymbol{AVy}=\boldsymbol{W}^T\boldsymbol{r}^{\left( 0 \right)}
$$

If $\boldsymbol{W}^T\boldsymbol{AV}$ is invertible, we get:

$$
\boldsymbol{y}=\left( \boldsymbol{W}^T\boldsymbol{AV} \right) ^{-1}\boldsymbol{W}^T\boldsymbol{r}^{\left( 0 \right)}
$$

where $\boldsymbol{W}^T\boldsymbol{AV}$ is a $n\times n$ matrix. Then:

$$
\boldsymbol{x}=\boldsymbol{x}^{\left( 0 \right)}+\boldsymbol{\delta }=\boldsymbol{x}^{\left( 0 \right)}+\boldsymbol{Vy}
$$

$$
=\boldsymbol{x}^{\left( 0 \right)}+\boldsymbol{V}\left( \boldsymbol{W}^T\boldsymbol{AV} \right) ^{-1}\boldsymbol{W}^T\boldsymbol{r}^{\left( 0 \right)}
$$

`Algorithm`( **Projection Algorithm** ):

Until convergence, do:

- Select a pair of subspaces $\mathbf{K},\mathbf{L}$;
- Choose bases $\boldsymbol{V}=\left[ \begin{matrix}
	\boldsymbol{v}_1&		\cdots&		\boldsymbol{v}_n\\
\end{matrix} \right] , \boldsymbol{W}=\left[ \begin{matrix}
	\boldsymbol{w}_1&		\cdots&		\boldsymbol{w}_n\\
\end{matrix} \right]$;
- $\boldsymbol{r}=\boldsymbol{b}-\boldsymbol{Ax}$;
- $\boldsymbol{y}=\left( \boldsymbol{W}^T\boldsymbol{AV} \right) ^{-1}\boldsymbol{W}^T\boldsymbol{r}$;
- $\boldsymbol{x}=\boldsymbol{x}+\boldsymbol{Vy}$;

End

## Examples

### Gauss-Seidel Method

`Example`: Gauss-Seidel is a Projection Method by taking $\mathbf{K}=\mathbf{L}=\mathrm{span}\left\{ \boldsymbol{e}_i \right\}$. Then:

$$
\boldsymbol{\delta }=\left[ \begin{array}{c}
	0\\
	\vdots\\
	0\\
	\delta _i\\
	0\\
	\vdots\\
	0\\
\end{array} \right] , \boldsymbol{x}=\boldsymbol{x}^{\left( 0 \right)}+\boldsymbol{\delta }=\left[ \begin{array}{c}
	{\boldsymbol{x}_1}^{\left( 0 \right)}\\
	\vdots\\
	{\boldsymbol{x}_i}^{\left( 0 \right)}\\
	\vdots\\
	{\boldsymbol{x}_n}^{\left( 0 \right)}\\
\end{array} \right] +\left[ \begin{array}{c}
	0\\
	\vdots\\
	0\\
	\delta _i\\
	0\\
	\vdots\\
	0\\
\end{array} \right] ;
$$

$$
\boldsymbol{r}=\boldsymbol{b}-\boldsymbol{Ax}=\boldsymbol{b}-\boldsymbol{Ax}^{\left( 0 \right)}-\boldsymbol{A\delta }=\boldsymbol{r}^{\left( 0 \right)}-\boldsymbol{A\delta };
$$

$$
\boldsymbol{r}\bot \mathbf{L}\Leftrightarrow \boldsymbol{r}\bot \boldsymbol{e}_i;
$$

$$
0={\boldsymbol{e}_i}^T\boldsymbol{r}=\left( \boldsymbol{b}-\boldsymbol{Ax}^{\left( 0 \right)}-\boldsymbol{A\delta } \right) _i=0;
$$

$$
b_i-\sum_{j=1}^m{a_{ij}{x_j}^{\left( 0 \right)}}-a_{ii}\delta _i=0,
$$

$$
\delta _i=\frac{1}{a_{ii}}\left( b_i-\sum_{j=1}^m{a_{ij}{x_j}^{\left( 0 \right)}} \right) 
$$

### Steepest Descent Method

#### Introduction to Steepest Descent

`Example`: Assume $\boldsymbol{A}$ is SPD, and $\boldsymbol{r}$ is the current residual. We pick $\mathbf{L}=\mathbf{K}=\mathrm{span}\left\{ \boldsymbol{r} \right\}$. Then:

$$
\boldsymbol{\delta }=\alpha \boldsymbol{r}; \boldsymbol{\delta }\in \mathbf{L}, \alpha \in \mathbb{R} 
$$

Also:

$$
\boldsymbol{x}^{new}=\boldsymbol{x}+\alpha \boldsymbol{r},
$$

$$
\boldsymbol{r}^{new}=\boldsymbol{b}-\boldsymbol{Ax}^{new}=\boldsymbol{r}-\alpha \boldsymbol{Ar}
$$

Note that:

$$
\boldsymbol{r}^{new}\bot \boldsymbol{r}\Leftrightarrow \boldsymbol{r}^{new}\bot \mathbf{L},
$$

$$
\left( \boldsymbol{r}^{new} \right) ^T\boldsymbol{r}=0, \left( \boldsymbol{r}-\alpha \boldsymbol{Ar} \right) ^T\boldsymbol{r}=0,
$$

$$
\boldsymbol{r}^T\boldsymbol{r}-\alpha \boldsymbol{r}^T\boldsymbol{Ar}=0
$$

We get:

$$
\alpha =\frac{\boldsymbol{r}^T\boldsymbol{r}}{\boldsymbol{r}^T\boldsymbol{Ar}}, \boldsymbol{x}^{new}=\boldsymbol{x}+\alpha \boldsymbol{r}
$$

`Algorithm`( **Steepest Descent Method** ):

For $k=1,2,\cdots$:

- $\boldsymbol{r}^{\left( k \right)}=\boldsymbol{b}^{\left( k \right)}-\boldsymbol{Ax}^{\left( k \right)}$;
- $\alpha ^{\left( k \right)}=\frac{\left( \boldsymbol{r}^{\left( k \right)} \right) ^T\boldsymbol{r}^{\left( k \right)}}{\left( \boldsymbol{r}^{\left( k \right)} \right) ^T\boldsymbol{Ar}^{\left( k \right)}}$;
- $\boldsymbol{x}^{\left( k+1 \right)}=\boldsymbol{x}^{\left( k \right)}+\alpha ^{\left( k \right)}\boldsymbol{r}^{\left( k \right)}$;

End

#### Discussion on Steepest Descent

`Question`: Why is it called Steepest Descent? Whys is $\alpha$ optimal?

##### First Perspective (from error viewpoint)

Assume $\boldsymbol{x}^*$ is the true solution. We define the error as $\boldsymbol{e}=\boldsymbol{x}-\boldsymbol{x}^*$. Also define $f\left( \boldsymbol{x} \right) =\boldsymbol{e}^T\boldsymbol{Ae}\triangleq \left\| \boldsymbol{e} \right\| _{\boldsymbol{A}}^{2}$ as the $\boldsymbol{A}$-norm.

Note that since $f\left( \boldsymbol{x} \right)$ is convex (quadratic function), there is a unique minimizer, $\boldsymbol{x}^*$ satisfying:

$$
\mathrm{arg}\min_{\boldsymbol{x}} f\left( \boldsymbol{x} \right) =\boldsymbol{x}^*
$$

Using gradient descent to find the minimize search along the direction of negative gradient:

$$
-\nabla f\left( \boldsymbol{x} \right) =-2\boldsymbol{A}\left( \boldsymbol{x}-\boldsymbol{x}^* \right) 
$$

$$
=-2\left( \boldsymbol{Ax}-\boldsymbol{Ax}^* \right) =-2\left( \boldsymbol{Ax}-\boldsymbol{b} \right) 
$$

$$
=2\left( \boldsymbol{b}-\boldsymbol{Ax} \right) =2\boldsymbol{r}
$$

This means that the residual and negative gradient have the same direction.

$$
\boldsymbol{x}^{new}=\boldsymbol{x}+\underset{\mathrm{step}\  \mathrm{size}}{\underbrace{\alpha }}\cdot \boldsymbol{r}
$$

What is the best step size? Actually we would like:

$$
\min_{\alpha} f\left( \boldsymbol{x}+\alpha \boldsymbol{r} \right) 
$$

Therefore, we want to find $\alpha$ such that $\frac{\mathrm{d}}{\mathrm{d}\alpha}\left( f\left( \boldsymbol{x}+\alpha \boldsymbol{r} \right) \right) =0$. It turns out that (HW: Proof):

$$
\alpha =\frac{\boldsymbol{r}^T\boldsymbol{r}}{\boldsymbol{r}^T\boldsymbol{Ar}}
$$

This is the optimal solution for $\alpha$.

##### Second Perspective (Quadratic Optimization)

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

We can examine the level set of $f\left( \boldsymbol{x} \right)$ or $g\left( \boldsymbol{x} \right)$ to get geometric properties.
