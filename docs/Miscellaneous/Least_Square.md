# Least Square Problems

Given $\boldsymbol{A}\in \mathbb{R} ^{m\times n}, m\geqslant n; \boldsymbol{b}\in \mathbb{R} ^m$, the **Least Square** problem is defined as: Find $\boldsymbol{x}^*\in \mathbb{R} ^n$, such that:

$$
\boldsymbol{x}^*=\mathrm{arg}\min_{\boldsymbol{x}} \left\| \boldsymbol{b}-\boldsymbol{Ax} \right\| _2
$$

where $\boldsymbol{r}=\boldsymbol{b}-\boldsymbol{Ax}$ is called **residual**. Residual $\boldsymbol{r}\ne \mathbf{0}$.

We will introduce three different methods to solve the Least Square problem.

## Normal Equation

### Traditional Method

We solve $\boldsymbol{A}^*\boldsymbol{Ax}=\boldsymbol{A}^*\boldsymbol{b}$, where $\boldsymbol{A}^*\boldsymbol{A}\in \mathbb{R} ^{n\times n}$. If $\boldsymbol{A}$ is full rank, then $\boldsymbol{A}^*\boldsymbol{A}$ is invertible.

The Least Square solution is:

$$
\boldsymbol{x}^*=\left( \boldsymbol{A}^*\boldsymbol{A} \right) ^{-1}\boldsymbol{A}^*\boldsymbol{b}
$$

`Question`: How about $\boldsymbol{A}$ is not full rank?

`Comment`: This method is not recommended is $\kappa \left( \boldsymbol{A} \right)$ is large.

Note: $\kappa \left( \boldsymbol{A}^*\boldsymbol{A} \right) =\left( \kappa \left( \boldsymbol{A} \right) \right) ^2$

Usually in practice, the number of matrix rows are way larger than that of columns.

## Two Recommended Methods

### QR Factorization

$$
\boldsymbol{A}=\boldsymbol{QR}\,\,\Longleftrightarrow \boldsymbol{Q}^*\boldsymbol{A}=\boldsymbol{R}
$$

$$
\left\| \boldsymbol{Ax}-\boldsymbol{b} \right\| _2=\left\| \boldsymbol{Q}^*\left( \boldsymbol{Ax}-\boldsymbol{b} \right) \right\| _2
$$

$$
=\left\| \boldsymbol{Q}^*\boldsymbol{Ax}-\boldsymbol{Q}^*\boldsymbol{b} \right\| _2=\left\| \boldsymbol{Rx}-\boldsymbol{Q}^*\boldsymbol{b} \right\| _2
$$

Let $\boldsymbol{c}=\boldsymbol{Q}^*\boldsymbol{b}$, then:

$$
\boldsymbol{x}^*=\boldsymbol{R}^{-1}\boldsymbol{c}=\boldsymbol{R}^{-1}\boldsymbol{Q}^*\boldsymbol{b}
$$

### SVD Method (Recommended)

Assume $\boldsymbol{A}=\boldsymbol{U}\mathbf{\Sigma }\boldsymbol{V}^*$ as *full SVD*. We can get:

$$
\boldsymbol{A}=\boldsymbol{U}\mathbf{\Sigma }\boldsymbol{V}^*\Leftrightarrow \boldsymbol{U}^*\boldsymbol{AV}=\mathbf{\Sigma };
$$

$$
\left\| \boldsymbol{Ax}-\boldsymbol{b} \right\| _2=\left\| \boldsymbol{U}^*\boldsymbol{AVV}^*\boldsymbol{x}-\boldsymbol{U}^*\boldsymbol{b} \right\| _2=\left\| \mathbf{\Sigma }\boldsymbol{V}^*\boldsymbol{x}-\boldsymbol{U}^*\boldsymbol{b} \right\| _2
$$

Define $\boldsymbol{V}^*\boldsymbol{x}=\boldsymbol{y},\ \boldsymbol{w}=\boldsymbol{U}^*\boldsymbol{b}$, then:

$$
\left\| \boldsymbol{Ax}-\boldsymbol{b} \right\| _2=\left\| \mathbf{\Sigma }\boldsymbol{y}-\boldsymbol{w} \right\| _2;
$$

$$
\mathrm{arg}\min_{\boldsymbol{y}} \left\| \mathbf{\Sigma }\boldsymbol{y}-\boldsymbol{w} \right\| _2
$$

$$
\Longrightarrow \boldsymbol{y}^*=\mathbf{\Sigma }^{-1}\boldsymbol{w}=\mathbf{\Sigma }^{-1}\boldsymbol{U}^*\boldsymbol{b},
$$

$$
\boldsymbol{x}^*=\boldsymbol{Vy}^*=\boldsymbol{V}\mathbf{\Sigma }^{-1}\boldsymbol{U}^*\boldsymbol{b}
$$

We get the formula:

$$
\boldsymbol{x}^*=\boldsymbol{V}\mathbf{\Sigma }^{-1}\boldsymbol{U}^*\boldsymbol{b}
$$

However,  $\mathbf{\Sigma }^{-1}$ may not exist! Assume:

$$
\mathbf{\Sigma }=\left[ \begin{matrix}
	\sigma _1&		&		&		&		&		\\
	&		\ddots&		&		&		&		\\
	&		&		\sigma _r&		&		&		\\
	&		&		&		0&		&		\\
	&		&		&		&		\ddots&		\\
	&		&		&		&		&		0\\
\end{matrix} \right] , \sigma _1\geqslant \sigma _2\geqslant \cdots \geqslant \sigma _r>\sigma _{r+1}=\cdots =\sigma _m=0
$$

The **pseudo inverse** $\mathbf{\Sigma }^{-1}$ is defined as:

$$
\mathbf{\Sigma }^{-1}=\left[ \begin{matrix}
	\frac{1}{\sigma _1}&		&		&		&		&		\\
	&		\ddots&		&		&		&		\\
	&		&		\frac{1}{\sigma _r}&		&		&		\\
	&		&		&		0&		&		\\
	&		&		&		&		\ddots&		\\
	&		&		&		&		&		0\\
\end{matrix} \right] 
$$

Then $\boldsymbol{x}^*=\boldsymbol{V}\mathbf{\Sigma }^{-1}\boldsymbol{U}^*\boldsymbol{b}$ can be applied to the reduced SVD.

We can rewrite another formula of the result:

$$
\boldsymbol{x}^*=\sum_{i=1}^r{\frac{{\boldsymbol{u}_i}^T\boldsymbol{b}}{\sigma _i}\boldsymbol{v}_i}
$$

where $\boldsymbol{u}_i, \boldsymbol{v}_i$ are columns of $\boldsymbol{U}, \boldsymbol{V}$ respectively.

`Comments`: If $\sigma _i$ is small, the error is enlarged by $\sigma _i$.

We can compute the approximate value ("the **best approximation** of $\boldsymbol{x}^*$ in $\mathbb{R} ^k$ "):

$$
\boldsymbol{y}^*=\sum_{i=1}^k{\frac{{\boldsymbol{u}_i}^T\boldsymbol{b}}{\sigma _i}\boldsymbol{v}_i}, k<r
$$

This is the idea behind [Principal Component Analysis](https://en.wikipedia.org/wiki/Principal_component_analysis) (where to cut off?).
