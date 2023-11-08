# Linear Algebra Review

## Vector Norm

The range of vector norm is $\boldsymbol{x}\in \mathbb{R} ^n, \left\| \boldsymbol{x} \right\| \in \mathbb{R} ^+\cup \left\{ \mathbf{0} \right\}$.

Some examples:

$L^2$-**Norm**:

$$
\left\| \boldsymbol{x} \right\| _2=\sqrt{\boldsymbol{x}^T\boldsymbol{x}}=\left< \boldsymbol{x},\boldsymbol{x} \right> ^{\frac{1}{2}}
$$

$L^1$-**Norm**:

$$
\left\| \boldsymbol{x} \right\| _1=\sum_{i=1}^n{\left| x_i \right|}
$$

$L^\infty$-**Norm**:

$$
\left\| \boldsymbol{x} \right\| _{\infty}=\max_{1\le i\le n} \left\| x_i \right\| 
$$

$L^p$-**Norm**:

$$
\left\| \boldsymbol{x} \right\| _p=\left( \sum_{i=1}^n{\left| x_i \right|^p} \right) ^{\frac{1}{p}}, p\in \left[ 1,\infty \right) 
$$

**Weighted Norm**: ( $\boldsymbol{W}$ is a given matrix)

$$
\boldsymbol{W}=\left[ \begin{matrix}
	w_1&		\cdots&		O\\
	\vdots&		\ddots&		\vdots\\
	O&		\cdots&		w_n\\
\end{matrix} \right] ,
\\
\left\| \boldsymbol{x} \right\| _{\boldsymbol{w}}=\left\| \boldsymbol{Wx} \right\| _2=\left( \sum_{i=1}^n{\left| w_ix_i \right|^2} \right) ^{\frac{1}{2}}
$$

Geometric Intuition:

Unit disk with different norms: $\left\{ \boldsymbol{x}\in \mathbb{R} ^n|\left\| \boldsymbol{x} \right\| \le 1 \right\}$

$L^2$ is like a circle. $L^1$ is like a diamond. $L^\infty$ is like a square. $L^p (1<p<2)$ is like the shape between circle and diamond. $L^p (p>2)$ is like the shape between circle and square.

## Matrix Norm

Two different ways to define the matrix norm for $\boldsymbol{A}\in \mathbb{R} ^{m\times n}$:

*Way 1*: View matrix as a vector (reshape):

$$
\left\| \boldsymbol{A} \right\| _F=\left( \sum_{i=1}^m{\sum_{j=1}^n{\left| a_{ij} \right|^2}} \right) ^{\frac{1}{2}}
$$

*Way 2*: **Induced Matrix Norm** (preferred way to define the matrix norm)

$$
\boldsymbol{A}\in \mathbb{R} ^{m\times n}, \boldsymbol{A}:\mathbb{R} ^n\longmapsto \mathbb{R} ^m\,\,\left( \forall \boldsymbol{x}\in \mathbb{R} ^n, \boldsymbol{Ax}\in \mathbb{R} ^m \right) 
$$

$$
\left\| \boldsymbol{A} \right\| _{\left( m,n \right)}=\mathop {\mathrm{sup}} \limits_{\left\| \boldsymbol{x} \right\| \ne 0}\frac{\left\| \boldsymbol{Ax} \right\| _v}{\left\| \boldsymbol{x} \right\| _v}
$$

Also:

$$
\left\| \boldsymbol{A} \right\| _{\left( m,n \right)}=\mathop {\mathrm{sup}} \limits_{\begin{array}{c}
	\boldsymbol{x}\in \mathbb{R} ^n,\\
	\left\| \boldsymbol{x} \right\| =1\\
\end{array}}\frac{\left\| \boldsymbol{Ax} \right\|}{\left\| \boldsymbol{x} \right\|}=\mathop {\mathrm{sup}} \limits_{\begin{array}{c}
	\boldsymbol{x}\in \mathbb{R} ^n,\\
	\left\| \boldsymbol{x} \right\| =1\\
\end{array}}\left\| \boldsymbol{Ax} \right\| 
$$

`Conclusion`:

If $\left\| \boldsymbol{x} \right\| _1$ is taken:

$$
\boldsymbol{A}=\left[ \begin{matrix}
	\boldsymbol{a}_1&		\boldsymbol{a}_2&		\cdots&		\boldsymbol{a}_n\\
\end{matrix} \right] , \left\| \boldsymbol{A} \right\| _1=\max_{1\le j\le n} \left\| \boldsymbol{a}_j \right\| _1
$$

If $\left\| \boldsymbol{x} \right\| _{\infty}$ is taken:

$$
\boldsymbol{A}=\left[ \begin{array}{c}
	{\boldsymbol{a}_1}^*\\
	{\boldsymbol{a}_2}^*\\
	\vdots\\
	{\boldsymbol{a}_m}^*\\
\end{array} \right] , \left\| \boldsymbol{A} \right\| _{\infty}=\max_{1\le i\le m} \left\| {\boldsymbol{a}_i}^* \right\| _1
$$

What about induced 2-norms? It turns out that we cannot find a formula for it!

The 2-norm(spectral norm) of a matrix $\boldsymbol{A}$ is the largest singular value of $\boldsymbol{A}$ (i.e., the square root of the largest eigenvalue of the matrix $\boldsymbol{A}^{*}\boldsymbol{A}$, where $\boldsymbol{A}^{*}$ denotes the conjugate transpose of $\boldsymbol{A}$ ):

$$
\left\| \boldsymbol{A} \right\| _2=\sqrt{\lambda _{\max}\left( \boldsymbol{A}^*\boldsymbol{A} \right)}=\sigma _{\max}\left( \boldsymbol{A} \right) 
$$

## Spectral Radius

The **Spectral Radius** of a square matrix is the maximum of the absolute values of its eigenvalues.

Let $\lambda _1,\lambda _2,\cdots ,\lambda _n$ be the eigenvalues of a matrix $\boldsymbol{A}\in \mathbb{C} ^{n\times n}$. Then the spectral radius of $\boldsymbol{A}$ is defined as:

$$
\rho \left( \boldsymbol{A} \right) =\max \left( \left| \lambda _1 \right|,\left| \lambda _2 \right|,\cdots ,\left| \lambda _n \right| \right) 
$$

If $\boldsymbol{A}$ satisfies $\boldsymbol{AA}^*=\boldsymbol{A}^*\boldsymbol{A}$ ([***Normal Matrix***](https://en.wikipedia.org/wiki/Normal_matrix)), then:

$$
\rho \left( \boldsymbol{A} \right) =\left\| \boldsymbol{A} \right\| _2
$$
