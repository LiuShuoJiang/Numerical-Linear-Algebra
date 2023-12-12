# Gaussian Elimination

**Lower–Upper (LU)** decomposition or factorization factors a matrix as the product of a lower triangular matrix and an upper triangular matrix. The product sometimes includes a permutation matrix as well. LU decomposition can be viewed as the matrix form of **Gaussian Elimination**.

## Introduction to LU Factorization

General steps of Gaussian Elimination:

$$
\left[ \begin{matrix}
	\times&		\times&		\times&		\times\\
	\times&		\times&		\times&		\times\\
	\times&		\times&		\times&		\times\\
	\times&		\times&		\times&		\times\\
\end{matrix} \right] \xrightarrow{\boldsymbol{L}_1}\left[ \begin{matrix}
	\times&		\times&		\times&		\times\\
	0&		\times&		\times&		\times\\
	0&		\times&		\times&		\times\\
	0&		\times&		\times&		\times\\
\end{matrix} \right] 
$$

$$
\xrightarrow{\boldsymbol{L}_2}\left[ \begin{matrix}
	\times&		\times&		\times&		\times\\
	0&		\times&		\times&		\times\\
	0&		0&		\times&		\times\\
	0&		0&		\times&		\times\\
\end{matrix} \right] \xrightarrow{\boldsymbol{L}_3}\left[ \begin{matrix}
	\times&		\times&		\times&		\times\\
	0&		\times&		\times&		\times\\
	0&		0&		\times&		\times\\
	0&		0&		0&		\times\\
\end{matrix} \right] ;
$$

$$
\Longrightarrow \boldsymbol{L}_3\boldsymbol{L}_2\boldsymbol{L}_1\boldsymbol{A}=\boldsymbol{U}
$$

`Question`: What is $\boldsymbol{L} _i$ ?

In essence:

$$
\boldsymbol{x}_k=\left[ \begin{array}{c}
	x_{1k}\\
	\vdots\\
	x_{kk}\\
	x_{k+1,k}\\
	\vdots\\
	x_{mk}\\
\end{array} \right] \xrightarrow{\boldsymbol{L}_k}\left[ \begin{array}{c}
	x_{1k}\\
	\vdots\\
	x_{kk}\\
	0\\
	\vdots\\
	0\\
\end{array} \right] 
$$

Define $l_{jk}\triangleq \frac{x_{jk}}{x_{kk}}$. Then we can do the $\left( j \right) -l_{jk}\cdot \left( k \right)$ row reduction. Also define:

$$
\boldsymbol{L}_k\triangleq \left[ \begin{matrix}
	1&		&		&		&		&		\\
	&		\ddots&		&		&		&		\\
	&		&		1&		&		&		\\
	&		&		-l_{k+1,k}&		\ddots&		&		\\
	&		&		\vdots&		&		\ddots&		\\
	&		&		-l_{m,k}&		&		&		1\\
\end{matrix} \right] 
$$

And:

$$
\boldsymbol{L}_k\cdot \left[ \begin{array}{c}
	x_{1k}\\
	\vdots\\
	x_{kk}\\
	x_{k+1,k}\\
	\vdots\\
	x_{mk}\\
\end{array} \right] \rightarrow \left[ \begin{array}{c}
	x_{1k}\\
	\vdots\\
	x_{kk}\\
	0\\
	\vdots\\
	0\\
\end{array} \right] ;
$$

$$
\boldsymbol{L}_k=\mathbf{I}-\boldsymbol{l}_k{\boldsymbol{e}_k}^T;\boldsymbol{l}_k=\left[ \begin{array}{c}
	0\\
	\vdots\\
	0\\
	l_{k+1,k}\\
	\vdots\\
	l_{mk}\\
\end{array} \right] ,\boldsymbol{e}_k=\left[ \begin{array}{c}
	0\\
	\vdots\\
	0\\
	1_{\left( k+1\ \mathrm{position} \right)}\\
	0\\
	\vdots\\
	0\\
\end{array} \right] 
$$

Then:

$$
\boldsymbol{L}_k=\mathbf{I}-\boldsymbol{l}_k{\boldsymbol{e}_k}^T,
$$

$$
\boldsymbol{L}_{m-1}\cdots \boldsymbol{L}_2\boldsymbol{L}_1\boldsymbol{A}=\boldsymbol{U}\Longrightarrow \boldsymbol{A}={\boldsymbol{L}_1}^{-1}{\boldsymbol{L}_2}^{-1}\cdots {\boldsymbol{L}_{m-1}}^{-1}\boldsymbol{U}
$$

Note that the inverse of a lower triangular matrix is also a lower triangular matrix. The product of lower triangular matrices is also a lower triangular matrix.

Define $\boldsymbol{L}\triangleq {\boldsymbol{L}_1}^{-1}{\boldsymbol{L}_2}^{-1}\cdots {\boldsymbol{L}_{m-1}}^{-1}$. Then:

$$
\boldsymbol{A}=\boldsymbol{LU}
$$

is the **LU factorization**.

Here are some useful conclusions:

$$
{\boldsymbol{L}_k}^{-1}=\mathbf{I}+\boldsymbol{l}_k{\boldsymbol{e}_k}^T
$$

$$
\Longleftarrow {\boldsymbol{L}_k}^{-1}\boldsymbol{L}_k=\left( \mathbf{I}+\boldsymbol{l}_k{\boldsymbol{e}_k}^T \right) \cdot \left( \mathbf{I}-\boldsymbol{l}_k{\boldsymbol{e}_k}^T \right) =\mathbf{I}-\boldsymbol{l}_k{\boldsymbol{e}_k}^T\boldsymbol{l}_k{\boldsymbol{e}_k}^T=\mathbf{I}
$$

Also note:

$$
{\boldsymbol{L}_k}^{-1}{\boldsymbol{L}_{k+1}}^{-1}=\left( \mathbf{I}+\boldsymbol{l}_k{\boldsymbol{e}_k}^T \right) \cdot \left( \mathbf{I}+\boldsymbol{l}_{k+1}{\boldsymbol{e}_{k+1}}^T \right) 
$$

$$
=\mathbf{I}+\boldsymbol{l}_k{\boldsymbol{e}_k}^T+\boldsymbol{l}_{k+1}{\boldsymbol{e}_{k+1}}^T+\boldsymbol{l}_k{\boldsymbol{e}_k}^T\boldsymbol{l}_{k+1}{\boldsymbol{e}_{k+1}}^T
$$

$$
=\mathbf{I}+\boldsymbol{l}_k{\boldsymbol{e}_k}^T+\boldsymbol{l}_{k+1}{\boldsymbol{e}_{k+1}}^T
$$

$$
=\left[ \begin{matrix}
	1&		&		&		&		&		&		\\
	&		\ddots&		&		&		&		&		\\
	&		&		1&		&		&		&		\\
	&		&		l_{k+1,k}&		1&		&		&		\\
	&		&		\vdots&		l_{k+2,k+1}&		\ddots&		&		\\
	&		&		\vdots&		\vdots&		&		\ddots&		\\
	&		&		l_{m,k}&		l_{m,k+1}&		&		&		1\\
\end{matrix} \right] 
$$

This is good for bookkeeping the solutions. Therefore:

$$
\boldsymbol{L}={\boldsymbol{L}_1}^{-1}{\boldsymbol{L}_2}^{-1}\cdots {\boldsymbol{L}_{m-1}}^{-1}=\left[ \begin{matrix}
	1&		&		&		&		&		\\
	l_{21}&		\ddots&		&		&		&		\\
	\vdots&		&		1&		&		&		\\
	\vdots&		&		l_{k+1,k}&		\ddots&		&		\\
	\vdots&		&		\vdots&		&		\ddots&		\\
	l_{m1}&		&		l_{m,k}&		&		&		1\\
\end{matrix} \right] ;
$$

$$
\boldsymbol{A}=\boldsymbol{LU}
$$

To solve $\boldsymbol{Ax}=\boldsymbol{b}$, we may solve $\boldsymbol{LUx}=\boldsymbol{b}$ instead. Define $\boldsymbol{y}=\boldsymbol{Ux}$. We can solve two linear systems:

$$
\boldsymbol{Ly}=\boldsymbol{b}, \boldsymbol{Ux}=\boldsymbol{y}
$$

## Steps of Gaussian Elimination

`Algorithm` ( **Simple Gaussian Elimination (GE)** ):

( Given $\boldsymbol{A}$, we output $\boldsymbol{L}, \boldsymbol{U}$. )

Let $\boldsymbol{U}=\boldsymbol{A},\boldsymbol{L}=\mathbf{I}$ first;

For $k=1,2,\cdots ,m-1$:

- For $j=k+1,k+2,\cdots ,m$:
    - $l_{j,k}=\frac{u_{j,k}}{u_{k,k}}$;
    - $u_{j,k:m}=u_{j,k:m}-l_{j,k}\cdot u_{k,k:m}$;
- End;

End

There are actually *three* loops for this algorithm.

`Question`: What is the cost of the algorithm?

We measure the computational complexity by the number of flops (floating point operations):

$$
\sum_{k=1}^{m-1}{\sum_{j=k+1}^m{2\left( m-k+1 \right) +1}}\approx \frac{2}{3}m^3+O\left( m^2 \right) 
$$

Therefore, the cost of LU factorization is $O\left( m^3 \right)$.

## Forward and Backward Substitution

Consider $\boldsymbol{Ly}=\boldsymbol{b}$, we do **Forward Substitution**:

For $j=1,2,\cdots ,m$:

$$
y_j=b_j-\sum_{k=1}^{j-1}{l_{kj}y_k}
$$

End. The cost for it is $O\left( m^2 \right)$.

Consider $\boldsymbol{Ux}=\boldsymbol{y}$, we do **Backward Substitution**:

For $j=m,m-1,\cdots ,1$:

$$
x_j=\frac{1}{r_{jj}}\left( y_j-\sum_{k=j+1}^m{x_kr_{jk}} \right) 
$$

## Evaluation

***Gaussian-Elimination (Native) is not stable***!

e.g., Matrix $\boldsymbol{A}=\left[ \begin{matrix}
	0&		1\\
	1&		1\\
\end{matrix} \right]$ fails at the first step.

e.g., For matrix

$$
\boldsymbol{A}=\left[ \begin{matrix}
	10^{-20}&		1\\
	1&		1\\
\end{matrix} \right] 
$$

Apply GE, we get:

$$
\boldsymbol{L}=\left[ \begin{matrix}
	1&		0\\
	10^{20}&		1\\
\end{matrix} \right] ,\boldsymbol{U}=\left[ \begin{matrix}
	10^{-20}&		1\\
	0&		1-10^{20}\\
\end{matrix} \right] 
$$

However, in computer representation, the actual matrix might be like:

$$
\tilde{\boldsymbol{L}}=\left[ \begin{matrix}
	1&		0\\
	10^{20}&		1\\
\end{matrix} \right] ,\tilde{\boldsymbol{U}}=\left[ \begin{matrix}
	10^{-20}&		1\\
	0&		-10^{20}\\
\end{matrix} \right] 
$$

Then:

$$
\tilde{\boldsymbol{L}}\tilde{\boldsymbol{U}}=\left[ \begin{matrix}
	10^{-20}&		1\\
	1&		0\\
\end{matrix} \right] \ne \boldsymbol{A}
$$

A small mistake can lead to a large perturbation.

Solution to this problem is partial pivoting.
