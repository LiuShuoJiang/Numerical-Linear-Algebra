# Introduction to Eigenvalue Problem

`Overall Task`: Given

$$
\boldsymbol{Ax}=\lambda \boldsymbol{x}
$$

We need to find $\lambda$ and $\boldsymbol{x}$.

What are the methods to compute eigenvalues and eigenvectors?

## Traditional Methods

### Traditional Methods Introduction

Traditional Method:

- *Step 1*: Calculate the character polynomial $p_{\boldsymbol{A}}\left( z \right) =\det \left( z\boldsymbol{I}-\boldsymbol{A} \right)$;
- *Step 2*: Find the roots of $p_{\boldsymbol{A}}\left( z \right)$: $\lambda _i$. They are eigenvalues;
- *Step 3*: Solve $\boldsymbol{Ax}=\lambda _i\boldsymbol{x}$ for $\boldsymbol{x}$. $\boldsymbol{x}$ is the eigenvalue corresponding to $\lambda _i$.

### Problem of Traditional Methods

`Example`: Let

$$
\boldsymbol{A}=\left[ \begin{matrix}
	1&		1\\
	0&		1\\
\end{matrix} \right] 
$$

Then $\lambda _{1,2}=1$, $p_{\boldsymbol{A}}\left( z \right) =\left( z-1 \right) ^2=z^2-2z+1$. The coefficients encodes in the computer are $\left[ \begin{matrix}
	1&		-2&		1\\
\end{matrix} \right]$. Assume that the input is $\left[ \begin{matrix}
	1&		-2&		1-10^{-16}\\
\end{matrix} \right]$. Then the polynomial becomes $z^2-2z+\left( 1-10^{-16} \right) =\left( z-\left( 1-10^{-8} \right) \right) \cdot \left( z-\left( 1+10^{-8} \right) \right) =0$. We get $\lambda _1=1-10^{-8}, \lambda _2=1+10^{-8}$. Therefore, when the input error is $10^{-16}$, the output error is $O(10^{-8})$. The root finding procedure amplifies the error!

Root finding is not stable, and algorithms using root finding procedures are not a stable algorithms.

### Root Finding V.S. Eigenvalue Problems

#### Comparison

Actually, root finding in computers is implemented as eigenvalue solving algorithms.

`Claim`: **Eigenvalue problem is equivalent to root finding for polynomials**.

Given $p\left( z \right) =z^m+a_{m-1}z^{m-1}+\cdots +a_1z+a_0$. Introduce a matrix:

$$
\boldsymbol{A}=\left[ \begin{matrix}
	0&		&		&		O&		-a_0\\
	1&		\ddots&		&		&		-a_1\\
	&		\ddots&		\ddots&		&		\\
	&		&		\ddots&		0&		-a_{m-2}\\
	O&		&		&		1&		-a_{m-1}\\
\end{matrix} \right] 
$$

$\boldsymbol{A}\in \mathbb{R} ^{m\times m}$ is called the **companion matrix** for $p(z)$. We can prove that:

$$
p_{\boldsymbol{A}}\left( z \right) =\det \left( z\boldsymbol{I}-\boldsymbol{A} \right) =p\left( z \right) 
$$

Then getting the roots for the original polynomial is euqivalent to getting the eigenvalues for the corresponding companion matrix.

#### Difficulty

There is **fundamental difficulty** in the eigenvalue computation: no explicit formula for a *general* matrix $\boldsymbol{A}_{m\times m}$ when $m\geqslant 5$.

`Theorem`: For any $m\geqslant 5$, there is a polynomial $p(z)$ of degree $m$ with rational coefficients that has a real root $p(r)=0$ with the property that $r$ can NOT be written using any expression involving rational number additions, subtractions, multiplications, divisions or $k$ -th roots.

> **All the methods for eigenvalues must be iterative!**

Even if working with exact arithmetic, there could be no computing program that would produce the exact roots of an arbitrary polynomial of degree $m\geqslant 5$ in a finite number of steps. All eigenvalue methods must be iterative.

## Numerical Methods

We want to convert the original matrix to a matrix that is easy to recognize the eigenvalues (upper/lower triangular or diagonal), and the transformation should not change the eigenvalues of matrices (this is the concept of **similar**).

Commonly used numerical methods for eigenvalue problem have ***two phases***:

### Phase One

**Phase 1**: A **direct method** to produce an upper **Hessenberg matrix** $\boldsymbol{H}$ (an *upper Hessenberg matrix* has zero entries below the first subdiagonal, and a lower Hessenberg matrix has zero entries above the first superdiagonal):

$$
\boldsymbol{H}=\left[ \begin{matrix}
	a_{11}&		a_{12}&		a_{13}&		a_{14}&		\cdots&		\cdots&		a_{1\left( n-1 \right)}&		a_{1n}\\
	a_{21}&		a_{22}&		a_{23}&		a_{24}&		\ddots&		\cdots&		a_{2\left( n-1 \right)}&		a_{2n}\\
	0&		a_{32}&		a_{33}&		a_{34}&		\ddots&		\cdots&		a_{3\left( n-1 \right)}&		a_{3n}\\
	0&		0&		a_{43}&		a_{44}&		\ddots&		\cdots&		a_{4\left( n-1 \right)}&		a_{4n}\\
	0&		0&		0&		a_{54}&		\ddots&		\ddots&		a_{5\left( n-1 \right)}&		a_{5n}\\
	\vdots&		\vdots&		\vdots&		0&		\ddots&		\vdots&		\vdots&		\vdots\\
	0&		0&		0&		0&		\ddots&		a_{\left( n-1 \right) \left( n-2 \right)}&		a_{\left( n-1 \right) \left( n-1 \right)}&		a_{\left( n-1 \right) n}\\
	0&		0&		0&		0&		\cdots&		0&		a_{n\left( n-1 \right)}&		a_{nn}\\
\end{matrix} \right] 
$$

The procedure can end within finite number of steps. This is done by **similar transforms**. In linear algebra, $\boldsymbol{A}$ and $\boldsymbol{B}$ are similar ( $\boldsymbol{A}\sim \boldsymbol{B}$ ), if there exists $\boldsymbol{X}\in \mathbb{R} ^{m\times m}$ invertible such that $\boldsymbol{XAX}^{-1}=\boldsymbol{B}$. We can let:

$$
0=\det \left( \lambda \boldsymbol{I}-\boldsymbol{B} \right) =\det \left( \lambda \boldsymbol{I}-\boldsymbol{XAX}^{-1} \right) =\det \left( \boldsymbol{X}\left( \lambda \boldsymbol{I}-\boldsymbol{A} \right) \boldsymbol{X}^{-1} \right) =\det \left( \boldsymbol{X} \right) \cdot \det \left( \lambda \boldsymbol{I}-\boldsymbol{A} \right) \cdot \det \left( \boldsymbol{X}^{-1} \right)
$$

Because $\det \left( \boldsymbol{X} \right) \ne 0, \det \left( \boldsymbol{X}^{-1} \right) \ne 0$, we know that the characteristic polynomials of $\boldsymbol{A}$ and $\boldsymbol{B}$ are the same, and $\boldsymbol{A}$ and $\boldsymbol{B}$ have the same set of eigenvalues.

### Phase Two

**Phase 2**: An **iterative procedure** that is also based on *similarity transforms* to produce a formally infinite sequence of upper Hessenberg matrices that converges to a triangular form.

Why we divide the algorithms into two phases? The cost concern is the reason for the two-phase strategy. Phase 1 has $O(m^3)$ flops, while in Phase 2, it normally converges to $O\left( \varepsilon _{machine} \right)$ within $O(m)$ steps, each step requiring at most $O(m^2)$ steps to finish.

Building blocks: [**QR Factorization**](./QR_Factorization.md).
