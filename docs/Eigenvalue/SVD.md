# Singular Value Decomposition (SVD)

Unitary Matrix

[**Unitary matrix**](https://en.wikipedia.org/wiki/Unitary_matrix): $\boldsymbol{Q}$ is a unitary matrix if:

$$
\boldsymbol{QQ}^*=\mathbf{I}, \boldsymbol{Q}^*\boldsymbol{Q}=\mathbf{I}
$$

where $\boldsymbol{Q}^*$ is the *conjugate transpose* of $\boldsymbol{Q}$.

A unitary matrix satisfies $\left\| \boldsymbol{QA} \right\| _2=\left\| \boldsymbol{A} \right\| _2$, which means that a unitary transform does not change the 2-norm of a matrix. Similarly, it can also be verified that $\left\| \boldsymbol{QA} \right\| _F=\left\| \boldsymbol{A} \right\| _F$, which means that a unitary transform does not change the Frobenius norm of a matrix.

Introduction to SVD

Example: Let $\boldsymbol{A}\in \mathbb{R} ^{2\times 2}, \boldsymbol{A}:\mathbb{R} ^2\mapsto \mathbb{R} ^2$. What is the image of the unit disk of $\boldsymbol{A}$ ?

We get:

$$
\boldsymbol{Av}_1=\sigma _1\boldsymbol{u}_1;\boldsymbol{Av}_2=\sigma _2\boldsymbol{u}_2
$$

$$
\boldsymbol{A}\left[ \begin{matrix}
	\boldsymbol{v}_1&		\boldsymbol{v}_2\\
\end{matrix} \right] =\left[ \begin{matrix}
	\sigma _1\boldsymbol{u}_1&		\sigma _2\boldsymbol{u}_2\\
\end{matrix} \right] =\left[ \begin{matrix}
	\boldsymbol{u}_1&		\boldsymbol{u}_2\\
\end{matrix} \right] \cdot \left[ \begin{matrix}
	\sigma _1&		0\\
	0&		\sigma _2\\
\end{matrix} \right] 
$$

Let:

$$
\boldsymbol{V}=\left[ \begin{matrix}
	\boldsymbol{v}_1&		\boldsymbol{v}_2\\
\end{matrix} \right] , \boldsymbol{U}=\left[ \begin{matrix}
	\boldsymbol{u}_1&		\boldsymbol{u}_2\\
\end{matrix} \right] , \mathbf{\Sigma }=\left[ \begin{matrix}
	\sigma _1&		0\\
	0&		\sigma _2\\
\end{matrix} \right] 
$$

We can get:

$$
\boldsymbol{AV}=\boldsymbol{U}\mathbf{\Sigma }
$$

where

$$
\boldsymbol{VV}^*=\mathbf{I}, \boldsymbol{UU}^*=\mathbf{I}
$$

Then $\boldsymbol{A}=\boldsymbol{AVV}^*=\boldsymbol{U}\mathbf{\Sigma }\boldsymbol{V}^*$.

The **Singular Value Decomposition (SVD)** is:

$$
\boldsymbol{A}=\boldsymbol{U}\mathbf{\Sigma }\boldsymbol{V}^*
$$

where $\boldsymbol{U}\in \mathbb{C} ^{m\times m}$ is a complex unitary matrix, $\mathbf{\Sigma }\in \mathbb{R} ^{m\times n}$ is a rectangular diagonal matrix with non-negative real numbers on the diagonal, and $\boldsymbol{V}\in \mathbb{C} ^{n\times n}$ is a complex unitary matrix. ($\boldsymbol{A}\in \mathbb{C} ^{m\times n}$ is a complex matrix)

The diagonal entries $\sigma _i=\Sigma _{ii}$ of $\mathbf{\Sigma }$ are uniquely determined by $\boldsymbol{A}$ and are called **singular values** of $\boldsymbol{A}$.

`Theorem`(Existence and Uniqueness): Every matrix $\boldsymbol{A}\in \mathbb{C} ^{m\times n}$ has a SVD. The singular values are uniquely determined. If $\boldsymbol{A}$ is square and $\sigma _j$ are distinct, the left and right singular vectors $\boldsymbol{u}_i,\boldsymbol{v}_i$ are uniquely determined up to a complex sign (complex number with modular one).

Properties of SVD:

From:

$$
\boldsymbol{A}^*\boldsymbol{A}=\left( \boldsymbol{U}\mathbf{\Sigma }\boldsymbol{V}^* \right) ^*\cdot \boldsymbol{U}\mathbf{\Sigma }\boldsymbol{V}^*
$$

$$
=\left( \boldsymbol{V}\mathbf{\Sigma }^*\boldsymbol{U}^* \right) \boldsymbol{U}\mathbf{\Sigma }\boldsymbol{V}^*=\boldsymbol{V}\mathbf{\Sigma }^2\boldsymbol{V}^*=\boldsymbol{V}\mathbf{\Sigma }^2\boldsymbol{V}^{-1}
$$

We know that ${\sigma _i}^2$ is the eigenvalue of $\boldsymbol{A}^*\boldsymbol{A}$.

We can also get $\mathrm{rank}\left( \boldsymbol{A} \right)$ is the number of non-zero elements of singular values.

Furthermore, $\left\| \boldsymbol{A} \right\| _2=\sigma _1$ (largest singular value) if $\sigma _1\geqslant \sigma _2\geqslant \cdots \geqslant \sigma _r$.

Also:

$$
\left\| \boldsymbol{A} \right\| _F=\left( \sum_{i=1}^r{{\sigma _i}^2} \right) ^{\frac{1}{2}}
$$

If $\boldsymbol{A}\in \mathbb{R} ^{m\times m}$, then:

$$
\left| \det \left( \boldsymbol{A} \right) \right|=\prod_{i=1}^m{\sigma _i}
$$


