# Eigenvalue Problems

[Introduction to Eigenvalue Problems](./Introduction.md)

[QR Factorization](./QR_Factorization.md)

[Classical Eigenvalue Algorithms](./Classical_Eigenvalue_Algorithms.md)

[QR Algorithm](./QR_Algorithm.md)

[Divide and Conquer](./Divide_and_Conquer.md)

[SVD](./SVD.md)

`Quick review`: **Gram-Schmidt Process**

If a vector group $\boldsymbol{\alpha }_1,\boldsymbol{\alpha }_2,\cdots ,\boldsymbol{\alpha }_k$ is linearly independent, let:

$$
\boldsymbol{\beta }_1=\boldsymbol{\alpha }_1,
$$

$$
\boldsymbol{\beta }_2=\boldsymbol{\alpha }_2-\frac{\left( \boldsymbol{\alpha }_2,\boldsymbol{\beta }_1 \right)}{\left( \boldsymbol{\beta }_1,\boldsymbol{\beta }_1 \right)}\boldsymbol{\beta }_1,
$$

$$
\boldsymbol{\beta }_3=\boldsymbol{\alpha }_3-\frac{\left( \boldsymbol{\alpha }_3,\boldsymbol{\beta }_1 \right)}{\left( \boldsymbol{\beta }_1,\boldsymbol{\beta }_1 \right)}\boldsymbol{\beta }_1-\frac{\left( \boldsymbol{\alpha }_2,\boldsymbol{\beta }_2 \right)}{\left( \boldsymbol{\beta }_2,\boldsymbol{\beta }_2 \right)}\boldsymbol{\beta }_2,
$$

$$
\vdots 
$$

$$
\boldsymbol{\beta }_k=\boldsymbol{\alpha }_k-\frac{\left( \boldsymbol{\alpha }_k,\boldsymbol{\beta }_1 \right)}{\left( \boldsymbol{\beta }_1,\boldsymbol{\beta }_1 \right)}\boldsymbol{\beta }_1-\frac{\left( \boldsymbol{\alpha }_k,\boldsymbol{\beta }_2 \right)}{\left( \boldsymbol{\beta }_2,\boldsymbol{\beta }_2 \right)}\boldsymbol{\beta }_2-\cdots -\frac{\left( \boldsymbol{\alpha }_k-\boldsymbol{\beta }_{k-1} \right)}{\left( \boldsymbol{\beta }_{k-1},\boldsymbol{\beta }_{k-1} \right)}\boldsymbol{\beta }_{k-1}
$$

Then $\boldsymbol{\beta }_1,\boldsymbol{\beta }_2,\cdots ,\boldsymbol{\beta }_k$ are orthogonal to each other. Unitize them, we get:

$$
\boldsymbol{\gamma }_1=\frac{\boldsymbol{\beta }_1}{\left\| \boldsymbol{\beta }_1 \right\|},\boldsymbol{\gamma }_2=\frac{\boldsymbol{\beta }_2}{\left\| \boldsymbol{\beta }_2 \right\|},\cdots ,\boldsymbol{\gamma }_k=\frac{\boldsymbol{\beta }_k}{\left\| \boldsymbol{\beta }_k \right\|}
$$

The process from $\boldsymbol{\alpha }_1,\boldsymbol{\alpha }_2,\cdots ,\boldsymbol{\alpha }_k$ to $\boldsymbol{\gamma }_1,\boldsymbol{\gamma }_2,\cdots ,\boldsymbol{\gamma }_k$ is called **Gram-Schmidt Orthogonalization**. We can also get:

$$
\boldsymbol{A}=\left[ \begin{matrix}
	\boldsymbol{\alpha }_1&		\boldsymbol{\alpha }_2&		\cdots&		\boldsymbol{\alpha }_k\\
\end{matrix} \right] 
$$

$$
=\left[ \begin{matrix}
	\boldsymbol{\gamma }_1&		\boldsymbol{\gamma }_2&		\cdots&		\boldsymbol{\gamma }_k\\
\end{matrix} \right]\cdot \left[ \begin{matrix}
	\left\| \boldsymbol{\beta }_1 \right\|&		\left( \boldsymbol{\alpha }_2,\boldsymbol{\gamma }_1 \right)&		\cdots&		\left( \boldsymbol{\alpha }_k,\boldsymbol{\gamma }_1 \right)\\
	0&		\left\| \boldsymbol{\beta }_2 \right\|&		\cdots&		\left( \boldsymbol{\alpha }_k,\boldsymbol{\gamma }_2 \right)\\
	\vdots&		\ddots&		\ddots&		\vdots\\
	0&		\cdots&		0&		\left\| \boldsymbol{\beta }_k \right\|\\
\end{matrix} \right] 
$$
