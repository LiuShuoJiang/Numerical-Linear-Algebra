# Other Versions of Gaussian Elimination

## Gaussian Elimination with Partial Pivoting

Partial pivoting only rearranges the rows of $\boldsymbol{A}$ and leaves the columns fixed.

Define the **permutation matrix**:

$$
\boldsymbol{P}\triangleq \left[ \begin{matrix}
	1&		&		&		&		&		&		&		&		&		&		\\
	&		\ddots&		&		&		&		&		&		&		&		&		\\
	&		&		1&		&		&		&		&		&		&		&		\\
	&		&		&		0_{\left( \mathrm{ith}\ \mathrm{row} \right)}&		&		\cdots&		&		1&		&		&		\\
	&		&		&		&		1&		&		&		&		&		&		\\
	&		&		&		\vdots&		&		\ddots&		&		\vdots&		&		&		\\
	&		&		&		&		&		&		1&		&		&		&		\\
	&		&		&		1&		&		\cdots&		&		0_{\left( \mathrm{jth}\ \mathrm{col} \right)}&		&		&		\\
	&		&		&		&		&		&		&		&		1&		&		\\
	&		&		&		&		&		&		&		&		&		\ddots&		\\
	&		&		&		&		&		&		&		&		&		&		1\\
\end{matrix} \right] 
$$

$\boldsymbol{PA}$ switches the $i$-th and $j$-th rows in $\boldsymbol{A}$.

We claim:

$$
\boldsymbol{L}_{m-1}\boldsymbol{P}_{m-1}\cdots \boldsymbol{L}_2\boldsymbol{P}_2\boldsymbol{L}_1\boldsymbol{P}_1\boldsymbol{A}=\boldsymbol{U}
$$

$$
\Longleftrightarrow \boldsymbol{PA}=\boldsymbol{LU}
$$

Define:

$$
\boldsymbol{L}_{m-1}\prime=\boldsymbol{L}_{m-1},
$$

$$
\boldsymbol{L}_{m-2}\prime=\boldsymbol{P}_{m-1}\boldsymbol{L}_{m-2}{\boldsymbol{P}_{m-1}}^{-1},
$$

$$
\boldsymbol{L}_{m-3}\prime=\boldsymbol{P}_{m-1}\boldsymbol{P}_{m-2}\boldsymbol{L}_{m-3}{\boldsymbol{P}_{m-2}}^{-1}{\boldsymbol{P}_{m-1}}^{-1},
$$

$$
\vdots 
$$

$$
\boldsymbol{L}_1\prime=\boldsymbol{P}_{m-1}\boldsymbol{P}_{m-2}\cdots \boldsymbol{P}_2\boldsymbol{L}_1{\boldsymbol{P}_2}^{-1}\cdots {\boldsymbol{P}_{m-2}}^{-1}{\boldsymbol{P}_{m-1}}^{-1}
$$

Then:

$$
\boldsymbol{L}_{m-1}\boldsymbol{P}_{m-1}\cdots \boldsymbol{L}_2\boldsymbol{P}_2\boldsymbol{L}_1\boldsymbol{P}_1=\boldsymbol{A},
$$

$$
\boldsymbol{L}_{m-1}\prime\underset{\boldsymbol{L}_{m-2}\prime}{\underbrace{\boldsymbol{P}_{m-1}\boldsymbol{L}_{m-2}{\boldsymbol{P}_{m-1}}^{-1}}}\boldsymbol{P}_{m-1}\boldsymbol{P}_{m-2}\boldsymbol{L}_{m-3}\cdots \boldsymbol{L}_2\boldsymbol{P}_2\boldsymbol{L}_1\boldsymbol{P}_1=\boldsymbol{A},
$$

$$
\underset{\boldsymbol{L}^{-1}}{\underbrace{\boldsymbol{L}_{m-1}\prime\boldsymbol{L}_{m-2}\prime\cdots \boldsymbol{L}_1}}\prime\underset{\boldsymbol{P}}{\underbrace{\boldsymbol{P}_{m-1}\cdots \boldsymbol{P}_1}}\boldsymbol{A}=\boldsymbol{U},
$$

$$
\Longrightarrow \boldsymbol{PA}=\boldsymbol{LU}
$$

`Questions`: How to write the algorithm? What is the cost?

Note: **Complete pivoting** (Full pivoting) rearranges both rows and columns. Although it is more stable, this method is a lot more expensive than partial pivoting regarding cost, and is not commonly used in practice.

## Cholesky Factorization

If $\boldsymbol{A}$ is symmetric positive definite (SPD), then $\boldsymbol{A}^T=\boldsymbol{A}$ and $\boldsymbol{A}$ has positive eigenvalues. We can also say that:

$$
\forall \boldsymbol{x}\in \mathbb{R} ^m,\boldsymbol{x}\ne \mathbf{0}:\ \boldsymbol{x}^T\boldsymbol{Ax}>0
$$

By Applying LU factorization to SPD matrix $\boldsymbol{A}$, we can find that $\boldsymbol{U}$'s diagonal parts are positive. Denote the diagonal part of $\boldsymbol{U}$ as $\boldsymbol{D}=\mathrm{diag}\left( \boldsymbol{U} \right)$. $\bar{\boldsymbol{U}}$ is diagonal 1 upper triangular matrix.

$$
\boldsymbol{A}=\boldsymbol{LU},\boldsymbol{U}=\boldsymbol{D}\bar{\boldsymbol{U}},
$$

$$
\boldsymbol{A}=\boldsymbol{LU}=\boldsymbol{LD}\bar{\boldsymbol{U}},
$$

$$
\boldsymbol{A}^T=\left( \boldsymbol{LD}\bar{\boldsymbol{U}} \right) ^T=\bar{\boldsymbol{U}}^T\boldsymbol{DL}^T=\boldsymbol{A}=\boldsymbol{LD}\bar{\boldsymbol{U}},
$$

$$
\Longrightarrow \bar{\boldsymbol{U}}^T=\boldsymbol{L},\boldsymbol{L}^T=\bar{\boldsymbol{U}},
$$

$$
\Longrightarrow \boldsymbol{A}=\boldsymbol{LDL}^T
$$

Define:

$$
\boldsymbol{D}^{\frac{1}{2}}=\left[ \begin{matrix}
	{d_{11}}^{\frac{1}{2}}&		&		\boldsymbol{O}\\
	&		\ddots&		\\
	\boldsymbol{O}&		&		{d_{mm}}^{\frac{1}{2}}\\
\end{matrix} \right] 
$$

We get:

$$
\boldsymbol{A}=\boldsymbol{LDL}^T=\boldsymbol{LD}^{\frac{1}{2}}\boldsymbol{D}^{\frac{1}{2}}\boldsymbol{L}^T=\left( \boldsymbol{LD}^{\frac{1}{2}} \right) \left( \boldsymbol{LD}^{\frac{1}{2}} \right) ^T
$$

Rename $\boldsymbol{LD}^{\frac{1}{2}}=\boldsymbol{R}$, then:

$$
\boldsymbol{A}=\boldsymbol{RR}^T
$$

This is **Cholesky factorization**.

Actually, there is no need for partial pivoting here!

$$
\boldsymbol{A}=\left[ \begin{matrix}
	a_{11}&		\boldsymbol{w}^T\\
	\boldsymbol{w}&		\boldsymbol{K}\\
\end{matrix} \right] ,a_{11}>0,\alpha =\sqrt{a_{11}};
$$

$$
\boldsymbol{A}=\underset{\boldsymbol{R}_1}{\underbrace{\left[ \begin{matrix}
	\alpha&		0\\
	\frac{\boldsymbol{w}}{\alpha}&		\mathbf{I}\\
\end{matrix} \right] }}\underset{\boldsymbol{A}_1}{\underbrace{\left[ \begin{matrix}
	1&		0\\
	0&		\boldsymbol{K}-\frac{\boldsymbol{ww}^T}{a_{11}}\\
\end{matrix} \right] }}\underset{{\boldsymbol{R}_1}^T}{\underbrace{\left[ \begin{matrix}
	\alpha&		\frac{\boldsymbol{w}^T}{\alpha}\\
	0&		\mathbf{I}\\
\end{matrix} \right] }}
$$

$$
=\boldsymbol{R}_1\boldsymbol{A}_1{\boldsymbol{R}_1}^T
$$

We can continue this step:

$$
\boldsymbol{A}=\boldsymbol{R}_1\boldsymbol{A}_1{\boldsymbol{R}_1}^T=\cdots =\underset{\boldsymbol{R}}{\underbrace{\boldsymbol{R}_1\cdots \boldsymbol{R}_m}}\mathbf{I}\underset{\boldsymbol{R}^T}{\underbrace{{\boldsymbol{R}_m}^T\cdots {\boldsymbol{R}_1}^T}}
$$

`Algorithm` ( **Cholesky Factorization** ):

Let initial $\boldsymbol{R}=\boldsymbol{A}$;

For $k=1,2,\cdots ,m$:

- For $j=k+1,\cdots ,m$:
    - $R_{j,j:m}=R_{j,j:m}-R_{j,j:m}\frac{R_{kj}}{R_{kk}}$;
- End;
- $R_{k,k:m}=\frac{R_{k,k:m}}{\sqrt{R_{kk}}}$;

End

The cost is $O\left( \frac{1}{3}m^3 \right) \approx O\left( m^3 \right)$.

