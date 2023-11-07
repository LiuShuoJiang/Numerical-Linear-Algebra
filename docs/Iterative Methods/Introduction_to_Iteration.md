# General Strategy of Iteration Methods

We will discuss the general **Splitting Strategy**.

For linear system $\boldsymbol{Ax}=\boldsymbol{b}$, let $\boldsymbol{A}=\boldsymbol{M}-\boldsymbol{N}$. Then:

$$
\boldsymbol{Ax}=\boldsymbol{b}\,\,\Longleftrightarrow \left( \boldsymbol{M}-\boldsymbol{N} \right) \boldsymbol{x}=\boldsymbol{b}\,\,\Longleftrightarrow \boldsymbol{Mx}=\boldsymbol{b}+\boldsymbol{Nx}
$$

We can have the iteration:

$$
\boldsymbol{Mx}^{\left( k+1 \right)}=\boldsymbol{Nx}^{\left( k \right)}+\boldsymbol{b}
$$

If $\boldsymbol{M}$ is invertible, we have:

$$
\boldsymbol{x}^{\left( k+1 \right)}=\boldsymbol{M}^{-1}\boldsymbol{Nx}^{\left( k \right)}+\boldsymbol{M}^{-1}\boldsymbol{b}
$$

Denote $\boldsymbol{M}^{-1}\boldsymbol{b}=\boldsymbol{f}, \boldsymbol{M}^{-1}\boldsymbol{N}=\boldsymbol{G}$, we get:

$$
\boldsymbol{x}^{\left( k+1 \right)}=\boldsymbol{Gx}^{\left( k \right)}+\boldsymbol{f}
$$

where $\boldsymbol{G}$ is called **iterative matrix**.

For example, Jacobi method: $\boldsymbol{M}=\boldsymbol{D}, \boldsymbol{N}=\boldsymbol{E}+\boldsymbol{F}$.

`Question`: Is it convergent? If so, how fast?

If $\boldsymbol{x}^*$ is a solution of $\boldsymbol{Ax}=\boldsymbol{b}$, then:

$$
\boldsymbol{Ax}^*=\boldsymbol{b}\,\,\Longleftrightarrow \left( \boldsymbol{M}-\boldsymbol{N} \right) \boldsymbol{x}^*=\boldsymbol{b}\Longleftrightarrow \boldsymbol{Mx}^*=\boldsymbol{Nx}^*+\boldsymbol{b};
$$

$$
\boldsymbol{x}^*=\boldsymbol{M}^{-1}\boldsymbol{Nx}^*+\boldsymbol{M}^{-1}\boldsymbol{b}=\boldsymbol{Gx}^*+\boldsymbol{f};
$$

$$
\boldsymbol{x}^*=\boldsymbol{Gx}^*+\boldsymbol{f}
$$

Therefore, $\boldsymbol{x}^*$ is a fixed point of the iteration.

Note that $\boldsymbol{x}^{\left( k+1 \right)}=\boldsymbol{Gx}^{\left( k \right)}+\boldsymbol{f}$, then:

$$
\boldsymbol{x}^{\left( k+1 \right)}-\boldsymbol{x}^*=\boldsymbol{G}\left( \boldsymbol{x}^{\left( k \right)}-\boldsymbol{x}^* \right) 
$$

Define $\boldsymbol{e}^{\left( k \right)}=\boldsymbol{x}^{\left( k+1 \right)}-\boldsymbol{x}^*$. We get:

$$
\boldsymbol{e}^{\left( k+1 \right)}=\boldsymbol{Ge}^{\left( k \right)}
$$

This is the ***error equation***. 

`Question`: Is $\boldsymbol{e}^{\left( k \right)}\rightarrow 0$ as $k\rightarrow +\infty$ ? If so, how fast?

`Theorem`: If the spectral radius $\rho \left( \boldsymbol{G} \right) <1$, it implies that $(\mathbf{I}-\boldsymbol{G})$ is invertible and $\boldsymbol{x}^{\left( k \right)}\rightarrow \boldsymbol{x}^*$.

The inverse is also true: If $\boldsymbol{x}^{\left( k \right)}$ converges for any $\boldsymbol{f}$ and $\boldsymbol{x}^{(0)}$, then $\rho \left( \boldsymbol{G} \right) <1$. (How to prove?)

`Remarks`: For an iterative method, if it is convergent:

1. $\boldsymbol{b}$ can be arbitrary;
2. $\boldsymbol{x}^{(0)}$ can be arbitrary;
3. The method always converges to $\boldsymbol{x}^*$.

`Question`: If $\rho \left( \boldsymbol{G} \right) >1$, there exists an initial guess that $\boldsymbol{x}^{(k)}$ does not converge. If $\rho \left( \boldsymbol{G} \right) =1$, what happens?

`Theorem`: If $\left\| \boldsymbol{G} \right\| <1$ (any matrix norm is fine), then $\rho \left( \boldsymbol{G} \right) <1$.

`Definition`: A matrix $\boldsymbol{A}\in \mathbb{R} ^{m\times m}$ is called **(weakly) diagonally dominant** if:

$$
\left| a_{ii} \right|\geqslant \sum_{\begin{array}{c}
	j=1\\
	j\ne i\\
\end{array}}^m{\left| a_{ij} \right|}; i=1,2,\cdots ,m
$$

and is called **(strongly) diagonally dominant** if:

$$
\left| a_{ii} \right|>\sum_{\begin{array}{c}
	j=1\\
	j\ne i\\
\end{array}}^m{\left| a_{ij} \right|}; i=1,2,\cdots ,m
$$

`Theorem`: If $\boldsymbol{A}$ is a strongly diagonally dominant matrix, then the associate Jacobi, Gauss-Seidel iterations converge for any $\boldsymbol{x}^{(0)}$.

`Theorem`: If $\boldsymbol{A}$ is symmetric with positive diagonal elements, and $\omega \in \left( 0,2 \right)$, then SOR converges for any $\boldsymbol{x}^{(0)}$ if and only if $\boldsymbol{A}$ is positive definite.

`Theorem`([**Gershgorin Theroem**](https://en.wikipedia.org/wiki/Gershgorin_circle_theorem)): Given $\boldsymbol{A}\in \mathbb{R} ^{m\times m}$, the eigenvalues of $\boldsymbol{A}$ are contained in the union of the following disks:

$$
D_i=\left\{ z\in \mathbb{C} \middle| \left| z-a_i \right|\leqslant \sum_{\begin{array}{c}
	j=1\\
	j\ne i\\
\end{array}}^m{\left| a_{ij} \right|} \right\} ; i=1,2,\cdots ,m
$$

Let:

$$
r_i=\sum_{\begin{array}{c}
	j=1\\
	j\ne i\\
\end{array}}^m{\left| a_{ij} \right|}
$$

We have:

$$
D_i=\left\{ z\in \mathbb{C} \middle| \left| z-a_i \right|\leqslant r_i \right\} ; i=1,2,\cdots ,m
$$
