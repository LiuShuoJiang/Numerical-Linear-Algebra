# QR Algorithm

**QR Algorithm** is the recommended algorithm for eigenvalue and eigenvector computation.

The matrix $\boldsymbol{A}$ below has been computed after phase 1.

## "Pure" QR Algorithm

**Pure QR Algorithm**:

Initialize: $\boldsymbol{A}^{\left( 0 \right)}=\boldsymbol{A}$;

For $k=1,2,\cdots$:

- Do QR factorization: $\boldsymbol{Q}^{\left( k \right)}\cdot \boldsymbol{R}^{\left( k \right)}=\boldsymbol{A}^{\left( k-1 \right)}$;
- $\boldsymbol{A}^{\left( k \right)}=\boldsymbol{R}^{\left( k \right)}\cdot \boldsymbol{Q}^{\left( k \right)}$;

End

### Similarity Transformation

**Task 1**: Verify that **QR algorithm performs similarity transformation**, that is, $\boldsymbol{A}^{\left( k \right)}\sim \boldsymbol{A}^{\left( k-1 \right)}$.

- *Step 1*: $\boldsymbol{A}^{\left( k-1 \right)}=\boldsymbol{Q}^{\left( k \right)}\cdot \boldsymbol{R}^{\left( k \right)}\Leftrightarrow \left( \boldsymbol{Q}^{\left( k \right)} \right) ^T\boldsymbol{A}^{\left( k-1 \right)}=\boldsymbol{R}^{\left( k \right)}$
- *Step 2*: $\boldsymbol{A}^{\left( k \right)}=\boldsymbol{R}^{\left( k \right)}\cdot \boldsymbol{Q}^{\left( k \right)}=\left( \boldsymbol{Q}^{\left( k \right)} \right) ^T\boldsymbol{A}^{\left( k-1 \right)}\boldsymbol{Q}^{\left( k \right)}\Rightarrow \boldsymbol{A}^{\left( k \right)}\sim \boldsymbol{A}^{\left( k-1 \right)}$. Then we can get:

$$
\boldsymbol{A}^{\left( k \right)}=\left( \boldsymbol{Q}^{\left( k \right)} \right) ^T\boldsymbol{A}^{\left( k-1 \right)}\boldsymbol{Q}^{\left( k \right)}=\left( \boldsymbol{Q}^{\left( k \right)} \right) ^T\left( \boldsymbol{Q}^{\left( k-1 \right)} \right) ^T\boldsymbol{A}^{\left( k-2 \right)}\boldsymbol{Q}^{\left( k \right)}=\cdots 
$$

$$
=\underset{\left( \bar{\boldsymbol{Q}}^{\left( k \right)} \right) ^T}{\underbrace{\left( \boldsymbol{Q}^{\left( k \right)} \right) ^T\left( \boldsymbol{Q}^{\left( k-1 \right)} \right) ^T\cdots \left( \boldsymbol{Q}^{\left( 1 \right)} \right) ^T}}\cdot \boldsymbol{A}^{\left( 0 \right)}\cdot \underset{\bar{\boldsymbol{Q}}^{\left( k \right)}}{\underbrace{\boldsymbol{Q}^{\left( 1 \right)}\cdots \boldsymbol{Q}^{\left( k-1 \right)}\boldsymbol{Q}^{\left( k \right)}}}
\\
=\left( \bar{\boldsymbol{Q}}^{\left( k \right)} \right) ^T\boldsymbol{A}\bar{\boldsymbol{Q}}^{\left( k \right)}
$$

If $\boldsymbol{A}^{\left( k \right)}$ converges to $\mathrm{diag}\left( \lambda _1,\cdots ,\lambda _n \right)$ as $k\rightarrow +\infty$, then $\bar{\boldsymbol{Q}}^{\left( k \right)}\rightarrow \boldsymbol{Q}$ which contains the eigenvectors of $\boldsymbol{A}$.

### About Phase One

**Task 2**: Why do we need **phase 1** *before* the algorithm performs?

Let us assume $\boldsymbol{A}$ is symmetric for simplicity, then $\boldsymbol{H}=\boldsymbol{QAQ}^T$ is tridiagonal. Starting from $\boldsymbol{H}$, the cost for each iteration is of $O(m)$ (using Given Rotation, for example).

If we don't perform phase 1, starting from $\boldsymbol{A}$, the cost for each iteration is $O(m^3)$ (using Householder QR factorization, for example).

Question: What is the cost per iteration if $\boldsymbol{H}$ is and upper Hessenberg matrix? (HW Practice)

### Convergence

**Task 3**: **QR algorithm is convergent**.

> `Theorem`: Let the pure QR algorithm be applied to a real *symmetric* matrix $\boldsymbol{A}$ whose eigenvalues satisfy $\left| \lambda _1 \right|>\left| \lambda _2 \right|>\cdots >\left| \lambda _m \right|$, and the corresponding eigenmatrix $\boldsymbol{Q}$ has all *nonsingular* **leading principal submatrices**. Then as $k\rightarrow +\infty$, $\boldsymbol{A}^{\left( k \right)}$ converges **linearly** with a constant rate given by $\max_j \frac{\left| \lambda _{j+1} \right|}{\left| \lambda _j \right|}$ to diagonal matrix $\mathrm{diag}\left( \lambda _1,\lambda _2,\cdots ,\lambda _n \right)$, and $\bar{\boldsymbol{Q}}^{\left( k \right)}$ (with signs of its columns adjusted as necessary) converges to $\boldsymbol{Q}$ at the same rate.

In the theorem above, $\boldsymbol{A}=\boldsymbol{Q}^T\mathbf{\Lambda }\boldsymbol{Q}$ where $\boldsymbol{Q}\in \mathbb{C} ^{m\times m}$ is unitary.

The sketch of the proof (in three steps):

#### Step One

***Step 1***: QR algorithm is essentially equivalent to the so-called simultaneous iteration.

**Simultaneous Iteration**:

Pick $\hat{\boldsymbol{Q}}^{\left( 0 \right)}\in \mathbb{R} ^{m\times n}$ with orthonormal columns;

For $k=1,2,\cdots$:

- $\boldsymbol{Z}^{\left( k \right)}=\boldsymbol{A}\hat{\boldsymbol{Q}}^{\left( k-1 \right)}$;
- Do QR factorization: $\hat{\boldsymbol{Q}}^{\left( k \right)}\hat{\boldsymbol{R}}^{\left( k \right)}=\boldsymbol{Z}^{\left( k \right)}$;

End

Claim: If we pick $\hat{\boldsymbol{Q}}^{\left( 0 \right)}=\mathbf{I}$, then the simultaneous iteration becomes the QR algorithm. (How to verify?)

#### Step Two

***Step 2***: Simultaneous iteration is a block Power Iteration.

- Power Iteration: $\boldsymbol{Z}^{\left( k \right)}=\boldsymbol{A}\hat{\boldsymbol{Q}}^{\left( k-1 \right)}$.
- Normalization and estimating eigenvalues: $\hat{\boldsymbol{Q}}^{\left( k \right)}\hat{\boldsymbol{R}}^{\left( k \right)}=\boldsymbol{Z}^{\left( k \right)}$.

#### Step Three

***Step 3***: The block Power Iteration is convergent with the rate given by

$$
c=\max_j \frac{\left| \lambda _{j+1} \right|}{\left| \lambda _j \right|}
$$

if $\left| \lambda _1 \right|>\left| \lambda _2 \right|>\cdots >\left| \lambda _j \right|>\left| \lambda _{j+1} \right|>\cdots >\left| \lambda _m \right|\geqslant 0$. Then we get $0<c<1$.

A special case for block Power Iteration to show the convergence (this is the essential idea):

Consider $\left[ \begin{matrix}
	{\boldsymbol{v}_1}^{\left( 0 \right)}&		{\boldsymbol{v}_{\left( 2 \right)}}^{\left( 0 \right)}\\
\end{matrix} \right]$ where ${\boldsymbol{v}_1}^{\left( 0 \right)}$ and ${\boldsymbol{v}_{\left( 2 \right)}}^{\left( 0 \right)}$ are unit vectors and orthogonal to each other. Then:

$$
\boldsymbol{A}\left[ \begin{matrix}
	{\boldsymbol{v}_1}^{\left( 0 \right)}&		{\boldsymbol{v}_{\left( 2 \right)}}^{\left( 0 \right)}\\
\end{matrix} \right] =\left[ \begin{matrix}
	{\boldsymbol{Av}_1}^{\left( 0 \right)}&		{\boldsymbol{Av}_{\left( 2 \right)}}^{\left( 0 \right)}\\
\end{matrix} \right] =\left[ \begin{matrix}
	\boldsymbol{w}_1&		\boldsymbol{w}_2\\
\end{matrix} \right] 
$$

Also:

$$
\boldsymbol{QR}=\left[ \begin{matrix}
	\boldsymbol{w}_1&		\boldsymbol{w}_2\\
\end{matrix} \right] ;
\\
\boldsymbol{Q}=\left[ \begin{matrix}
	{\boldsymbol{v}_1}^{\left( 1 \right)}&		{\boldsymbol{v}_{\left( 2 \right)}}^{\left( 1 \right)}\\
\end{matrix} \right] ;\cdots 
$$

$\boldsymbol{Q}$ spans the space generated by the eigenvectors $\boldsymbol{q}_1, \boldsymbol{q}_2$. We can further prove that ${\boldsymbol{v}_1}^{\left( k \right)}\rightarrow \boldsymbol{q}_1, {\boldsymbol{v}_2}^{\left( k \right)}\rightarrow \boldsymbol{q}_2$ as $k\rightarrow +\infty$.

## "Practical" QR Algorithm

### Algorithm Steps

We assume $\boldsymbol{A}$ is SPD below.

Phase 1 is the same as before: $\left( \boldsymbol{Q}^{\left( 0 \right)} \right) ^T\boldsymbol{A}^{\left( 0 \right)}\boldsymbol{Q}^{\left( 0 \right)}=\boldsymbol{A}$ where $\boldsymbol{A}^{\left( 0 \right)}$ is a tridiagonal matrix;

For $k=1,2,\cdots$:

- Pick a shift $\mu ^{\left( k \right)}$;
- Do QR factorization: $\boldsymbol{Q}^{\left( k \right)}\boldsymbol{R}^{\left( k \right)}=\boldsymbol{A}^{\left( k \right)}-\mu ^{\left( k \right)}\mathbf{I}$;
- $\boldsymbol{A}^{\left( k \right)}=\boldsymbol{R}^{\left( k \right)}\boldsymbol{Q}^{\left( k \right)}+\mu ^{\left( k \right)}\mathbf{I}$;
- If any off-diagonal element ${\boldsymbol{A}_{j,j+1}}^{\left( k \right)}$ is sufficiently small (close to zero is absolute value sense), we set ${\boldsymbol{A}_{j,j+1}}^{\left( k \right)}={\boldsymbol{A}_{j+1,j}}^{\left( k \right)}=0$ and obtain: $\boldsymbol{A}^{\left( k \right)}=\left[ \begin{matrix}
	\boldsymbol{A}_1&		\boldsymbol{O}\\
	\boldsymbol{O}&		\boldsymbol{A}_2\\
\end{matrix} \right]$ where $\boldsymbol{A}_1$ and $\boldsymbol{A}_2$ are tridiagonal matrices with smaller size;
- Apply QR algorithm to $\boldsymbol{A}_1$ and $\boldsymbol{A}_2$ in a recursive manner;

End

Two key modifications:

1. Divide the big matrix into smaller matrices;
2. Use the shift.

### About "Shift"

How can we pick the shift?

#### Method One

**Method 1**: we can pick $\mu ^{\left( k \right)}={\boldsymbol{A}_{m,m}}^{\left( k \right)}$ as the last element of $\boldsymbol{A}^k$.

$$
{\boldsymbol{A}_{m,m}}^{\left( k \right)}={\boldsymbol{e}_m}^T\boldsymbol{A}^{\left( k \right)}\boldsymbol{e}_m={\boldsymbol{e}_m}^T\left( \boldsymbol{Q}^{\left( k \right)} \right) ^T\boldsymbol{AQ}^{\left( k \right)}\boldsymbol{e}_m
$$

$$
=\left( {\boldsymbol{q}_m}^{\left( k \right)} \right) ^T{\boldsymbol{Aq}_m}^{\left( k \right)}=\frac{\left( {\boldsymbol{q}_m}^{\left( k \right)} \right) ^T{\boldsymbol{Aq}_m}^{\left( k \right)}}{\left( {\boldsymbol{q}_m}^{\left( k \right)} \right) ^T{\boldsymbol{q}_m}^{\left( k \right)}}
$$

As $k\rightarrow +\infty$, $\boldsymbol{q}^{\left( k \right)}\rightarrow \boldsymbol{q}_m, \mu ^{\left( k \right)}\rightarrow \lambda _m$

This strategy may fail in some situations. For example, if $\boldsymbol{A}=\left[ \begin{matrix}
	0&		1\\
	1&		0\\
\end{matrix} \right]$, we get $\lambda _{1,2}=\pm 1$. $\boldsymbol{AI}=\boldsymbol{A},\boldsymbol{RQ}=\boldsymbol{A}$ remains unchanged. The shift is zero.

#### Wilkerson's Shift

**Method 2** (***Wilkerson's Shift***):

Let:

$$
\boldsymbol{B}=\left[ \begin{matrix}
	a_{m-1}&		b_{m-1}\\
	b_{m-1}&		a_m\\
\end{matrix} \right] 
$$

as a lower-right 2 by 2 matrix. Pick the eigenvalue of $\boldsymbol{B}$, closer to $a_m$, as the shift $\mu ^{\left( k \right)}$, In the case of a tie, just pick and arbitrary one:

$$
\mu ^{\left( k \right)}=a_m-\mathrm{sign}\left( \delta \right) \frac{{b_{m-1}}^2}{\left| \delta \right|+\sqrt{\delta ^2+{b_{m-1}}^2}}
$$

where $\delta =\frac{a_{m-1}-a_m}{2}$. If $\delta =0$, set $\delta =\pm 1$ arbitrarily.

### Properties

`Claims`:

- QR algorithm with Wilkerson's Shift always converges in exact arithmetic.
- QR algorithm is backward stable.
- Overall cost of QR algorithm $O\left( \frac{4}{3}m^3 \right)$ flops (mainly the cost of phase 1).

