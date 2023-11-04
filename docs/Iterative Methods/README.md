# Iterative Methods

We introduce some **iterative methods** for solving the linear system $\boldsymbol{Ax}=\boldsymbol{b}$ in this chapter.

Why do we need iterative methods?

- Reduce the cost of direct methods: $O\left( m^3 \right) \rightarrow O\left( m^2 \right) \xrightarrow{\mathrm{sometimes}}O\left( m \right)$;
- Take advantage of matrix-vector multiplications:
  - Matrix-vector multiplications may be inexpensive if the matrix is structured (e.g., banned, sparse);
  - Iterative methods are especially useful if $\boldsymbol{A}$ is not known explicitly (e.g., black box, PDE).

[Classical Methods](./Classical_Iteration.md)
