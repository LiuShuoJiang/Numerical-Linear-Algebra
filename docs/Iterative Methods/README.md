# Iterative Methods

## General Introduction

We introduce some **iterative methods** for solving the linear system $\boldsymbol{Ax}=\boldsymbol{b}$ in this chapter.

Why do we need iterative methods?

- Reduce the cost of direct methods: $O\left( m^3 \right) \rightarrow O\left( m^2 \right) \xrightarrow{\mathrm{sometimes}}O\left( m \right)$;
- Take advantage of matrix-vector multiplications:
  - Matrix-vector multiplications may be inexpensive if the matrix is structured (e.g., banned, sparse);
  - Iterative methods are especially useful if $\boldsymbol{A}$ is not known explicitly (e.g., black box, PDE).

## Contents in This Chapter

- [Classical Methods](./Classical_Iteration.md)
- [General Strategy of Iteration Methods](./Introduction_to_Iteration.md)
- [Projection Method](./Projection_Method.md)
- [Steepest Descent](./Steepest_Descent.md)
- [Conjugate Gradient 1](./Conjugate_Gradient.md)
- [Conjugate Gradient 2](./Conjugate_Gradient_2.md)
- [GMRES](./GMRES.md)
- [Preconditioning](./Preconditioning.md)
