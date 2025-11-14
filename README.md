# nnls_based_qp_solver
Quadratic programming (QP) solver based on non-negative least squares algorithm

This project presents solver for dense quadratic problems:

min 0.5 * x^T * H * x + c^T * x
s.t. A * x <= b, lw <= x <= up

x - n x 1 - vector of decision variables 
H - n x n - must be positive definite or positive semidefinite 
c - n x 1

A - m x n - matrix of constraints
    First m_eq rows are equality constraints, rows from m_eq + 1 to m are inequality constraints
b - m x 1
lw - n x 1 - lower bounds of decision variables
up - n x 1 - upper bounds of decision variables
    lw[i] must be < up[i] for any i : 0 <= i < n
    If any x[i], 0 <= i < n have no lower or upper bound (or both) set correspoding bound with -10^20 for lower and 10^20 for upper

See example in smoke/main.cpp

The main idea of the algorithm is taken from paper:

"A Quadratic Programming Algorithm Based on
Nonnegative Least Squares With Applications
to Embedded Model Predictive Control"
Alberto Bemporad, Fellow, IEEE
IEEE TRANSACTIONS ON AUTOMATIC CONTROL, VOL. 61, NO. 4, APRIL 2016 

Many impovements have been introduced to achive good results for different types of QP problems,
including problems with bad condition number of H, and "unbalanced" problems with ||H|| << ||c|| 