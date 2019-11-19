To run this code your MATLAB directory should be set to the main folder ``LossLinearizationCode''. Also, unzip the file ``matpower6.0.zip'' in this folder. 

This code uses the MATPOWER toolbox, the NESTA test case archive, and the MATLAB Optimization Toolbox.
Note: The original code used mosek to solve most of the optimization problems.  To simplify the setup procedures, this code only uses the MATLAB Optimization Toolbox.  For this reason, the results might be slightly different than in the paper. 
Note: The results in Table I are not provided in this code.  To produce these results I used an interior point algorithm from the Knitro software package, which requires a license to be purchased.  For this reason I removed that code. 

`maincode.m' is the main file that should be run.  It produces the results used in Table II of the IEEE TPWRS paper entitled ``Approximating Economic Dispatch by Linearizing Transmission Losses.''

`data.m' is a function I created to produce standard test case data including admittance matrices.

`constraints.m' is a function that represents the non-linear constraints of the economic dispatch problems.

`LossFunc.m' is a function representing the loss function.

`MPFFunc.m' is the mid-line power flow function as specified in the paper.

`myhessian.m' is a function that represents the analytical hessian of the Lagrangian of each economic dispatch problem.

`cost.m' is the cost function of the economic dispatch problems
 