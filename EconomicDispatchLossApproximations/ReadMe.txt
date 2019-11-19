To run this code your MATLAB directory should be set to the main folder ``EconomicDispatchLossApproximations''.  Also, unzip the file ``matpower6.0.zip'' in this folder. 

This code uses the MATPOWER toolbox, the NESTA test case archive, and other standard MATLAB functions.

`maincode.m' is the main file that should be run.  It produces the results used in Table I of the ACC paper entitled ``A General Economic Dispatch Problem with Marginal Losses.''

`data.m' is a function I created to produce standard test case data including admittance matrices.

`constraints.m' is a function that represents the non-linear constraints of the economic dispatch problems.

`deltathetalimits.m' is a function that produces the current magnitude limits of each line as discussed in Section III-C2 in the paper.

`ISFuncFrom.m' and `ISFuncTo.m' represent the current magnitude functions at both ends of the transmission line.  

`ISFuncFromSolve.m' and `ISFuncToSolve.m' are functions used to solve for the current magnitude limits.

`LossFunc.m' is a function representing the loss function.

`MPFFunc.m' is the mid-line power flow function as specified in the paper.

`myhessian.m' is a function that represents the analytical hessian of the Lagrangian of each economic dispatch problem.

`cost.m' is the cost function of the economic dispatch problems
 