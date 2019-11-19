# TransmissionLossModels
There are two main folders in this repository. Both folders stand alone and their code should be run separately. 

The folder "EconomicDispatchLossApproximations" is based on the ACC paper entitled "A General Economic Dispatch Problem with Marginal Losses."  This code specifically produces the results from table I of this paper.  This folder should be set as the MATLAB directory when running this code.  Also, the file "matpower6.0.zip" must be unzipped.

The folder "LossLinearizationCode" is based on the IEEE Trans. on Power Systems paper entitled "Approximating Economic Dispatch by Linearizing Transmission Losses."  This code specifically produces the results from table II of this paper. This folder should be set as the MATLAB directory when running this code.  Also, the file "matpower6.0.zip" must be unzipped.

The code in both folders has been simplified.  The original code used more efficient toolboxes such as Knitro and MOSEK.  To simplify the setup, this code only uses the standard Optimization Toolbox in MATLAB.  As a result, some of the results in the papers may differ from the results in this code.
