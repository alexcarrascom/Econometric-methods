# SVAR package

Read SVAR.pdf for a nice application: reduced form estimation, structural estimation, and confidence intervals for the impulse response functions were all computing in 
MATLAB. Basically, I have created five script functions that can be used to a wide range of structural VAR models with short and 
long-run identification restrictions:

1. mlSVAR.m: Computes quasi maximum likelihood estimator for a structural VAR with A-model residual structure whenever identification 
             restrictions can be expressed as $\mathbf{R_A}\pmb{\gamma}_A=vec(\mathbf{A})$.
2. VARest.m: Computes OLS estimators for parameters in a VAR(p) reduced form model.
3. VAR_boots.m: Implements model-based bootstrap method for dependent data.  
4. VAROptLag.m: Selects the lag order by minimization of information statistics.
5. VARimpulse.m: Computes impulse response functions using companion form of the estimated VAR(p) model. This functions needs an 
                 estimation for matrix A.
6. uroottest.m: Simplifies the process of executing ADF and Phillips-Perron unit root test with different lag order for univariate 
                time series, selecting the best model based on BIC statistic, and reporting results. Econometric Toolbox must be 
                added to MATLAB search path.

By Alex Carrasco/ alex.carmar93@gmail.com

