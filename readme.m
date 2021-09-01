% *--------------------------------*
% |    DSGE-SVt MATLAB Toolbox     |
% |                                |
% |   Programmed by Fei Tan        |
% |   Last modified: 01-Feb-2020   |
% *--------------------------------*
%
% I. GENERAL INFORMATION
%
% The MATLAB(R) programs in this package performs the Bayesian analysis of
% high-dimensional DSGE models with Student-t shocks and stochastic volatility.
% It simulates the posterior distribution by the tailored randomized block
% Metropolis-Hastings (TaRB-MH) algorithm of Chib & Ramamurphy (2010), calculates
% the marginal likelihood by the method of Chib (1995) and Chib & Jeliazkov (2001),
% and includes various post-estimation tools. See our accompanying paper,
% Chib, Shin & Tan (2020): 'DSGE-SVt: An Econometric Toolkit for High-Dimenionsal
% DSGE Models with SV and t Errors.'
% 
% II. CONTENTS
%
% The DSGE-SVt/ directory contains a set of .m programs.
% See also:
%    tarb.m
%    sample_prior.m
%    sample_post.m
%    marg_lik.m
%    find_mode.m
%    demo.m
%
% The directory also contains a set of subfolders:
%    mex/          - C source MEX functions.
%    user/         - data file and a set of user input .m programs.
%    utils/        - a set of routines called by the above .m programs.
% 
% III. EXAMPLE
%
%    1. Set MATLAB working directory to the DSGE-SVt/ folder.
%
%    2. Include data (.txt file with ROW observations) and user input .m
%       programs in the user/ subfolder.
%
%    3. Run the demo file 'demo.m'.
%
%    4. Note: results will automatically be saved under the user/ subfolder.
% 
% IV. CONTACT
%
% If you experience any difficulty using this package or have any comment
% for improvement, please feel free to contact me via e-mail:
%
% Fei Tan
% tanf@slu.edu
% 
% V. REVISION NOTES
%
% 01-Feb-2020 update:
%    1. This is the latest version. Please check back for new updates.
%
% *-------------------- END --------------------*
clc
help readme