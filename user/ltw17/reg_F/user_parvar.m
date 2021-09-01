% Function USER_PARVAR
%
% Purpose:    User input - parameters & variables
%
% Format:     user_parvar
%
% Input:      user
%
% Output:     para      model parameters
%             svp       stochastic volatility parameters
%             ssp       steady state/implied parameters
%             mvar      model variables
%             shock     structural shocks
%             fore      forecast errors
%             data      observables
%
% Written by Fei Tan, Saint Louis University
% Updated: November 27, 2017

%% User Input Start Here

% Note: parameter names in (para,svp,ssp) must be unique
% Parameter name, optimization initializer, lower & upper bound, prior type & parameter
para = {
           'gam100 ', 0.27, 0, 1e5, 'N', 0.40, 0.05     % [1] 100*s.s. growth rate of technology
           'xi     ', 2.25, 0, 1e5, 'G', 2.00, 0.50     % [2] inverse Frisch labor elasticity
           'theta  ', 0.96, 0, 1,   'B', 0.50, 0.20     % [3] habit formation
           'alpha_g',-.38,-1.75,1.75,'U',-1.75, 1.75     % [4] subs. of private/public consumption
           'psi    ', 0.31, 0, 1,   'B', 0.60, 0.15     % [5] capital utilization
           's      ', 3.47, 0, 1e5, 'N', 6.00, 1.50     % [6] investment adjustment cost
           'omega_p', 0.95, 0, 1,   'B', 0.50, 0.20     % [7] price stickiness
           'omega_w', 0.74, 0, 1,   'B', 0.50, 0.20     % [8] wage stickiness
           'chi_p  ', 0.11, 0, 1,   'B', 0.50, 0.20     % [9] price partial indexation
           'chi_w  ', 0.06, 0, 1,   'B', 0.50, 0.20     % [10] wage partial indexation
           'phi_p  ', 0.19, 0, 1,   'B', 0.50, 0.15     % [11] interest rate response to inflation
           'phi_y  ', 0.21, 0, 1e5, 'N', 0.125,0.05     % [12] interest rate response to output
           'rho_r  ', 0.37, 0, 1,   'B', 0.50, 0.20     % [13] interest rate lagged response
           'rho_g  ', 0.95, 0, 1,   'B', 0.50, 0.20     % [14] G lagged response
           'rho_z  ', 0.97, 0, 1,   'B', 0.50, 0.20     % [15] Z lagged response
           'rhos_a ', 0.30, 0, 1,   'B', 0.50, 0.20     % [16] technology shock AR coeff.
           'rhos_b ', 0.22, 0, 1,   'B', 0.50, 0.20     % [17] preference shock AR coeff.
           'rhos_i ', 0.47, 0, 1,   'B', 0.50, 0.20     % [18] investment shock AR coeff.
           'rhos_p ', 0.61, 0, 1,   'B', 0.50, 0.20     % [19] price markup shock AR coeff.
           'rhos_w ', 0.20, 0, 1,   'B', 0.50, 0.20     % [20] wage markup shock AR coeff.
           'rhos_m ', 0.87, 0, 1,   'B', 0.50, 0.15     % [21] monetary shock AR coeff.
           'rhos_g ', 0.29, 0, 1,   'B', 0.50, 0.15     % [22] G shock AR coeff.
           'rhos_z ', 0.90, 0, 1,   'B', 0.50, 0.15     % [23] Z shock AR coeff.
%            'sigma_a', 1.18, 0, 1e5, 'I1',2.006,0.057    % [24] 100*technology shock s.d.
%            'sigma_b',16.95, 0, 1e5, 'I1',2.006,0.057    % [25] 100*preference shock s.d.
%            'sigma_i', 1.30, 0, 1e5, 'I1',2.006,0.057    % [26] 100*investment shock s.d.
%            'sigma_p', 0.13, 0, 1e5, 'I1',2.006,0.057    % [27] 100*price markup shock s.d.
%            'sigma_w', 0.23, 0, 1e5, 'I1',2.006,0.057    % [28] 100*wage markup shock s.d.
%            'sigma_m', 0.21, 0, 1e5, 'I1',2.006,0.057    % [29] 100*monetary shock s.d.
%            'sigma_g', 2.05, 0, 1e5, 'I1',2.006,0.057    % [30] 100*G shock s.d.
%            'sigma_z', 0.78, 0, 1e5, 'I1',2.006,0.057    % [31] 100*Z shock s.d.
           'lbar   ',468.98,0, 1e5, 'N', 468,  5.00     % [32] s.s. hours worked
           'pbar   ', 0.75, 0, 1e5, 'N', 0.75, 0.25     % [33] s.s. inflation
          };

% Stochastic volatility parameters (svp = {} if constant volatility)
svp = {
           'mu_a   ', 0.50,-50,50,  'N', 0.00, 10.0     % [24]
           'mu_b   ', 6.00,-50,50,  'N', 0.00, 10.0     % [25]
           'mu_i   ', 0.50,-50,50,  'N', 0.00, 10.0     % [26]
           'mu_p   ', -3.0,-50,50,  'N', 0.00, 10.0     % [27]
           'mu_w   ', -3.0,-50,50,  'N', 0.00, 10.0     % [28]
           'mu_m   ', -3.0,-50,50,  'N', 0.00, 10.0     % [29]
           'mu_g   ', 1.50,-50,50,  'N', 0.00, 10.0     % [30]
           'mu_z   ', -0.5,-50,50,  'N', 0.00, 10.0     % [31]
           'phiv_a ', 0.90, 0, 1,   'N', 0.90, 0.05     % [34]
           'phiv_b ', 0.90, 0, 1,   'N', 0.90, 0.05     % [35]
           'phiv_i ', 0.90, 0, 1,   'N', 0.90, 0.05     % [36]
           'phiv_p ', 0.90, 0, 1,   'N', 0.90, 0.05     % [37]
           'phiv_w ', 0.90, 0, 1,   'N', 0.90, 0.05     % [38]
           'phiv_m ', 0.90, 0, 1,   'N', 0.90, 0.05     % [39]
           'phiv_g ', 0.90, 0, 1,   'N', 0.90, 0.05     % [40]
           'phiv_z ', 0.90, 0, 1,   'N', 0.90, 0.05     % [41]
           'sig2_a ', 0.05, 0, 1e5, 'I2',2.00, 0.05     % [42]
           'sig2_b ', 0.05, 0, 1e5, 'I2',2.00, 0.05     % [43]
           'sig2_i ', 0.05, 0, 1e5, 'I2',2.00, 0.05     % [44]
           'sig2_p ', 0.05, 0, 1e5, 'I2',2.00, 0.05     % [45]
           'sig2_w ', 0.05, 0, 1e5, 'I2',2.00, 0.05     % [46]
           'sig2_m ', 0.05, 0, 1e5, 'I2',2.00, 0.05     % [47]
           'sig2_g ', 0.05, 0, 1e5, 'I2',2.00, 0.05     % [48]
           'sig2_z ', 0.05, 0, 1e5, 'I2',2.00, 0.05     % [49]
          };

% Steady state/implied parameters
ssp = {
           % Calibration
           'mat'        % average maturity of government debt
           'beta'       % discount factor
           'delta'      % private capital depreciation rate
           'alpha'      % share of capital in the prod. function
           'eta_p'      % elasticity of substitution b/w intermediate goods
           'eta_w'      % elasticity of substitution b/w labor talents
           'sgy'        % steady state government consumption to model Y ratio
           'sby'        % steady state government debt to model Y ratio
           'tau_c'      % steady state consumption tax rate
           'tau_k'      % steady state capital tax rate
           'tau_l'      % steady state labor tax rate
           % Fixed
           'mu_h'       % fraction of non-savers in population
           'gam_g'      % G response to debt
           'gam_k'      % K response to debt
           'gam_l'      % L response to debt
           'gam_z'      % Z response to debt
           'rho_k'      % K lagged response
           'rho_l'      % L lagged response
           'rho_c'      % consumption tax lagged response
           % Steady state
           'rho'
           'gam'
           'expg'
           'pi_ss'
           'R_ss'       % real & nominal rate (pi = 1 in s.s.)
           'R_bar'
           'p_b'
           'R_k'
           'r_k'
           'psi1'
           'mc'
           'w'
           'kl'
           'omegl'
           'yl'
           'il'
           'cl'
           'zl'
           'c_nl'
           'c_sl'
           'c_starl'
           'l'
           'c_s'
           'c_n'
           'y'
           'k'
           'omeg'
           'c'
           'inv'
           'z'
           'b'
           'g'
           'ky'
           'cy'
           'ly'
           'S'
           'lamp'
           'lamw'
          };

% Note: variable names in (mvar,shock,fore,data) must be unique
% Model variables
mvar = {
           % Endogenous variables: x(t)
           'c_s'        % 1 consumption: Savers
           'c_n'        % 2 consumption: Non-Savers
           'R'          % 3 nominal interest rate
           'i'          % 4 investment
           'k'          % 5 effective capital
           'v'          % 6 capital utilization rate
           'l'          % 7 labor
           'y'          % 8 output 
           'g'		    % 9 govt consumption
           'c'          % 10 aggregate consumtpion
           'q'          % 11 Lagragian multiplier for investment goods				
           'r_k'        % 12 real return for private k 
           'w'          % 13 real wage 
           'pi'         % 14 inflation 
           'b'          % 15 govt debt
           'sby'        % 16 b/y ratio
           'tau_k'      % 17 tau_k
           'tau_l'      % 18 tau_l
           'tau_c'      % 19 tau_c
           'r'          % 20 real interest rate 
           'z'          % 21 transfers 
           'mc'         % 22 real marginal cost
           'k_bar'      % 23 private capital
           'lambda'     % 24 household Lagragian multiplier from budget constraint
           'p_b'        % 25 price of bonds
           'c_star'     % 26 consumption in utility function
           'S'          % 27 primary surplus
           % Expectations: E_tx(t+1)
           'E_lambda'   % 28
           'E_pi'       % 29
           'E_i'        % 30
           'E_q'        % 31
           'E_r_k'      % 32
           'E_tau_k'    % 33
           'E_w'        % 34
           'E_p_b'      % 35
           % Shocks
           'u_a'        % 36 technology shock
           'u_b'        % 37 general preference shock
           'u_i'        % 38 investment shock in adjustment costs
           'u_p'        % 39 price markup shock
           'u_w'        % 40 wage markup shock
           'u_m'        % 41 monetary policy shock
           'u_g'        % 42 G shock 
           'u_z'        % 43 Z shock
           % Dummy variables: x(t-1)
           'd_c'        % 44
           'd_i'        % 45
           'd_w'        % 46
           'd_g'        % 47
           'd_b'        % 48
          };

% Structural shocks
shock = {     
           'e_a'        % 1
           'e_b'        % 2
           'e_i'        % 3
           'e_p'        % 4
           'e_w'        % 5
           'e_m'        % 6
           'e_g'        % 7
           'e_z'        % 8
         };

% Forecast errors
fore = {
           'f_lambda'   % 1
           'f_pi'       % 2
           'f_i'        % 3
           'f_q'        % 4
           'f_r_k'      % 5
           'f_tau_k'    % 6
           'f_w'        % 7
           'f_p_b'      % 8
         };

% Observables
data = {
           'cobs'       % 1 consumption
           'iobs'       % 2 investment
           'wobs'       % 3 wages
           'gobs'       % 4 govt spending
           'bobs'       % 5 debt
           'lobs'       % 6 hours worked
           'pobs'       % 7 inflation
           'robs'       % 8 interest rate
          };

%% User Input End Here