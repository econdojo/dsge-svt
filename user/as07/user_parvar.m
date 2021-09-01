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
% Updated: December 15, 2016

%% User Input Start Here

% Note: parameter names in (para,svp,ssp) must be unique
% Parameter name, optimization initializer, lower & upper bound, prior type & parameter
para = {
           'tau    ', 2.00, 0, 1e5, 'G', 2.00, 0.50     % [1]
           'kappa  ', 0.15, 0, 1e5, 'G', 0.20, 0.10     % [2]
           'psi_1  ', 1.50, 1, 1e5, 'G', 1.50, 0.25     % [3]
           'psi_2  ', 0.50, 0, 1e5, 'G', 0.50, 0.25     % [4]
           'r_A    ', 0.40, 0, 1e5, 'G', 0.50, 0.25     % [5]
           'pi_A   ', 4.00, 0, 1e5, 'G', 7.00, 2.00     % [6]
           %'gamma_Q', 0.50, 0, 1e5, 'N', 0.40, 0.20     % [7]
           'gamma_Q', 0.50, 0, 1e5, 't', 0.40, 15.0     % [7]
           'rho_R  ', 0.70, 0, 1,   'B', 0.50, 0.20     % [8]
           'rho_G  ', 0.90, 0, 1,   'B', 0.50, 0.20     % [9]
           'rho_Z  ', 0.80, 0, 1,   'B', 0.50, 0.20     % [10]
           % constant volatility (100*s.d.)
%            'sigma_R', 0.20, 0, 1e5, 'I1',2.00, 0.10     % [11]
%            'sigma_G', 0.80, 0, 1e5, 'I1',2.00, 0.10     % [12]
%            'sigma_Z', 0.50, 0, 1e5, 'I1',2.00, 0.10     % [13]
          };

% Stochastic volatility parameters (svp = {} if constant volatility)
svp = {
           % log(100*s.d.)^2
           'mu_R   ', -3.0,-50,50,  'N', 0.00, 10.0     % [11]
           'mu_G   ', -0.5,-50,50,  'N', 0.00, 10.0     % [12]
           'mu_Z   ', -1.5,-50,50,  'N', 0.00, 10.0     % [13]
           'phi_R  ', 0.90, 0, 1,   'N', 0.90, 0.05     % [14]
           'phi_G  ', 0.90, 0, 1,   'N', 0.90, 0.05     % [15]
           'phi_Z  ', 0.90, 0, 1,   'N', 0.90, 0.05     % [16]
           % variance
           'sig2_R ', 0.20, 0, 1e5, 'I2',2.00, 0.05     % [17]
           'sig2_G ', 0.20, 0, 1e5, 'I2',2.00, 0.05     % [18]
           'sig2_Z ', 0.20, 0, 1e5, 'I2',2.00, 0.05     % [19]
          };

% Steady state/implied parameters
ssp = {
           'gamma'
           'beta'
           'pi_s'
          };

% Note: variable names in (mvar,shock,fore,data) must be unique
% Model variables
mvar = {
           'y'          % y(t)
           'pi'
           'R'
           'g'
           'z'
           'E_y'        % E_ty(t+1)
           'E_pi'
           'd_y'        % dummy y(t-1)
          };

% Structural shocks
shock = {
           'eps_R'      % eps_R(t)
           'eps_G'
           'eps_Z'
         };

% Forecast errors
fore = {
           'fe_y'       % y(t)-E_{t-1}y(t)
           'fe_pi'
         };

% Observables
data = {
           'YGR'        % per capita real GDP growth rate
           'INF'        % annualized inflation rate
           'INT'        % annualized nominal interest rate
          };

%% User Input End Here