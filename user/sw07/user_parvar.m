% Function USER_PARVAR
%
% Purpose:    User input - parameters & variables; see
%             Smets & Wouters (2007) for model details
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
% Updated: June 10, 2022

%% User Input Start Here

% Note: parameter names in (para,svp,ssp) must be unique
% Parameter name, optimization initializer, lower & upper bound, prior type & parameter
para = {
           % endogenous parameters
           'varphi   ', 4.00, 0, 1e5, 'N', 4.00, 1.50   % [1]
           'sigma_c  ', 1.50, 0, 1e5, 'N', 1.50, 0.37   % [2]
           'h        ', 0.70, 0, 1,   'B', 0.70, 0.10   % [3]
           'xi_w     ', 0.50, 0, 1,   'B', 0.50, 0.10   % [4]
           'sigma_l  ', 2.00, 0, 1e5, 'N', 2.00, 0.75   % [5]
           'xi_p     ', 0.50, 0, 1,   'B', 0.50, 0.10   % [6]
           'iota_w   ', 0.50, 0, 1,   'B', 0.50, 0.15   % [7]
           'iota_p   ', 0.50, 0, 1,   'B', 0.50, 0.15   % [8]
           'psi      ', 0.50, 0, 1,   'B', 0.50, 0.15   % [9]
           'Phi      ', 1.25, 1, 1e5, 'N', 1.25, 0.12   % [10]
           'r_pi     ', 1.50, 1, 1e5, 'N', 1.50, 0.25   % [11]
           'rho      ', 0.75, 0, 1,   'B', 0.75, 0.10   % [12]
           'r_y      ', 0.12, 0, 1e5, 'N', 0.12, 0.05   % [13]
           'r_delta_y', 0.12, 0, 1e5, 'N', 0.12, 0.05   % [14]
           'pi_bar   ', 0.62, 0, 1e5, 'G', 0.62, 0.10   % [15]
           'beta_hat ', 0.25, 0, 1e5, 'G', 0.25, 0.10   % [16]
           'l_bar    ', 0.00,-10,10,  'N', 0.00, 2.00   % [17]
           'gamma_bar', 0.40, 0, 1e5, 'N', 0.40, 0.10   % [18]
           'alpha    ', 0.30, 0, 1,   'N', 0.30, 0.05   % [19]
           % exogenous parameters
           'rho_a    ', 0.50, 0, 1,   'B', 0.50, 0.20   % [20]
           'rho_b    ', 0.50, 0, 1,   'B', 0.50, 0.20   % [21]
           'rho_g    ', 0.50, 0, 1,   'B', 0.50, 0.20   % [22]
           'rho_i    ', 0.50, 0, 1,   'B', 0.50, 0.20   % [23]
           'rho_r    ', 0.50, 0, 1,   'B', 0.50, 0.20   % [24]
           'rho_p    ', 0.50, 0, 1,   'B', 0.50, 0.20   % [25]
           'rho_w    ', 0.50, 0, 1,   'B', 0.50, 0.20   % [26]
           'rho_ga   ', 0.50, 0, 1,   'B', 0.50, 0.20   % [27]
           'mu_p     ', 0.50, 0, 1,   'B', 0.50, 0.20   % [28]
           'mu_w     ', 0.50, 0, 1,   'B', 0.50, 0.20   % [29]
           'sigma_a  ', 0.18, 0, 1e5, 'I1',2.00, 0.10   % [30]
           'sigma_b  ', 0.18, 0, 1e5, 'I1',2.00, 0.10   % [31]
           'sigma_g  ', 0.18, 0, 1e5, 'I1',2.00, 0.10   % [32]
           'sigma_i  ', 0.18, 0, 1e5, 'I1',2.00, 0.10   % [33]
           'sigma_r  ', 0.18, 0, 1e5, 'I1',2.00, 0.10   % [34]
           'sigma_p  ', 0.18, 0, 1e5, 'I1',2.00, 0.10   % [35]
           'sigma_w  ', 0.18, 0, 1e5, 'I1',2.00, 0.10   % [36]
           % free parameters
          };

% Stochastic volatility parameters (svp = {} if constant volatility)
svp = {
           % log(100*s.d.)^2
%            'mu_R   ', -3.0,-50,50,  'N', 0.00, 10.0     % [11]
%            'mu_G   ', -0.5,-50,50,  'N', 0.00, 10.0     % [12]
%            'mu_Z   ', -1.5,-50,50,  'N', 0.00, 10.0     % [13]
%            'phi_R  ', 0.90, 0, 1,   'N', 0.90, 0.05     % [14]
%            'phi_G  ', 0.90, 0, 1,   'N', 0.90, 0.05     % [15]
%            'phi_Z  ', 0.90, 0, 1,   'N', 0.90, 0.05     % [16]
           % variance
%            'sig2_R ', 0.20, 0, 1e5, 'I2',2.00, 0.05     % [17]
%            'sig2_G ', 0.20, 0, 1e5, 'I2',2.00, 0.05     % [18]
%            'sig2_Z ', 0.20, 0, 1e5, 'I2',2.00, 0.05     % [19]
          };

% Steady state/implied parameters
ssp = {
           'delta'
           'lambda_w'
           'g_y'
           'eps_w'
           'eps_p'
           'beta'
           'gamma'
           'pi_star'
           'r_bar'
           'r_k_ss'
           'w_ss'
           'i_k'
           'l_k'
           'k_y'
           'i_y'
           'c_y'
           'z_y'
           'wl_c'
          };

% Note: variable names in (mvar,shock,fore,data) must be unique
% Model variables
mvar = {
           'y'          % y(t)
           'c'
           'i'
           'z'
           'l'
           'r'
           'r_k'
           'q'
           'k'
           'k_s'
           'mu_p'
           'mu_w'
           'w'
           'pi'
           'y_star'     % y*(t)
           'c_star'
           'i_star'
           'z_star'
           'l_star'
           'r_star'
           'r_k_star'
           'q_star'
           'k_star'
           'k_s_star'
           'w_star'
           'eps_a'      % e_a(t)
           'eps_b'
           'eps_g'
           'eps_i'
           'eps_r'
           'eps_p'
           'eps_w'
           'E_c'        % E_tc(t+1)
           'E_i'
           'E_l'
           'E_pi'
           'E_q'
           'E_r_k'
           'E_w'
           'E_c_star'   % E_tc*(t+1)
           'E_i_star'
           'E_l_star'
           'E_q_star'
           'E_r_k_star'
           'd_eta_p'    % dummy eta_p(t)
           'd_eta_w'
           'd_y'        % dummy y(t-1)
           'd_c'
           'd_i'
           'd_w'
          };

% Structural shocks
shock = {
           'eta_a'
           'eta_b'
           'eta_g'
           'eta_i'
           'eta_r'
           'eta_p'
           'eta_w'
         };

% Forecast errors
fore = {
           'fe_c'       % c(t)-E_{t-1}c(t)
           'fe_i'
           'fe_l'
           'fe_pi'
           'fe_q'
           'fe_r_k'
           'fe_w'
           'fe_c_star'
           'fe_i_star'
           'fe_l_star'
           'fe_q_star'
           'fe_r_k_star'
         };

% Observables
data = {
           'YGR'        % per capita real output growth
           'CGR'        % per capita real consumption growth
           'IGR'        % per capita real investment growth
           'WGR'        % per capita real wage growth
           'HOUR'       % per capita hours index
           'INF'        % inflation rate
           'FFR'        % federal funds rate
          };

%% User Input End Here