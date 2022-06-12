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
%             ssp       steady state/implied parameters
%             mvar      model variables
%             mshock    model structural shocks
%             mfore     model forecast errors
%             dvar      data variables
%
% Written by Fei Tan, Saint Louis University
% Updated: December 15, 2016

%% User Input Start Here

% Parameter names, optimization initializer, lower & upper bounds
% Prior types & parameters, plot or not
% Original prior
para = {
% Endogenous parameters
           'varphi   ', 4.00, 0, 1e5, 'N', 4.00, 1.50, 1   % [1]
           'sigma_c  ', 1.50, 0, 1e5, 'N', 1.50, 0.37, 1   % [2]
           'h        ', 0.70, 0, 1,   'B', 0.70, 0.10, 1   % [3]
           'xi_w     ', 0.50, 0, 1,   'B', 0.50, 0.10, 1   % [4]
           'sigma_l  ', 2.00, 0, 1e5, 'N', 2.00, 0.75, 1   % [5]
           'xi_p     ', 0.50, 0, 1,   'B', 0.50, 0.10, 1   % [6]
           'iota_w   ', 0.50, 0, 1,   'B', 0.50, 0.15, 1   % [7]
           'iota_p   ', 0.50, 0, 1,   'B', 0.50, 0.15, 1   % [8]
           'psi      ', 0.50, 0, 1,   'B', 0.50, 0.15, 1   % [9]
           'Phi      ', 1.25, 1, 1e5, 'N', 1.25, 0.12, 1   % [10]
           'r_pi     ', 1.50, 1, 1e5, 'N', 1.50, 0.25, 1   % [11]
           'rho      ', 0.75, 0, 1,   'B', 0.75, 0.10, 1   % [12]
           'r_y      ', 0.12, 0, 1e5, 'N', 0.12, 0.05, 1   % [13]
           'r_delta_y', 0.12, 0, 1e5, 'N', 0.12, 0.05, 1   % [14]
           'pi_bar   ', 0.62, 0, 1e5, 'G', 0.62, 0.10, 1   % [15]
           'beta_hat ', 0.25, 0, 1e5, 'G', 0.25, 0.10, 1   % [16]
           'l_bar    ', 0.00,-10,10,  'N', 0.00, 2.00, 1   % [17]
           'gamma_bar', 0.40, 0, 1e5, 'N', 0.40, 0.10, 1   % [18]
           'alpha    ', 0.30, 0, 1,   'N', 0.30, 0.05, 1   % [19]
% Exogenous parameters
           'rho_a    ', 0.50, 0, 1,   'B', 0.50, 0.20, 0   % [20]
           'rho_b    ', 0.50, 0, 1,   'B', 0.50, 0.20, 0   % [21]
           'rho_g    ', 0.50, 0, 1,   'B', 0.50, 0.20, 0   % [22]
           'rho_i    ', 0.50, 0, 1,   'B', 0.50, 0.20, 0   % [23]
           'rho_r    ', 0.50, 0, 1,   'B', 0.50, 0.20, 0   % [24]
           'rho_p    ', 0.50, 0, 1,   'B', 0.50, 0.20, 0   % [25]
           'rho_w    ', 0.50, 0, 1,   'B', 0.50, 0.20, 0   % [26]
           'rho_ga   ', 0.50, 0, 1,   'B', 0.50, 0.20, 0   % [27]
           'mu_p     ', 0.50, 0, 1,   'B', 0.50, 0.20, 0   % [28]
           'mu_w     ', 0.50, 0, 1,   'B', 0.50, 0.20, 0   % [29]
           'sigma_a  ', 0.18, 0, 1e5, 'I1',2.00, 0.10, 0   % [30]
           'sigma_b  ', 0.18, 0, 1e5, 'I1',2.00, 0.10, 0   % [31]
           'sigma_g  ', 0.18, 0, 1e5, 'I1',2.00, 0.10, 0   % [32]
           'sigma_i  ', 0.18, 0, 1e5, 'I1',2.00, 0.10, 0   % [33]
           'sigma_r  ', 0.18, 0, 1e5, 'I1',2.00, 0.10, 0   % [34]
           'sigma_p  ', 0.18, 0, 1e5, 'I1',2.00, 0.10, 0   % [35]
           'sigma_w  ', 0.18, 0, 1e5, 'I1',2.00, 0.10, 0   % [36]
          };

% Diffuse prior = original prior S.D. x 3 (except for *)
% para = {
% % Endogenous parameters
%            'varphi   ', 4.00, 0, 1e5, 'N', 4.00, 4.50, 1   % [1]
%            'sigma_c  ', 1.50, 0, 1e5, 'N', 1.50, 1.11, 1   % [2]
%            'h        ', 0.70, 0, 1,   'U', 0.00, 1.00, 1   % [3]
%            'xi_w     ', 0.50, 0, 1,   'U', 0.00, 1.00, 1   % [4]
%            'sigma_l  ', 2.00, 0, 1e5, 'N', 2.00, 2.25, 1   % [5]
%            'xi_p     ', 0.50, 0, 1,   'U', 0.00, 1.00, 1   % [6]
%            'iota_w   ', 0.50, 0, 1,   'B', 0.50, 0.15, 1   % [7]  *
%            'iota_p   ', 0.50, 0, 1,   'U', 0.00, 1.00, 1   % [8]
%            'psi      ', 0.50, 0, 1,   'U', 0.00, 1.00, 1   % [9]
%            'Phi      ', 1.25, 1, 1e5, 'N', 1.25, 0.36, 1   % [10]
%            'r_pi     ', 1.50, 1, 1e5, 'N', 1.50, 0.75, 1   % [11]
%            'rho      ', 0.75, 0, 1,   'U', 0.00, 1.00, 1   % [12]
%            'r_y      ', 0.12, 0, 1e5, 'N', 0.12, 0.15, 1   % [13]
%            'r_delta_y', 0.12, 0, 1e5, 'N', 0.12, 0.15, 1   % [14]
%            'pi_bar   ', 0.62, 0, 1e5, 'G', 0.62, 0.30, 1   % [15]
%            'beta_hat ', 0.25, 0, 1e5, 'G', 0.25, 0.10, 1   % [16] *
%            'l_bar    ', 0.00,-10,10,  'N', 0.00, 6.00, 1   % [17]
%            'gamma_bar', 0.40, 0, 1e5, 'N', 0.40, 0.30, 1   % [18]
%            'alpha    ', 0.30, 0, 1,   'N', 0.30, 0.15, 1   % [19]
% % Exogenous parameters
%            'rho_a    ', 0.50, 0, 1,   'U', 0.00, 1.00, 0   % [20]
%            'rho_b    ', 0.50, 0, 1,   'U', 0.00, 1.00, 0   % [21]
%            'rho_g    ', 0.50, 0, 1,   'U', 0.00, 1.00, 0   % [22]
%            'rho_i    ', 0.50, 0, 1,   'U', 0.00, 1.00, 0   % [23]
%            'rho_r    ', 0.50, 0, 1,   'U', 0.00, 1.00, 0   % [24]
%            'rho_p    ', 0.50, 0, 1,   'U', 0.00, 1.00, 0   % [25]
%            'rho_w    ', 0.50, 0, 1,   'U', 0.00, 1.00, 0   % [26]
%            'rho_ga   ', 0.50, 0, 1,   'U', 0.00, 1.00, 0   % [27]
%            'mu_p     ', 0.50, 0, 1,   'U', 0.00, 1.00, 0   % [28]
%            'mu_w     ', 0.50, 0, 1,   'U', 0.00, 1.00, 0   % [29]
%            'sigma_a  ', 0.18, 0, 1e5, 'I1',2.00, 0.10, 0   % [30] *
%            'sigma_b  ', 0.18, 0, 1e5, 'I1',2.00, 0.10, 0   % [31] *
%            'sigma_g  ', 0.18, 0, 1e5, 'I1',2.00, 0.10, 0   % [32] *
%            'sigma_i  ', 0.18, 0, 1e5, 'I1',2.00, 0.10, 0   % [33] *
%            'sigma_r  ', 0.18, 0, 1e5, 'I1',2.00, 0.10, 0   % [34] *
%            'sigma_p  ', 0.18, 0, 1e5, 'I1',2.00, 0.10, 0   % [35] *
%            'sigma_w  ', 0.18, 0, 1e5, 'I1',2.00, 0.10, 0   % [36] *
%           };

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

% Model structural shocks
mshock = {
           'eta_a'
           'eta_b'
           'eta_g'
           'eta_i'
           'eta_r'
           'eta_p'
           'eta_w'
         };

% Model forecast errors
mfore = {
           'c'          % c(t)-E_{t-1}c(t)
           'i'
           'l'
           'pi'
           'q'
           'r_k'
           'w'
           'c_star'
           'i_star'
           'l_star'
           'q_star'
           'r_k_star'
         };

% Data variables
dvar = {
           'YGR'        % per capita real output growth
           'CGR'        % per capita real consumption growth
           'IGR'        % per capita real investment growth
           'WGR'        % per capita real wage growth
           'HOUR'       % per capita hours index
           'INF'        % inflation rate
           'FFR'        % federal funds rate
          };

%% User Input End Here