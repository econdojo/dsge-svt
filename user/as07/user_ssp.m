% Function USER_SSP
%
% Purpose:    User input - steady state/implied parameters
%
% Format:     user_ssp
%
% Input:      P         parameter structure
%             para      model parameters
%
% Output:     ssp       steady state/implied parameters
%
% Written by Fei Tan, Saint Louis University
% Updated: December 15, 2016

%% User Input Start Here

% Steady state parameters
ssp(P.gamma) = 1+para(P.gamma_Q)/100;
ssp(P.beta) = 1/(1+para(P.r_A)/400);
ssp(P.pi_s) = 1+para(P.pi_A)/400;

%% User Input End Here