function [p,m,v] = mix10
% Function MIX10
%
% Purpose:    10-component mixture normal distribution; see
%             Omori Chib Shephard & Nakajima (2007), Stochastic Volatility
%             with Leverage: Fast and Efficient Likelihood Inference, JOE
%
% Format:     [p,m,v] = mix10
%
% Input:      none
%
% Output:     p(k)      probability of component-k
%             m(k)      mean of component-k
%             v(k)      variance of component-k
%
% Written by Fei Tan, Saint Louis University
% Updated: January 6, 2020

% Probability
p = [0.00609 0.04775 0.13057 0.20674 0.22715 0.18842 0.12047 0.05591 0.01575 0.00115];

% Mean
m = [1.92677 1.34744 0.73504 0.02266 -0.85173 -1.97278 -3.46788 -5.55246 -8.68384 -14.65000];

% Variance
v = [0.11265 0.17788 0.26768 0.40611 0.62699 0.98583 1.57469 2.54498 4.16591 7.33342];