function [mu,Sigma,fmax] = TailorProp(fun,mu,Sigma,sa_spec,cs_spec,varargin)
% Function MOVEPROB
%
% Purpose:    Tailor proposal density; see
%             Chib & Ramamurthy (2010), TaRB MCMC methods with application
%             to DSGE models, Journal of Econometrics
%
% Format:     [mu,Sigma,fmax] = TailorProp(fun,mu,Sigma,sa_spec,cs_spec,varargin)
%
% Input:      fun       function handle naming (-)log posterior kernel
%             mu        initial proposal mean vector (1 x dim)
%             Sigma     initial proposal scaling matrix (dim x dim)
%             sa_spec   simulated annealing specification
%             cs_spec   csminwel specification
%             varargin  additional inputs required by fun.m
%
% Output:     mu        optimal proposal mean vector
%             Sigma     optimal proposal scaling matrix
%             fmax      posterior kernel at mu
%
% Written by Fei Tan, Saint Louis University
% Updated: February 1, 2020

%% -------------------------------------------
%           Tailor Proposal Density
%---------------------------------------------

% Simulated annealing; first attempt (optional)
if ~isempty(sa_spec)
    mu = simulannealbnd(@(theta) fun(theta,varargin{:}),...
        mu,[],[],sa_spec);
end

% Sims' BFGS quasi-Newton; second attempt (required)
[fmax,mu,~,Sigma,~,~,~] = csminwel(fun,mu,Sigma,...
    cs_spec{1},cs_spec{2},cs_spec{3},cs_spec{4},varargin{:});
fmax = -fmax;

%-------------------- END --------------------