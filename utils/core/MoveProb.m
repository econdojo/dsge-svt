function [pk_last,pk_next,alpha] = MoveProb(fun,x,y,mu,Sigma,pdof,pk_last,pk_next,varargin)
% Function MOVEPROB
%
% Purpose:    Compute Metropolis-Hastings probability of move; see
%             Chib & Greenberg (1995), Understanding Metropolis-Hastings
%             algorithm, American Statistician
%
% Format:     [pk_last,pk_next,alpha] = MoveProb(x,y,mu,Sigma,pdof,pk_last,pk_next,varargin)
%
% Input:      fun       string naming (-)log posterior kernel density
%             x         current state vector (1 x dim)
%             y         candidate state vector (1 x dim)
%             mu        proposal mean vector (1 x dim)
%             Sigma     proposal scaling matrix (dim x dim)
%             pdof      proposal t degrees of freedom
%             pk_last   log current posterior kernel
%             pk_next   log candidate posterior kernel
%             varargin  additional inputs required by PostKer.m
%
% Output:     pk_last   log current posterior kernel
%             pk_next   log candidate posterior kernel
%             alpha     M-H probability of move
%
% Written by Fei Tan, Saint Louis University
% Updated: December 15, 2016

%% -------------------------------------------
%           M-H Probability of Move
%---------------------------------------------

% Current posterior kernel
if isempty(pk_last)
    pk_last = -feval(fun,x,varargin{:});
end

% Candidate posterior kernel
if isempty(pk_next)
    pk_next = -feval(fun,y,varargin{:});
end

% Probability of move
r = exp(pk_next-pk_last+mvt_pdf_mex(x,mu,Sigma,pdof)-mvt_pdf_mex(y,mu,Sigma,pdof));
alpha = min([r 1]);

%-------------------- END --------------------
