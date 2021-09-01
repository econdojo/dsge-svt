function index = Resample(weight,strategy)
% Function RESAMPLE
%
% Purpose:    Resampling in particle filter; see
%             Arnaud Doucet & Adam Johansen (2014), A tutorial on particle
%             filtering and smoothing: fifteen years later
%
% Format:     index = Resample(weight,strategy)
%
% Input:      weight    filtered weights
%             strategy  resampling strategy
%
% Output:     index     resampled indices
%
% Programmed by Diego Andres Alvarez Marin
% Modified by Fei Tan, Saint Louis University
% Updated: December 15, 2016

%% -------------------------------------------
%             Resampling Strategy
%---------------------------------------------

% Initialization
N = length(weight);                    % number of particles

% Choose resampling strategy
switch strategy
    case 'multinomial'
        index = randsample(1:N,N,true,weight);
%{
        This is equivalent to:
        edge = min([0 cumsum(weight)'],1);  % protect against accumulated round-off
        edge(end) = 1;                 % get the upper edge exact
        [~,index] = histc(sort(rand(N,1)),edge);
%}
    case 'systematic'
        % Latin hypercube sampling on weight
        edge = min([0 cumsum(weight)'],1);  % protect against accumulated round-off
        edge(end) = 1;                 % get the upper edge exact
        [~,index] = histc(rand/N+(0:1/N:1-1/N),edge);
    %case 'regularized'                 % to be added
    %case 'stratified'                  % to be added
    %case 'residual'                    % to be added
    otherwise
        error('Resampling strategy does not exist.')
end

%-------------------- END --------------------
