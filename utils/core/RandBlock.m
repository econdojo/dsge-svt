function block = RandBlock(ind,prob_nb)
% Function RANDBLOCK
%
% Purpose:    Generate random blocks of parameters; see
%             Chib & Ramamurthy (2010), TaRB MCMC methods with application
%             to DSGE models, Journal of Econometrics
%
% Format:     block = RandBlock(ind,prob_nb)
%
% Input:      ind       parameter indices
%             prob_nb   blocking probability; number of blocks if >=1 integer
%
% Output:     block     random blocks of parameters (cell w/ row vecs)
%
% Written by Fei Tan, Saint Louis University
% Updated: June 15, 2017

%% -------------------------------------------
%             Parameter Blocking
%---------------------------------------------

% Initialization
if ~isempty(ind)
    npara = length(ind);          % number of parameters
    ind = ind(randperm(npara));   % random permutation
else
    block = [];
    return
end

% Choose blocking method
if prob_nb<1                      % random blocks
    j = 1;
    block{j} = ind(1);
    for k = 2:npara
        if rand<=prob_nb
            block{j} = cat(2,block{j},ind(k));
            continue
        end
        j = j+1;
        block{j} = ind(k);
    end
else                              % almost equal-sized blocks
    block = cell(prob_nb,1);
    bi = ones(1,prob_nb)*floor(npara/prob_nb);
    bi(1:mod(npara,prob_nb)) = bi(1:mod(npara,prob_nb))+1;
    bi = [1 1+cumsum(bi)];
    for k = 1:prob_nb
        block{k} = ind(bi(k):(bi(k+1)-1));
    end
end

%-------------------- END --------------------
