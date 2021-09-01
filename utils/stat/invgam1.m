function invgam1(nu0,mu,sigma)
% Function INVGAM1
%
% Purpose:    Compute Inverse-Gamma type-1 distribution hyperparameters
%
% Format:     invgam1(nu0,mu,sigma)
%
% Input:      nu0       initial guess; nu0>2
%             mu        mean
%             sigma     standard deviation
%
% Output:     none
%
% Written by Fei Tan, Saint Louis University
% Updated: April 23, 2017

%% -------------------------------------------
%            Inverse-Gamma Type-1
%---------------------------------------------

% Check validity of input
if nu0<=2
    error('Variance not well-defined for nu <= 2.')
end

% Find (nu,s) as in An & Schorfheide(2007) with Sims' csolve.m
[nu,~] = csolve(@ig1fun,nu0,[],1e-5,500,1,mu,sigma);
s = sqrt((mu^2+sigma^2)*(1-2/nu));

% Display results
fprintf('\n');
disp('nu       s')
fprintf('%.3f    %.3f\n',nu,s);

% Subroutine
function err = ig1fun(nu,mu,sigma)
% Copyright (C) 2011 Dynare Team
%
% This file is part of Dynare.
% Modified by Fei Tan, 4/23/2017
%
% Dynare is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% Dynare is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with Dynare.  If not, see <http://www.gnu.org/licenses/>.

% Well-defined variance
if nu<=2
    err = 1e5;
else
    err = log(2*mu^2)-log((mu^2+sigma^2)*(nu-2))+2*(gammaln(nu/2)-gammaln((nu-1)/2));
end

%-------------------- END --------------------
