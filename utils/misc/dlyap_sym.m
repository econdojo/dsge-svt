function [x,u]=dlyap_sym(a,b,qz_criterium)
% Taken from Dynare, Version 4
% Solves the Lyapunov equation x-a*x*a' = b, for b (and then x) symmetrical
% if a has some unit roots, the function computes only the solution of the stable subsystem
%  
% INPUTS:
%   a:      coefficient matrix (n x n)  
%   b:      coefficient square matrix (n x n)
%
% OUTPUTS
%   x:      solution matrix (m x m)
%   ns_var: vector of indices of non-stationary variables (p x 1)
%           (m + p = n)
%   u:      Schur vectors associated with unit roots  
%
% ALGORITHM
%   uses reordered Schur decomposition
%
% SPECIAL REQUIREMENTS
%   needs Matlab >= 7.0.1 for ordeig function

% Copyright (C) 2006-2008 Dynare Team
%
% This file is part of Dynare.
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

if size(a,1) == 1
    x=b/(1-a*a);
    return
end
[u,t] = schur(a);
e1 = abs(ordeig(t)) > 2-qz_criterium;
k = sum(e1);

% selects stable roots
[u,t] = ordschur(u,t,e1); 
n = length(e1)-k;
b=u(:,k+1:end)'*b*u(:,k+1:end);
t = t(k+1:end,k+1:end);

x=zeros(n,n);
i = n;
while i >= 2
    if t(i,i-1) == 0
        if i == n
            c = zeros(n,1);
        else
            c = t(1:i,:)*(x(:,i+1:end)*t(i,i+1:end)')+...
            t(i,i)*t(1:i,i+1:end)*x(i+1:end,i);
        end
        q = eye(i)-t(1:i,1:i)*t(i,i);
        x(1:i,i) = q\(b(1:i,i)+c);
        x(i,1:i-1) = x(1:i-1,i)';
        i = i - 1;
    else
        if i == n
            c = zeros(n,1);
            c1 = zeros(n,1);
        else
            c = t(1:i,:)*(x(:,i+1:end)*t(i,i+1:end)')+...
            t(i,i)*t(1:i,i+1:end)*x(i+1:end,i)+...
            t(i,i-1)*t(1:i,i+1:end)*x(i+1:end,i-1);
            c1 = t(1:i,:)*(x(:,i+1:end)*t(i-1,i+1:end)')+...
            t(i-1,i-1)*t(1:i,i+1:end)*x(i+1:end,i-1)+...
            t(i-1,i)*t(1:i,i+1:end)*x(i+1:end,i);
        end
        q = [eye(i)-t(1:i,1:i)*t(i,i) -t(1:i,1:i)*t(i,i-1);...
	      -t(1:i,1:i)*t(i-1,i) eye(i)-t(1:i,1:i)*t(i-1,i-1)];
        z =  q\[b(1:i,i)+c;b(1:i,i-1)+c1];
        x(1:i,i) = z(1:i);
        x(1:i,i-1) = z(i+1:end);
        x(i,1:i-1)=x(1:i-1,i)';
        x(i-1,1:i-2)=x(1:i-2,i-1)';
        i = i - 2;
    end
end

if i == 1
    c = t(1,:)*(x(:,2:end)*t(1,2:end)')+t(1,1)*t(1,2:end)*x(2:end,1);
    x(1,1)=(b(1,1)+c)/(1-t(1,1)*t(1,1));
end
x = u(:,k+1:end)*x*u(:,k+1:end)';
u = u(:,1:k);