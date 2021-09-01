function stat = PostStat(chain,name)
% Function POSTSTAT
%
% Purpose:    Compute posterior statistics
%
% Format:     stat = PostStat(chain,P)
%
% Input:      chain     markov chain
%             name      parameter names (cell)
%
% Output:     stat      posterior statistics
%
% Written by Fei Tan, Saint Louis University
% Updated: January 6, 2020

%% -------------------------------------------
%            Posterior Statistics
%---------------------------------------------

% Initialization
np = size(chain,2);
stat = zeros(np,4);
prob = 0.9;

% Posterior statistics
stat(:,1) = mean(chain)';                   % 1st column: post mean
stat(:,2:3) = ProbBand(chain,prob,1);       % 2nd & 3rd columns: prob band

disp('Para             Mean          90% Interval          Ineff')
for k = 1:np
    stat(k,4) = IneffFactor(chain(:,k),0.1);% 4th column: inefficiency factor
    fprintf('%s          %.3f         [%.3f, %.3f]        %.1f\n',...
        name{k},stat(k,1),stat(k,2),stat(k,3),stat(k,4));
end

%-------------------- END --------------------