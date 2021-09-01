function irfplot(SSR,plotvar,V,T)
% Function IRFPLOT
%
% Purpose:    Plot impulse response function
%
% Format:     irfplot(SSR,plotvar,V,T)
%
% Input:      SSR       state space representation (structure)
%             plotvar   variables to plot
%             V         variable (structure)
%             T         number of periods (quarters)
%
% Output:     none
%
% Programmed by Alex Richter & Nate Throckmorton
% Modified by Fei Tan, Saint Louis University
% Updated: December 15, 2016

%% -------------------------------------------
%            Plot Impulse Response           
%---------------------------------------------

% Initialization
shock = diag(sqrt(SSR.Sigma_e(:,1)));  % one s.d. shock
period = 0:(T-1);                      % use (0:(T-1))/4 for yearly model
IRF = zeros(V.mod.nvar,V.mod.nshock,T);% impulse response function
nplot = length(plotvar);               % number of plot variables
plotpos = zeros(nplot,1);              % plot position

% Percentage deviation from steady state
IRF(:,:,1) = SSR.M*shock;
for k = 2:T
   IRF(:,:,k)= SSR.G*IRF(:,:,k-1);
end
IRF(abs(IRF)<1e-10) = 0;          % eliminate very small responses

% Determine plot position
for k = 1:nplot
    plotpos(k) = strmatch(plotvar(k),V.mod.var,'exact'); %#ok<MATCH3>
end    

% Determine plot dimension
if nplot<=3
    dim1 = 1;
elseif nplot>=4 && nplot<=8
    dim1 = 2;
elseif nplot>=9
    dim1 = 3;
end
dim2 = ceil(nplot/dim1);

% Plot setting
set(0,'DefaultAxesFontName', 'Times New Roman')
set(0,'DefaultAxesFontSize',12)

% Plot
for j = 1:V.mod.nshock
    figure(j)
    for k = 1:nplot
        subplot(dim1,dim2,k)
        hold on; grid on; box on;
        IRFper(1,:) = IRF(plotpos(k),j,:);
        plot(period,IRFper,'-.r','MarkerSize',3,'linewidth',2)
        title(plotvar(k),'fontsize',12,'fontname','times')
        ylabel('%')
    end
    legend(V.mod.shock(j),'fontsize',12,'fontname','times')
end

%-------------------- END --------------------
