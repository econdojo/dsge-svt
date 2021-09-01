function DiagPlot(chain,name,savepath)
% Function DIAGPLOT
%
% Purpose:    Plot convergence diagnostic
%
% Format:     DiagPlot(chain,P,savepath)
%
% Input:      chain     MCMC posterior draws
%             name      parameter names (cell)
%             savepath  save path
%
% Output:     none
%
% Written by Fei Tan, Saint Louis University
% Updated: December 15, 2016

%% -------------------------------------------
%         Convergence Diagnostic Plot
%---------------------------------------------

% Initialization
[M,np] = size(chain);                  % sample size & number of parameters

% Determine plot dimension
if np<=4
    dim1 = 1;
elseif np>=5 && np<=8
    dim1 = 2;
elseif np>=9 && np<=12
    dim1 = 3;
elseif np>=13
    dim1 = 4;
end
dim2 = ceil(np/dim1);

% Plot settings
set(0,'DefaultAxesFontName', 'Times New Roman')
set(0,'DefaultAxesFontSize',12)

% Plot autocorrelation function
figure(1)
for k = 1:np
    subplot(dim1,dim2,k);
    autocorr(chain(:,k),500)
    set(gca,'xlabel',[],'ylabel',[])
    title(strtrim(name{k}),'fontsize',12,'fontname','times');
end
saveas(gcf,[savepath filesep 'Fig_ACF.fig']);

% Traceplot
figure(2)
for k = 1:np
    subplot(dim1,dim2,k);
    hold on; box on;
    plot(1:M,chain(:,k),'linewidth',1.5);
    title(strtrim(name{k}),'fontsize',12,'fontname','times');
end
saveas(gcf,[savepath filesep 'Fig_Trace.fig']);

%-------------------- END --------------------
