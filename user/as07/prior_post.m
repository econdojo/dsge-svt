clear
close all
clc

savepath = ['user' filesep 'as07'];
load([savepath filesep '1' filesep 'tarb_full.mat'],'chain_para')
chain_para1 = chain_para;
load([savepath filesep '2' filesep 'tarb_full.mat'],'chain_para')
chain_para2 = chain_para;
load([savepath filesep '2' filesep 'tarb_prior.mat'],'prior_para')

para = {
           '$\tau$', 2.00, 0, 4, 'G', 2.00, 0.50     % [1]
           '$\kappa$', 0.15, 0, 2, 'G', 0.20, 0.10     % [2]
           '$\psi_1$', 1.50, 1, 3, 'G', 1.50, 0.25     % [3]
           '$\psi_2$', 1.00, 0, 2, 'G', 0.50, 0.25     % [4]
           '$r^{(A)}$', 0.40, 0, 2, 'G', 0.50, 0.25     % [5]
           '$\pi^{(A)}$', 4.00, 0, 8, 'G', 7.00, 2.00     % [6]
           '$\gamma^{(Q)}$', 0.50, 0, 2, 'N', 0.40, 0.20     % [7]
           '$\rho_R$', 0.60, 0, 1,   'B', 0.50, 0.20     % [8]
           '$\rho_g$', 0.95, 0, 1,   'B', 0.80, 0.10     % [9]
           '$\rho_z$', 0.65, 0, 1,   'B', 0.66, 0.15     % [10]
          };

% Plot
figure('Name','Prior-Posterior')
[ha,~] = tight_subplot(3,4,[.05 .05],[.03 .03],[.03 .02]);
delete(ha(11));
delete(ha(12));
set(gcf,'color','w');

for i = 1:10
    axes(ha(i));
    pd = fitdist(chain_para1(:,i),'Kernel','Kernel','epanechnikov');
    x = linspace(para{i,3},para{i,4});
    d = pdf(pd,x);
    ciplot(zeros(1,100),d,x,[0.65 0.65 0.65]);
    hold on
    pd = fitdist(chain_para2(:,i),'Kernel','Kernel','epanechnikov');
    x = linspace(para{i,3},para{i,4});
    d = pdf(pd,x);
    ciplot(zeros(1,100),d,x,[0.35 0.35 0.35]);
    xline(para{i,2},'Color','black','LineStyle','-','linewidth',0.8);
    pd = fitdist(prior_para(:,i),'Kernel','Kernel','epanechnikov');
    d = pdf(pd,x);
    plot(x,d,'--r','LineWidth',0.8);
    xlim([para{i,3},para{i,4}]);
    title(para{i,1},'Interpreter','latex','FontSize',12);
end

x0=50;
y0=100;
width = 650;
height = 620;
set(gcf,'units','points','position',[x0,y0,width,height]);
saveas(gca,[savepath filesep 'fig_para.png']);