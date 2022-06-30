%% Dimensions
npara = 100;
nvar = 100;
nshock = 20;
ndata = 20;
period = 1000;

%% Input types
x = coder.typeof(0,[1 npara],[false true]);
mu = coder.typeof(0,[1 npara],[false true]);
Sigma = coder.typeof(0,[npara npara],[true true]);
v = coder.typeof(0);
N = coder.typeof(0);

stake = coder.typeof(0);
A = coder.typeof(1i,[nvar nvar],[true true]);
B = coder.typeof(1i,[nvar nvar],[true true]);
Q = coder.typeof(1i,[nvar nvar],[true true]);
Z = coder.typeof(1i,[nvar nvar],[true true]);

Y = coder.typeof(0,[period ndata],[true true]);
SSR.Sigma_u = coder.typeof(0,[ndata 1],[true false]);
SSR.Z = coder.typeof(0,[ndata nvar],[true true]);
SSR.D = coder.typeof(0,[ndata 1],[true false]);
SSR.G = coder.typeof(0,[nvar nvar],[true true]);
SSR.C = coder.typeof(0,[nvar 1],[true false]);
SSR.M = coder.typeof(0,[nvar nshock],[true true]);
SSR.sdof = coder.typeof(0);
SSR.sv = coder.typeof(true);
SSR.Sigma_v = coder.typeof(0,[nshock 1],[true false]);
SSR.B = coder.typeof(0,[nshock nshock],[true true]);
SSR.A = coder.typeof(0,[nshock 1],[true false]);
SSR.Sigma_e = coder.typeof(0,[nshock period],[true true]);
fs_init = coder.typeof(0,[nvar 1],[true false]);
Omega_fs_init = coder.typeof(0,[nvar nvar],[true true]);
T = coder.typeof(0);
state = coder.typeof(0,[nvar 1],[true false]);
log_vol = coder.typeof(0,[nshock 1],[true false]);

Y2 = coder.typeof(0,[period 1],[true false]);
for k = 1:10
    SSR2.C{k} = coder.typeof(0);
    SSR2.G{k} = coder.typeof(0);
    SSR2.Sigma_e{k} = coder.typeof(0);
    SSR2.D{k} = coder.typeof(0);
    SSR2.Z{k} = coder.typeof(0);
    SSR2.Sigma_u{k} = coder.typeof(0);
end
P = coder.typeof(0,[10 10],[false false]);
R = coder.typeof(0,[1 period],[false true]);
draw_S = coder.typeof(true);
draw_R = coder.typeof(true);

% %% mvt_rnd.m
% fprintf('Generating mvt_rnd_mex... ');
% codegen -d mex mvt_rnd -args {mu,Sigma,v,N}
% movefile mvt_rnd_mex* mex
% fprintf('Done!\n');
% 
% %% mvt_pdf.m
% fprintf('Generating mvt_pdf_mex... ');
% codegen -d mex mvt_pdf -args {x,mu,Sigma,v}
% movefile mvt_pdf_mex* mex
% fprintf('Done!\n');
% 
% %% qzdiv.m
% fprintf('Generating qzdiv_mex... ');
% codegen -d mex qzdiv -args {stake,A,B,Q,Z}
% movefile qzdiv_mex* mex
% fprintf('Done!\n');
% 
% %% MixKalman.m
% fprintf('Generating MixKalman_mex... ');
% codegen -d mex MixKalman -args {Y,SSR,fs_init,Omega_fs_init,N}
% movefile MixKalman_mex* mex
% fprintf('Done!\n');

%% DistSmoother.m
fprintf('Generating DistSmoother_mex... ');
codegen -d mex DistSmoother -args {Y,SSR,fs_init,Omega_fs_init}
movefile DistSmoother_mex* mex
fprintf('Done!\n');

%% GibbsSmoother.m
fprintf('Generating GibbsSmoother_mex... ');
codegen -d mex GibbsSmoother -args {Y2,SSR2,P,R,draw_S,draw_R}
movefile GibbsSmoother_mex* mex
fprintf('Done!\n');

%% SimuData.m
fprintf('Generating SimuData_mex... ');
codegen -d mex SimuData -args {SSR,T,state,log_vol}
movefile SimuData_mex* mex
fprintf('Done!\n');
rmdir(['mex' filesep 'interface'],'s')