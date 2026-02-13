%
%   elliptic1DMGCCDiff.m
%
% Calling script for solving the 1D elliptic equation
%
%   -(D u')' + C(x) u = f
%
% on an interval, where D is a possibly non-constant,
% everywhere positive diffusion coefficient, and C(x) is a 
% non-negative function. The boundary conditions may be 
% taken to be (i) homogeneous Dirichlet,
% (ii) homogeneous Neumann, or (iii) periodic. 
% See the function applyBCs.
% 
% Discretization is accomplished using the cell-centered 
% finite difference method. The cell-centered multigrid
% method is used as the solver.
%
% If C(x) is zero everywhere, then homogeneous Dirichlet
% boundary conditions must be chosen.
%
%% Initialization and Parameter

% Clear command line window and workspace.
clear; clc;

% Tolerance of residual to stop multigrid iterations.
tol = 1.0e-11;

% Number of interior cells at finest level.
nL = 1024;

% Number of multigrid levels.
L = 10;

% Define the domain of the equation.
xLower = 0.0;
xUpper = 1.0;

hL = (xUpper-xLower)/nL;

% Parameters for the multigrid solver.
pCycle = 1;
m1     = 3;
m2     = 3;
omega  = 1/2;
kMax   = 100;

%% Setting up the multigrid solver

% Check that nL is divisible by 2^L.
nl = nL;
for k = L-1:-1:0
  if any(mod(nl,2) ~= 0)
    error('nL not divisible by 2^L.  Coarsening error.');
  end
  nl = nl/2;
end

MGParam.nL     = nL;
MGParam.xLower = xLower;
MGParam.xUpper = xUpper;
MGParam.L      = L;
MGParam.pCycle = pCycle;
MGParam.m1     = m1;
MGParam.m2     = m2;
MGParam.omega  = omega;
MGParam.kMax   = kMax;
MGParam.tol    = tol;

% Allocate solution arrays (with ghost layers).
u      = zeros(nL+2,1);
uExact = zeros(nL+2,1);
f      = zeros(nL  ,1);
xCC    = zeros(nL  ,1);
DeW    = zeros(nL+1,1);
C      = zeros(nL  ,1);

% Build an exact solution and define f accordingly, so we
% can test the true algebraic error.
for i = 1:nL
  x = (i-0.5)*hL+xLower;
  xArg = 3.0*pi*x/xUpper;
% uExact(i+1) = exp(cos(xArg));
  uExact(i+1) = exp(sin(xArg))-1.0;
  xCC(i) = x;
end

% We next build the finest-grid diffusion array DeW:
for i = 1:nL+1
  xEdge = (i-1.0)*hL+xLower;
  xArg = 2.0*pi*xEdge/xUpper;

% Evaluate D at the east-west edge point:
% DeW(i) = 1.0+0.5*sin(xArg);
  DeW(i) = 1.0;
end

% Define the reaction coefficient C(x) at the cell centers.
C = 1.5 + sin(2.0*pi*xCC/xUpper);

% Store C in the MGParam struct
MGParam.C = C;

% Apply the BCs to uExact, then compute f by calling 
% FDOperator.
uExact = applyBCs(uExact);
f = FDOperator(uExact,DeW,C,hL);

%% Call the Multigrid solver.

tic;
[u,errVals,kStop] = multiGridSolver(u,f,DeW,C,hL,...
                                    MGParam,uExact);
toc;

%% Plotting the numerical solution and error plots.

uPlot(1:nL) = u(2:nL+1);

figure(1)
clf

plot(xCC,uPlot,'b-','LineWidth',1.5);
title('Numerical Solution u');
xlabel('x');
ylabel('u');

figure(2)
clf
semilogy(errVals(1:kStop,3),'bo','LineWidth',1.5)
hold on
semilogy(errVals(1:kStop,2),'rs','LineWidth',1.5)
hold on
semilogy(errVals(1:kStop,1),'kd','LineWidth',1.5)
hold on

rate = 1.0;
if kStop >= 4
  kv = kStop-3:kStop;
  le = log(errVals(kStop-3:kStop,3));
  p1 = polyfit(kv,le,1);
  rate = exp(p1(1));
  p1k = polyval(p1,1:kStop);
  semilogy(1:kStop,exp(p1k),'-k')
else
  semilogy(errVals(1:kStop,1),'.','LineWidth',1.5)
end

xlabel('$k$','Interpreter','latex');
title('Multigrid Iteration Errors',...
  'Interpreter','latex');
legend(...
  '$\left\|{\bf u}_L^{\rm E}-{\bf u}_L^k\right\|_L$',...
  '$\left\|{\bf u}_L^k-{\bf u}_L^{k-1} \right\|_L$',...
  '$\left\|{\bf r}_L^k\right\|_L$',...
  'log-linear fit','Interpreter','latex');
text(1.5,128*MGParam.tol,strcat(...
  '$\gamma_{\rm comp} =\hspace{.1cm}$',...
  num2str(rate,'%10.5e')),'FontSize',14,...
  'Interpreter','latex')
text(1.5,32*MGParam.tol,strcat('$m_1 =\hspace{.1cm}$',...
  num2str(MGParam.m1)),'FontSize',14,...
  'Interpreter','latex')
text(1.5,8*MGParam.tol,strcat('$m_2 =\hspace{.1cm}$',...
  num2str(MGParam.m2)),'FontSize',14,...
  'Interpreter','latex')
text(1.5,2*MGParam.tol,strcat('$p =\hspace{.1cm}$',...
  num2str(MGParam.pCycle)),'FontSize',14,...
  'Interpreter','latex')
text(1.5,0.5*MGParam.tol,...
  strcat('$\omega =\hspace{.1cm}$',...
  num2str(MGParam.omega)),'FontSize',14,...
  'Interpreter','latex')
printstr = strcat('Err_nL_',num2str(MGParam.nL(1)),...
  '_m1_',num2str(MGParam.m1),'.pdf');
exportgraphics(gca, printstr)
hold off