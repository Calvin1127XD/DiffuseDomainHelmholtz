%
%   Elliptic2DMGCCDiff.m
%
%   Solves a 2D Diffuse Domain Method (DDM) approximation for a two-sided
%   problem with transmission-type boundary conditions on a complex interface.
%
%   -----------------------------------------------------------------------
%   1. GENERIC SOLVER FORMULATION
%   -----------------------------------------------------------------------
%   The underlying Multigrid solver targets the generic linear elliptic PDE:
%
%       - div( D(x) * grad(u) ) + C(x) * u = F(x)    in Omega
%
%   -----------------------------------------------------------------------
%   2. DDM FORMULATION
%   -----------------------------------------------------------------------
%   The DDM approximates the interface problem by introducing a phase-field
%   function phi and reformulating the equation on a regular domain:
%
%       - div( D_eps * grad(u) ) + c_eps * u + (kappa * u + g) * |grad(phi)| = f_eps
%
%   Where:
%       D_eps = alpha + (1 - alpha) * phi
%       c_eps = beta * (1 - phi) + gamma * phi
%       f_eps = h * (1 - phi) + q * phi
%       |grad(phi)| approx (1 / 2*eps) * sech^2(dist / eps)
%
%   -----------------------------------------------------------------------
%   3. ALGEBRAIC REARRANGEMENT
%   -----------------------------------------------------------------------
%   To map the DDM equation into the generic solver form, the boundary
%   source terms associated with the interface delta function (|grad(phi)|)
%   are distributed as follows:
%
%   A. The linear term (kappa * u * |grad(phi)|) is absorbed into the
%      reaction coefficient C(x).
%   B. The source term (g * |grad(phi)|) is moved to the RHS forcing F(x).
%
%   Effective Solver Coefficients:
%       D(x) = D_eps
%       C(x) = c_eps + kappa * |grad(phi)|
%       F(x) = f_eps - g * |grad(phi)|
%

clear; clc;

%% 1. Simulation Parameters

% Grid Resolution
L = 10; 
nL = [2^L, 2^L]; 

% Box Domain
xyLower = [-2.0, -2.0];
xyUpper = [ 2.0,  2.0];

hMesh = (xyUpper(1) - xyLower(1)) / nL(1);

% Problem Parameters (PDF Eq 4.6)
epsilon = 0.05; 
alpha   = 3.0;
beta    = 2.0;
gamma   = 1.0;
kappa   = 0.01;

% Check if the level is smaller than or equal to 13
if L > 13
    error("The grid is too fine.")
end

% Multigrid Parameters
MGParam.nL     = nL;
MGParam.xLower = xyLower;
MGParam.xUpper = xyUpper;
MGParam.L      = L;
MGParam.pCycle = 1;      
MGParam.m1     = 3;      
MGParam.m2     = 3;      
MGParam.omega  = 0.66;    
MGParam.kMax   = 100;
MGParam.tol    = 1.0e-12;

%% 2. Grid Setup

xCenter = linspace(xyLower(1) + hMesh/2, xyUpper(1) - hMesh/2, nL(1));
yCenter = linspace(xyLower(2) + hMesh/2, xyUpper(2) - hMesh/2, nL(2));
[XCC, YCC] = ndgrid(xCenter, yCenter);

xEdgePts = linspace(xyLower(1), xyUpper(1), nL(1)+1);
yEdgePts = linspace(xyLower(2), xyUpper(2), nL(2)+1);

[XEW, YEW] = ndgrid(xEdgePts, yCenter);
[XNS, YNS] = ndgrid(xCenter, yEdgePts);

%% 3. Geometry and Distance Function (Interactive Loading)

% Changed filename to avoid conflict with Starfish data
sdfFileName = sprintf('SDF_Data_N%d.mat', nL(1));
dataLoaded = false;

if exist(sdfFileName, 'file')
    msg = sprintf('Found existing SDF data "%s". Load it? [y/n]: ', sdfFileName);
    userResp = input(msg, 's');
    if strcmpi(userResp, 'y')
        fprintf('Loading SDF data from file...\n');
        load(sdfFileName, 'distCenter', 'distEW', 'distNS', 'polyX', 'polyY');
        dataLoaded = true;
    end
end

if ~dataLoaded
    fprintf('Generating Boundary (Ellipse)...\n');
    [polyX, polyY] = generateBoundary(3000);
    
    fprintf('Computing Signed Distance Function (Fast Marching)...\n');
    distCenter = computeSignedDistance(XCC, YCC, polyX, polyY, hMesh);
    distEW = computeSignedDistance(XEW, YEW, polyX, polyY, hMesh);
    distNS = computeSignedDistance(XNS, YNS, polyX, polyY, hMesh);
    
    fprintf('Saving SDF data to "%s"...\n', sdfFileName);
    save(sdfFileName, 'distCenter', 'distEW', 'distNS', 'polyX', 'polyY');
end

%% 4. Construct Analytical Solution (u_0) for Error Calculation

fprintf('Constructing Analytical Solution u_0...\n');

% Equation (4.5):
% u1 = -x1^2 - 4x2^2 + 6   (Inside: x1^2 + 4x2^2 < 4)
% u2 = 2                   (Outside)

uRefInterior = zeros(nL(1), nL(2));

% Check ellipse condition: x^2 + 4y^2 < 4
% Note: distCenter > 0 is inside, but using coordinate check is exact.
isInside = (XCC.^2 + 4*YCC.^2) < 4;

uRefInterior(isInside)  = -XCC(isInside).^2 - 4*YCC(isInside).^2 + 6;
uRefInterior(~isInside) = 2.0;

% Fill uRef with ghost layers (Neumann BCs for consistency in plotting/norm)
uRef = zeros(nL(1)+2, nL(2)+2);
uRef(2:end-1, 2:end-1) = uRefInterior;
uRef = applyBCs(uRef);
hasRef = true;

%% 5. Pre-Compute Coefficients

fprintf('Computing DDM Coefficients...\n');

% --- Phase Field Arrays ---
argEW = distEW ./ epsilon;
phiEW = 0.5 .* (1.0 + tanh(argEW));
diffCoefEW = alpha + (1.0 - alpha) .* phiEW;

argNS = distNS ./ epsilon;
phiNS = 0.5 .* (1.0 + tanh(argNS));
diffCoefNS = alpha + (1.0 - alpha) .* phiNS;

argCC = distCenter ./ epsilon;
phiCC = 0.5 .* (1.0 + tanh(argCC));

% Gradient Magnitude
absGradPhi = (1.0 / (2.0 * epsilon)) .* (sech(argCC).^2);

% Reaction: C = c_eps + kappa * |grad phi|
cEps = beta .* (1.0 - phiCC) + gamma .* phiCC;
reacCoef = cEps + kappa .* absGradPhi;

% --- Forcing Terms (Eq 4.7) ---
% q(x) = -x1^2 - 4x2^2 + 16
qFunc = -XCC.^2 + 15.0;
% h(x) = 2
hFunc = 2.5*sin(XCC)+exp(cos(YCC));

fEps = hFunc .* (1.0 - phiCC) + qFunc .* phiCC;

% --- Boundary Forcing g(x) (Eq 4.8) ---
% g(x) = x1^2 + 4x2^2 - 6 + 2*sqrt(x1^2 + 16x2^2)
gFunc = XCC - XCC + 4;

% Total RHS = f_eps - g(x) * |grad phi|
rhsVector = fEps - gFunc .* absGradPhi;

%% 6. Run Multigrid Solver

fprintf('Starting Multigrid Solver (N=%dx%d)...\n', nL(1), nL(2));

uInit = zeros(nL(1)+2, nL(2)+2);

tic;
% We pass uRef so the solver calculates ||u_k - u_0|| at every step
[uSol, errVals, kStop] = multiGridSolver(uInit, rhsVector, diffCoefEW, ...
                                         diffCoefNS, reacCoef, hMesh, ...
                                         MGParam, uRef);
solveTime = toc;

fprintf('Solver converged in %d iterations (%.4f seconds).\n', kStop, solveTime);

% Calculate final L2 error against analytical u0
uInterior = uSol(2:end-1, 2:end-1);
diff = uInterior - uRefInterior;
finalL2Err = sqrt(sum(diff(:).^2) / numel(diff));
fprintf('Final L2 Error ||u_eps - u_0|| = %.4e\n', finalL2Err);

%% 7. Plotting

% --- Figure 1: 3D Solution Surface ---
figure(1); clf;
set(gcf, 'Name', 'DDM Solution 3D', 'Color', 'w');

% Use surf for 3D visualization
surf(XCC, YCC, uInterior, 'EdgeColor', 'none', 'FaceColor', 'interp');
hold on;

% Interpolate solution onto the boundary curve to plot the interface in 3D
F_interp = griddedInterpolant(XCC, YCC, uInterior, 'linear', 'none');
u_on_boundary = F_interp(polyX, polyY);

% Plot boundary curve slightly lifted to prevent Z-fighting
plot3(polyX, polyY, u_on_boundary + 0.05, 'r-', 'LineWidth', 3);
hold off;

% Visual settings
axis tight;
view(-45, 30); 
colormap('jet');
colorbar;
camlight('headlight'); lighting phong; material dull;

title(['DDM Solution u_\epsilon (\epsilon = ' num2str(epsilon) ')'], 'FontSize', 14);
xlabel('x_1', 'FontSize', 12); 
ylabel('x_2', 'FontSize', 12);
zlabel('u', 'FontSize', 12);
set(gca, 'FontSize', 12);
rotate3d on;

% --- Figure 2: Convergence History ---
figure(2); clf;
set(gcf, 'Name', 'Convergence History', 'Color', 'w');

% Plot Norms
semilogy(1:kStop, errVals(1:kStop, 1), 'k-o', 'LineWidth', 1.2, 'MarkerSize', 6);
hold on;
semilogy(1:kStop, errVals(1:kStop, 2), 'r-s', 'LineWidth', 1.2, 'MarkerSize', 6);

legendLabels = {'Residual Norm $||r_k||$', 'Correction Norm $||u_k - u_{k-1}||$'};

if hasRef
    semilogy(1:kStop, errVals(1:kStop, 3), 'b-d', 'LineWidth', 1.2, 'MarkerSize', 6);
    legendLabels{end+1} = 'Analytic Error $||u_k - u_0||$';
end

% --- Log-Linear Fit for Correction Norm ---
gamma_comp = 0;
if kStop >= 4
    fitRange = kStop-3:kStop;
    yData = log(errVals(fitRange, 2));
    xData = fitRange';
    
    p = polyfit(xData, yData, 1);
    gamma_comp = exp(p(1));
    
    xFit = (kStop-4):kStop;
    yFit = exp(polyval(p, xFit));
    
    plot(xFit, yFit, 'g--', 'LineWidth', 2.0);
    legendLabels{end+1} = 'Log-Linear Fit';
    
    txtStr = sprintf('$\\gamma_{\\rm comp} = %.4f$', gamma_comp);
    text(kStop - 0.5, errVals(kStop, 2) * 2, txtStr, ...
         'Interpreter', 'latex', 'FontSize', 14, ...
         'HorizontalAlignment', 'right', 'VerticalAlignment', 'bottom', ...
         'BackgroundColor', 'w', 'EdgeColor', 'k');
end

grid on;
legend(legendLabels, 'Location', 'NorthEast', 'Interpreter', 'latex', 'FontSize', 12);
title('Multigrid Convergence History', 'FontSize', 14);
xlabel('Iteration $k$', 'Interpreter', 'latex', 'FontSize', 12); 
ylabel('Norm', 'Interpreter', 'latex', 'FontSize', 12);
set(gca, 'FontSize', 12);