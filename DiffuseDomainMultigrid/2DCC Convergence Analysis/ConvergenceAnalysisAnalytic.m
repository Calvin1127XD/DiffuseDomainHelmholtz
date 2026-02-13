%
%   runConvergenceAnalysis.m
%
%   Performs convergence analysis of the 2D DDM solver against an 
%   ANALYTICAL EXACT SOLUTION defined by a polynomial level set.
%   Includes SMART MEMORY MANAGEMENT to prevent uuter freezing.
%
%   Problem:
%       - div(grad u) + gamma*u = q   in Omega_1 (Interior)
%       - alpha * div(grad u)   = 0   in Omega_2 (Exterior)
%       Boundary Cond: u = 2 on interface.
%
%   Domain: [-3, 3]^2 (Enlarged to fit the polynomial tip at r~2.16)
%

clear; clc;

%% 1. Configuration

% --- Analysis Mode ---
% 1: Fixed Alpha, Vary Epsilon
% 2: Fixed Epsilon, Vary Alpha
% 3: Coupled Alpha = Epsilon^p (Recommended for Analytic Test)
mode = 3; 

% --- Grid Settings ---
L = 13; % N=512
nL = [2^L, 2^L]; 
xyLower = [-3.0, -3.0];
xyUpper = [ 3.0,  3.0];

% --- Multigrid Settings ---
MGParam.L      = L;
MGParam.xLower = xyLower;
MGParam.xUpper = xyUpper;
MGParam.pCycle = 2;
MGParam.m1     = 2;
MGParam.m2     = 2;
MGParam.omega  = 0.66;
MGParam.kMax   = 100;
MGParam.tol    = 1.0e-12;

% --- Parameter Ranges ---
numPoints = 20; 

% Default setups
epsRange_1   = logspace(log10(0.01), log10(0.006), numPoints);
alphaFixed_1 = 1.0;

alphaRange_2 = logspace(log10(1.5), log10(0.008), numPoints);
epsFixed_2   = 0.05;

epsRange_3   = logspace(log10(0.01), log10(0.002), numPoints);
p_exponent   = 1.0;

% --- Fixed Physics Constants for the Analytic Test ---
% These must match the derived exact solution.
beta = 0.0; 
gamma = 2.0; 
kappa = 1.0;

%% 2. Pre-computation

fprintf('Setting up Geometry and Exact Solution (N=%d)...\n', nL(1));
hMesh = (xyUpper(1) - xyLower(1)) / nL(1);

% Grid Generation
xCenter = linspace(xyLower(1) + hMesh/2, xyUpper(1) - hMesh/2, nL(1));
yCenter = linspace(xyLower(2) + hMesh/2, xyUpper(2) - hMesh/2, nL(2));
[XCC, YCC] = ndgrid(xCenter, yCenter);

xEdgePts = linspace(xyLower(1), xyUpper(1), nL(1)+1);
yEdgePts = linspace(xyLower(2), xyUpper(2), nL(2)+1);
[XEW, YEW] = ndgrid(xEdgePts, yCenter);
[XNS, YNS] = ndgrid(xCenter, yEdgePts);

% SDF Caching
sdfFileName = sprintf('SDF_PolyStar_Large_N%d.mat', nL(1));
if exist(sdfFileName, 'file')
    fprintf('  Loading pre-computed SDF...\n');
    load(sdfFileName, 'distCenter', 'distEW', 'distNS', 'polyX', 'polyY');
else
    fprintf('  Computing SDF via C++ Fast Marching...\n');
    [polyX, polyY] = generateBoundary(4000); % Uses generateBoundary.m
    distCenter = computeSignedDistance(XCC, YCC, polyX, polyY, hMesh);
    distEW = computeSignedDistance(XEW, YEW, polyX, polyY, hMesh);
    distNS = computeSignedDistance(XNS, YNS, polyX, polyY, hMesh);
    save(sdfFileName, 'distCenter', 'distEW', 'distNS', 'polyX', 'polyY');
end

% --- Construct Analytic Forcing and Reference Solution ---
fprintf('  Constructing Analytic Solution and Forcing terms...\n');

% Polynomial P(x,y) = (x^2+y^2)^2 + 2(x^3 - 3xy^2) - 2.5
R2 = XCC.^2 + YCC.^2;
HarmonicPart = 2.0 .* (XCC.^3 - 3.0 .* XCC .* YCC.^2);
P_val = R2.^2 + HarmonicPart - 2.5;

% q(x) = 16(x^2+y^2) + gamma*u1
% Note: u1 = 2 - P. So q = 16*r^2 + gamma*(2-P).
qFunc = 16.0 .* R2 + gamma .* (2.0 - P_val);

% g(x) = |grad P| - 2*kappa
Px = 4.0 .* XCC .* R2 + 6.0 .* (XCC.^2 - YCC.^2);
Py = 4.0 .* YCC .* R2 - 12.0 .* XCC .* YCC;
GradP_mag = sqrt(Px.^2 + Py.^2);
gFunc = GradP_mag - 2.0 * kappa;

% h(x) = 0
hFunc = zeros(size(XCC));

% Pack Data
precompData.distCenter = distCenter;
precompData.distEW = distEW;
precompData.distNS = distNS;
precompData.qFunc = qFunc;
precompData.hFunc = hFunc;
precompData.gFunc = gFunc;

gridInfo.nL = nL;
gridInfo.hMesh = hMesh;

% Clear heavy variables
clear XCC YCC XEW YEW XNS YNS distCenter distEW distNS qFunc hFunc gFunc P_val R2 Px Py;

%% 3. Experiment Setup

paramList = repmat(struct('epsilon',0,'alpha',0,'beta',beta,'gamma',gamma,'kappa',kappa), numPoints, 1);
xPlot = zeros(numPoints, 1); 

switch mode
    case 1 % Fixed Alpha, Vary Epsilon
        xLabelStr = '$\varepsilon$';
        varName = 'epsilon';
        for i = 1:numPoints
            paramList(i).epsilon = epsRange_1(i); 
            paramList(i).alpha = alphaFixed_1; 
            xPlot(i) = epsRange_1(i); 
        end
        
    case 2 % Fixed Epsilon, Vary Alpha
        xLabelStr = '$\alpha$';
        varName = 'alpha';
        for i = 1:numPoints
            paramList(i).epsilon = epsFixed_2; 
            paramList(i).alpha = alphaRange_2(i); 
            xPlot(i) = alphaRange_2(i); 
        end
        
    case 3 % Coupled
        xLabelStr = '$\varepsilon$';
        varName = 'epsilon';
        for i = 1:numPoints
            epsVal = epsRange_3(i); 
            paramList(i).epsilon = epsVal; 
            paramList(i).alpha = 0.1*epsVal^p_exponent; 
            xPlot(i) = epsVal; 
        end
end

%% 4. Smart Parallel Pool Configuration

% --- Memory Estimation ---
bytesPerDouble = 8;
elements = nL(1) * nL(2);
memPerWorkerGB = (45 * elements * bytesPerDouble) / 1024^3; 

fprintf('\n--- Resource Management ---\n');
fprintf('Estimated RAM per Worker: %.2f GB\n', memPerWorkerGB);

% --- Detect System RAM ---
try
    if ismac || isunix
        [~, memStr] = system('sysctl -n hw.memsize'); 
        if isempty(memStr), [~, memStr] = system('grep MemTotal /proc/meminfo | awk ''{print $2}'''); end 
        totalRAM_GB = str2double(strtrim(memStr)) / 1024^3;
        if isunix && ~ismac, totalRAM_GB = totalRAM_GB / 1024; end 
    elseif ispc
        [~, sys] = memory;
        totalRAM_GB = sys.PhysicalMemory.Total / 1024^3;
    else
        totalRAM_GB = 16; 
    end
catch
    totalRAM_GB = 16; 
end

% --- Calculate Safe Limit ---
availableRAM = max(0, totalRAM_GB - 4.0); 
maxWorkersRAM = floor(availableRAM / memPerWorkerGB);
maxWorkers = max(1, min(feature('numcores'), maxWorkersRAM));

fprintf('Total System RAM:       %.2f GB\n', totalRAM_GB);
fprintf('Safe RAM for Workers:   %.2f GB\n', availableRAM);
fprintf('Max Workers (RAM-bound): %d\n', maxWorkersRAM);
fprintf('<strong>Configuring Pool size to: %d</strong>\n', maxWorkers);

% --- Reconfigure Pool ---
currentPool = gcp('nocreate');
if ~isempty(currentPool)
    if currentPool.NumWorkers ~= maxWorkers
        delete(currentPool);
        parpool(maxWorkers);
    end
else
    parpool(maxWorkers);
end

%% 5. Parallel Solver Execution

fprintf('\nStarting Parallel Execution...\n');
solutions = cell(numPoints, 1);

tic;
parfor i = 1:numPoints
    fprintf('  Job %d running... (Eps=%.4f, Alpha=%.1e)\n', i, paramList(i).epsilon, paramList(i).alpha);
    solutions{i} = computeDiffuseDomainSolution2D(gridInfo, paramList(i), precompData, MGParam);
end
totalTime = toc;
fprintf('All jobs finished in %.2f seconds.\n', totalTime);

%% 6. Convergence Analysis (Against Exact Solution)

exactErrorsL2 = zeros(numPoints, 1);
exactErrorsInf = zeros(numPoints, 1);

fprintf('\nComputing Exact Errors:\n');
for i = 1:numPoints
    uNum = solutions{i};
    
    % Reconstruct uExact locally
    xC_ = linspace(xyLower(1)+hMesh/2, xyUpper(1)-hMesh/2, nL(1));
    yC_ = linspace(xyLower(2)+hMesh/2, xyUpper(2)-hMesh/2, nL(2));
    [XC_, YC_] = ndgrid(xC_, yC_);
    
    % Polynomial P
    R2_ = XC_.^2 + YC_.^2;
    Harm_ = 2.0 .* (XC_.^3 - 3.0 .* XC_ .* YC_.^2);
    P_val_ = R2_.^2 + Harm_ - 2.5;
    
    % Exact Solution: u = 2 - P inside, 2 outside.
    uRef = 2.0 * ones(size(XC_));
    maskIn = (P_val_ < 0); 
    uRef(maskIn) = 2.0 - P_val_(maskIn);
    
    diff = abs(uNum - uRef);
    
    exactErrorsInf(i) = max(diff(:));
    
    fprintf('  Point %d: L_inf Error = %.4e\n', i, exactErrorsInf(i));
end

%% 7. Plotting Convergence

figure(1); clf;
set(gcf, 'Color', 'w');

% Main Log-Log Plot
loglog(xPlot, exactErrorsInf, 'ko-', 'LineWidth', 1.5, 'MarkerSize', 10, 'MarkerFaceColor', 'b');
grid on;
hold on;

slope = NaN;
if length(xPlot) >= 2
    % Linear Fit on Log Data
    p = polyfit(log(xPlot), log(exactErrorsInf), 1);
    slope = p(1);
    yFit = exp(polyval(p, log(xPlot)));
    
    % Plot Fit Line
    loglog(xPlot, yFit, 'r--', 'LineWidth', 2.5);
    
    % Legend
    legend('Exact Error $||u_{num} - u_{exact}||_{L^\infty}$', sprintf('Fit (Slope = %.2f)', slope), ...
           'Location', 'Best', 'Interpreter', 'latex', 'FontSize', 12);
end

% Labels
xlabel(xLabelStr, 'Interpreter', 'latex', 'FontSize', 16);
ylabel('$L^\infty$ Error', 'Interpreter', 'latex', 'FontSize', 16);

% Title
titleStr = sprintf('Convergence Analysis (Mode %d)', mode);
if mode == 3, titleStr = [titleStr, sprintf(', $p=%.1f$', p_exponent)]; end
title(titleStr, 'Interpreter', 'latex', 'FontSize', 16);

% Axis Ticks formatting
set(gca, 'FontSize', 12, 'TickLabelInterpreter', 'latex');

% Annotation
dim = [0.15 0.15 0.2 0.1]; 
str = {sprintf('Points: %d', numPoints), ...
       sprintf('Slope: %.2f', slope)};

annotation('textbox', dim, 'String', str, ...
           'FitBoxToText', 'on', ...
           'BackgroundColor', 'w', ...
           'EdgeColor', 'k', ...
           'FontSize', 12, ...
           'Interpreter', 'latex');

fprintf('\nEstimated Convergence Rate: %.4f\n', slope);

%% 8. Comparative Visualization (4 Subplots)

figure(3); clf;
set(gcf, 'Color', 'w', 'Name', 'Comparison');

% Reconstruct Exact Solution Grid for Plotting
xC_plot = linspace(xyLower(1)+hMesh/2, xyUpper(1)-hMesh/2, nL(1));
yC_plot = linspace(xyLower(2)+hMesh/2, xyUpper(2)-hMesh/2, nL(2));
[XC_plot, YC_plot] = ndgrid(xC_plot, yC_plot);
R2_plot = XC_plot.^2 + YC_plot.^2;
Harm_plot = 2.0 .* (XC_plot.^3 - 3.0 .* XC_plot .* YC_plot.^2);
P_val_plot = R2_plot.^2 + Harm_plot - 2.5;

uExactPlot = 2.0 * ones(size(XC_plot));
maskInPlot = (P_val_plot < 0);
uExactPlot(maskInPlot) = 2.0 - P_val_plot(maskInPlot);

% Identify Indices
idxLarge = 1;                   
idxMid   = round(numPoints / 2);
idxSmall = numPoints;           

% Setup Titles based on variable
if strcmp(varName, 'epsilon')
    lblLarge = sprintf('Largest $\\varepsilon = %.5f$', paramList(idxLarge).epsilon);
    lblMid   = sprintf('Medium $\\varepsilon = %.5f$', paramList(idxMid).epsilon);
    lblSmall = sprintf('Smallest $\\varepsilon = %.5f$', paramList(idxSmall).epsilon);
else
    lblLarge = sprintf('Largest $\\alpha = %.5f$', paramList(idxLarge).alpha);
    lblMid   = sprintf('Medium $\\alpha = %.5f$', paramList(idxMid).alpha);
    lblSmall = sprintf('Smallest $\\alpha = %.5f$', paramList(idxSmall).alpha);
end

% --- Subplot 1: Exact Solution ---
subplot(2,2,1);
surf(XC_plot, YC_plot, uExactPlot);
shading interp; lighting phong; camlight; axis tight;
title('Exact Sharp Interface Solution', 'Interpreter', 'latex', 'FontSize', 12);
xlabel('x'); ylabel('y'); zlabel('u');
view(2); colorbar;

% --- Subplot 2: Largest Parameter ---
subplot(2,2,2);
surf(XC_plot, YC_plot, solutions{idxLarge});
shading interp; lighting phong; camlight; axis tight;
title(['Numerical: ' lblLarge], 'Interpreter', 'latex', 'FontSize', 12);
xlabel('x'); ylabel('y'); zlabel('u');
view(2); colorbar;

% --- Subplot 3: Medium Parameter ---
subplot(2,2,3);
surf(XC_plot, YC_plot, solutions{idxMid});
shading interp; lighting phong; camlight; axis tight;
title(['Numerical: ' lblMid], 'Interpreter', 'latex', 'FontSize', 12);
xlabel('x'); ylabel('y'); zlabel('u');
view(2); colorbar;

% --- Subplot 4: Smallest Parameter ---
subplot(2,2,4);
surf(XC_plot, YC_plot, solutions{idxSmall});
shading interp; lighting phong; camlight; axis tight;
title(['Numerical: ' lblSmall], 'Interpreter', 'latex', 'FontSize', 12);
xlabel('x'); ylabel('y'); zlabel('u');
view(2); colorbar;