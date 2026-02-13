%
%   runConvergenceAnalysis.m
%
%   Performs convergence analysis of the 2D DDM solver.
%   Includes SMART MEMORY MANAGEMENT to prevent computer freezing.
%

clear; clc;

%% 1. Configuration

% --- Analysis Mode ---
% 1: Fixed Alpha, Vary Epsilon
% 2: Fixed Epsilon, Vary Alpha
% 3: Coupled Alpha = Epsilon^p
mode = 3; 

% --- Grid Settings ---
L = 13; 
nL = [2^L, 2^L]; 
xyLower = [-2.0, -2.0];
xyUpper = [ 2.0,  2.0];

% --- Multigrid Settings ---
MGParam.L      = L;
MGParam.xLower = xyLower;
MGParam.xUpper = xyUpper;
MGParam.pCycle = 2;
MGParam.m1     = 2;
MGParam.m2     = 2;
MGParam.omega  = 0.66;
MGParam.kMax   = 200;
MGParam.tol    = 1.0e-12;

% --- Parameter Ranges ---
numPoints = 20; 

% Default setups
epsRange_1   = logspace(log10(0.01), log10(0.0006), numPoints);
alphaFixed_1 = 3.0;

alphaRange_2 = logspace(log10(1.5), log10(0.008), numPoints);
epsFixed_2   = 0.05;

epsRange_3   = logspace(log10(0.01), log10(0.0006), numPoints);
p_exponent   = 2.0; 

% Fixed Physics Constants
beta = 0.0; gamma = 1.0; kappa = 2.01;

%% 2. Pre-computation

fprintf('Setting up Geometry and Pre-computing Fields (N=%d)...\n', nL(1));
hMesh = (xyUpper(1) - xyLower(1)) / nL(1);

xCenter = linspace(xyLower(1) + hMesh/2, xyUpper(1) - hMesh/2, nL(1));
yCenter = linspace(xyLower(2) + hMesh/2, xyUpper(2) - hMesh/2, nL(2));
[XCC, YCC] = ndgrid(xCenter, yCenter);

xEdgePts = linspace(xyLower(1), xyUpper(1), nL(1)+1);
yEdgePts = linspace(xyLower(2), xyUpper(2), nL(2)+1);
[XEW, YEW] = ndgrid(xEdgePts, yCenter);
[XNS, YNS] = ndgrid(xCenter, yEdgePts);

% SDF Caching
sdfFileName = sprintf('SDF_Data_N%d.mat', nL(1));
if exist(sdfFileName, 'file')
    fprintf('  Loading pre-computed SDF...\n');
    load(sdfFileName, 'distCenter', 'distEW', 'distNS', 'polyX', 'polyY');
else
    fprintf('  Computing SDF via C++ Fast Marching...\n');
    [polyX, polyY] = generateBoundary(2500);
    distCenter = computeSignedDistance(XCC, YCC, polyX, polyY, hMesh);
    distEW = computeSignedDistance(XEW, YEW, polyX, polyY, hMesh);
    distNS = computeSignedDistance(XNS, YNS, polyX, polyY, hMesh);
    save(sdfFileName, 'distCenter', 'distEW', 'distNS', 'polyX', 'polyY');
end

% Pre-compute forcing
fprintf('  Pre-computing forcing fields...\n');
qFunc = -XCC.^2 +  15.0;
% hFunc = 2.5 .* sin(XCC) + exp(cos(YCC));
hFunc = 0.0 * ones(nL(1), nL(2));
gFunc = sin(XCC) + cos(XCC).^2 + 5;

% Pack Data
precompData.distCenter = distCenter;
precompData.distEW = distEW;
precompData.distNS = distNS;
precompData.qFunc = qFunc;
precompData.hFunc = hFunc;
precompData.gFunc = gFunc;

gridInfo.nL = nL;
gridInfo.hMesh = hMesh;

% Clear heavy variables to free Main RAM before forking
clear XCC YCC XEW YEW XNS YNS distCenter distEW distNS qFunc hFunc gFunc;

%% 3. Experiment Setup

paramList = repmat(struct('epsilon',0,'alpha',0,'beta',beta,'gamma',gamma,'kappa',kappa), numPoints, 1);
xPlot = zeros(numPoints, 1); 

switch mode
    case 1, xLabelStr = '$\varepsilon$';
        for i = 1:numPoints, paramList(i).epsilon = epsRange_1(i); paramList(i).alpha = alphaFixed_1; xPlot(i) = epsRange_1(i); end
    case 2, xLabelStr = '$\alpha$';
        for i = 1:numPoints, paramList(i).epsilon = epsFixed_2; paramList(i).alpha = alphaRange_2(i); xPlot(i) = alphaRange_2(i); end
    case 3, xLabelStr = '$\varepsilon$';
        for i = 1:numPoints, epsVal = epsRange_3(i); paramList(i).epsilon = epsVal; paramList(i).alpha = epsVal^p_exponent; xPlot(i) = epsVal; end
end

%% 4. Smart Parallel Pool Configuration (THE FIX)

% --- Memory Estimation ---
% A double is 8 bytes.
% MG Solver stores: u, rhs, diffEW, diffNS, reac (5 arrays).
% Plus internal MG levels (approx 1.33x overhead).
% Plus temporary smoothers/residual vectors (approx 3x overhead).
% Rough heuristic: ~30-40 vectors of size N*N.
bytesPerDouble = 8;
elements = nL(1) * nL(2);
memPerWorkerGB = (40 * elements * bytesPerDouble) / 1024^3; 

fprintf('\n--- Resource Management ---\n');
fprintf('Estimated RAM per Worker: %.2f GB\n', memPerWorkerGB);

% --- Detect System RAM (Mac/Linux/Windows) ---
try
    if ismac || isunix
        [~, memStr] = system('sysctl -n hw.memsize'); % Mac
        if isempty(memStr), [~, memStr] = system('grep MemTotal /proc/meminfo | awk ''{print $2}'''); end % Linux
        totalRAM_GB = str2double(strtrim(memStr)) / 1024^3;
        if isunix && ~ismac, totalRAM_GB = totalRAM_GB / 1024; end % Linux returns kB
    elseif ispc
        [~, sys] = memory;
        totalRAM_GB = sys.PhysicalMemory.Total / 1024^3;
    else
        totalRAM_GB = 16; % Fallback assumption
    end
catch
    totalRAM_GB = 16; % Fallback if detection fails
end

% --- Calculate Safe Limit ---
% Leave 4GB for OS + MATLAB Main Thread
availableRAM = max(0, totalRAM_GB - 4.0); 

% Calculate max workers that fit in RAM
maxWorkersRAM = floor(availableRAM / memPerWorkerGB);

% Cap at physical cores (usually 8 or 10 on M1/M2/M3) or 1 if RAM is tight
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
% parfor will now distribute jobs to the limited pool size.
% If you have 8 jobs but pool size is 2, it runs 2 at a time (safe!).
parfor i = 1:numPoints
    fprintf('  Job %d running... (Eps=%.2e)\n', i, paramList(i).epsilon);
    solutions{i} = computeDiffuseDomainSolution2D(gridInfo, paramList(i), precompData, MGParam);
end
totalTime = toc;
fprintf('All jobs finished in %.2f seconds.\n', totalTime);

%% 6. Convergence Analysis

cauchyErrors = zeros(numPoints-1, 1);
xMid = zeros(numPoints-1, 1);

fprintf('\nComputing L-inf Cauchy Differences:\n');
for i = 1:numPoints-1
    u1 = solutions{i};
    u2 = solutions{i+1};
    
    diff = u1 - u2;
    linf_diff = max(abs(diff(:))); % Changed from L2 (RMS) to L-infinity (Max Abs)
    
    cauchyErrors(i) = linf_diff;
    xMid(i) = xPlot(i+1); 
end

%% 7. Plotting

figure(1); clf;
set(gcf, 'Color', 'w');

% Main Log-Log Plot
loglog(xMid, cauchyErrors, 'ko-', 'LineWidth', 1.5, 'MarkerSize', 10, 'MarkerFaceColor', 'b');
grid on;
hold on;

slope = NaN;
if length(xMid) >= 2
    % Linear Fit on Log Data
    p = polyfit(log(xMid), log(cauchyErrors), 1);
    slope = p(1);
    yFit = exp(polyval(p, log(xMid)));
    
    % Plot Fit Line (Thick dashed red)
    loglog(xMid, yFit, 'r--', 'LineWidth', 2.5);
    
    % Legend in Top Left
    legend('Cauchy Diff $||u(\cdot;\varepsilon_{k+1},\alpha_{k+1})-u(\cdot;\varepsilon_{k},\alpha_{k})||_{L^\infty}$', sprintf('Fit (Slope = %.2f)', slope), ...
           'Location', 'NorthWest', 'Interpreter', 'latex', 'FontSize', 12);
end

% Labels
xlabel(xLabelStr, 'Interpreter', 'latex', 'FontSize', 16);
ylabel('$L^\infty$ Difference', 'Interpreter', 'latex', 'FontSize', 16);

% Title
titleStr = sprintf('Convergence (Mode %d)', mode);
if mode == 3, titleStr = [titleStr, sprintf(', $p=%.1f$', p_exponent)]; end
title(titleStr, 'Interpreter', 'latex', 'FontSize', 16);

% Axis Ticks formatting
set(gca, 'FontSize', 12, 'TickLabelInterpreter', 'latex');

% Annotation Box in Bottom Right
% [x y w h] normalized coordinates (0,0 is bottom left)
dim = [0.65 0.15 0.2 0.1]; 
str = {sprintf('Points: %d', numPoints), ...
       sprintf('Slope: %.2f', slope)};

annotation('textbox', dim, 'String', str, ...
           'FitBoxToText', 'on', ...
           'BackgroundColor', 'w', ...
           'EdgeColor', 'k', ...
           'FontSize', 12, ...
           'Interpreter', 'latex');

fprintf('\nEstimated Convergence Rate: %.4f\n', slope);