% Benchmark_FMM.m
% Compares C++ MEX Fast Marching vs Pure MATLAB Implementation

clear; clc;

% 1. Setup
N_values = [256, 512, 1024, 2048]; % Resolutions to test
numTests = length(N_values);

xyLower = [-3.0, -3.0];
xyUpper = [ 3.0,  3.0];

fprintf('==========================================================\n');
fprintf('  FMM BENCHMARK: C++ MEX vs MATLAB\n');
fprintf('==========================================================\n');
fprintf('|   N   |  Time C++ (s) | Time ML (s) | Speedup | L2 Diff |\n');
fprintf('|-------|---------------|-------------|---------|---------|\n');

for k = 1:numTests
    n = N_values(k);
    
    % Grid
    h = (xyUpper(1) - xyLower(1)) / n;
    x = linspace(xyLower(1)+h/2, xyUpper(1)-h/2, n);
    y = linspace(xyLower(2)+h/2, xyUpper(2)-h/2, n);
    [X, Y] = ndgrid(x, y);
    
    % Boundary
    [polyX, polyY] = generateBoundary(2000);
    
    % --- Run C++ MEX ---
    % Warmup (first run can handle memory allocation overhead)
    phi_cpp = mex_computeSignedDistance(X, Y, polyX, polyY, h); 
    
    tic;
    phi_cpp = mex_computeSignedDistance(X, Y, polyX, polyY, h);
    t_cpp = toc;
    
    % --- Run Pure MATLAB ---
    % Warning: large grids in MATLAB loop are slow
    tic;
    phi_ml = computeSignedDistance_MATLAB(X, Y, polyX, polyY, h);
    t_ml = toc;
    
    % --- Comparison ---
    % Compute L2 difference relative to grid size
    diff = phi_cpp - phi_ml;
    l2_err = sqrt(sum(diff(:).^2) / numel(diff));
    
    speedup = t_ml / t_cpp;
    
    fprintf('| %5d | %13.4f | %11.4f | %7.1fx | %7.2e |\n', ...
        n, t_cpp, t_ml, speedup, l2_err);
        
    % Store for plotting
    results(k).N = n;
    results(k).t_cpp = t_cpp;
    results(k).t_ml = t_ml;
end
fprintf('==========================================================\n');

% Validation Plot for the largest grid
figure(1); clf;
tiledlayout(1,3);

nexttile;
contourf(X, Y, phi_cpp, 20, 'LineColor', 'none');
axis equal tight; title('C++ Result'); colorbar;

nexttile;
contourf(X, Y, phi_ml, 20, 'LineColor', 'none');
axis equal tight; title('MATLAB Result'); colorbar;

nexttile;
% Plot difference
contourf(X, Y, abs(phi_cpp - phi_ml), 20, 'LineColor', 'none');
axis equal tight; title('Absolute Difference'); colorbar;
colormap(gca, 'jet');