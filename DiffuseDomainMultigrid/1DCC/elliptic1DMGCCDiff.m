%
%   DDMSolver1D.m
%
%   Solves the 1D Diffuse Domain Approximation problem corresponding to 
%   Figure 4.1 in the PDF.
%
%   Equation:
%   -(D u')' + C u = F
%
%   Where:
%   D = alpha + (1-alpha)phi
%   C = c_eps + kappa * phi'
%   F = f_eps - lambda * phi'
%
%   Naming Convention: CamelCase
%   Solver: Multigrid (Cell-Centered)
%

clear; clc;

%% 1. Simulation Parameters

% Use a power of 2 for Multigrid compatibility
% 2^18 = 262144
numLevels = 18; 
numCells  = 2^numLevels;

% Domain
xLower = -1.0;
xUpper =  1.0;
hMesh  = (xUpper - xLower) / numCells;

% Problem Parameters (Matching PDF Section 4.1 / Figure 4.1)
epsilon = 0.0025; 
alpha   = 0.5;     
beta    = 1.0;
gamma   = 1.0;
kappa   = 1.6; 

%% 2. Grid and Coordinate Setup

% Cell Centers
xCenter = linspace(xLower + hMesh/2, xUpper - hMesh/2, numCells)';

% Cell Edges
xEdge = linspace(xLower, xUpper, numCells+1)';

%% 3. Construct Sharp Interface Solution (uRef)
% Used for calculating the Physical Error ||u_eps - u_0||

uRefCenter = zeros(numCells, 1);

for i = 1:numCells
    xVal = xCenter(i);
    if xVal <= 0
        % Left domain: u_L(x) = 8(x+1)^2 - 2
        uRefCenter(i) = 8 * (xVal + 1)^2 - 2;
    else
        % Right domain: u_R(x) = (4x^2 - 8x + 6)cos(4pi x)
        uRefCenter(i) = (4 * xVal^2 - 8 * xVal + 6) * cos(4 * pi * xVal);
    end
end

% Fill ghost layers for BCs (Neumann) in uRef
uRef = zeros(numCells + 2, 1);
uRef(2:numCells + 1) = uRefCenter;
uRef = applyBCs(uRef);

%% 4. Calculate Lambda Explicitly
% Condition: u'_R(0) - alpha*u'_L(0) = kappa*u(0) + lambda

uZero      = 6.0;
duLeftZero = 16.0;   % u_L'(0)
duRightZero= -8.0;   % u_R'(0)

lambda = (duRightZero - alpha * duLeftZero) - kappa * uZero;

fprintf('Calculated Lambda: %.4f (Expected -25.6)\n', lambda);

%% 5. Pre-Compute Coefficient Arrays (No Inline Functions)

% --- Phase Field Computations ---
% phi(x) = 0.5 * (1 + tanh(3x/eps))
% phi'(x) = (1.5 * sech(3x/eps)^2) / eps

% Compute at Edges (for Diffusion Coefficient)
argEdge = 3.0 * xEdge / epsilon;
phiEdge = 0.5 * (1.0 + tanh(argEdge));

% Compute at Centers (for Reaction and RHS)
argCenter = 3.0 * xCenter / epsilon;
phiCenter = 0.5 * (1.0 + tanh(argCenter));
phiPrimeCenter = (1.5 * sech(argCenter).^2) / epsilon;

% --- Diffusion Coefficient D (at Edges) ---
% D = alpha + (1-alpha)phi
diffCoefEdge = alpha + (1.0 - alpha) * phiEdge;

% --- Reaction Coefficient C (at Centers) ---
% c_eps = beta(1-phi) + gamma*phi
cEps = beta * (1.0 - phiCenter) + gamma * phiCenter;

% Operator C = c_eps + kappa * phi'
reacCoefCenter = cEps + kappa * phiPrimeCenter;

% --- Forcing Functions h(x) and q(x) (at Centers) ---

% h(x)
hFunc = -alpha * 16.0 + beta * (8.0 * (xCenter + 1.0).^2 - 2.0);

% q(x)
term1 = (1.0 - 4.0 * pi^2 * (2.0 * xCenter.^2 - 4.0 * xCenter + 3.0)) .* cos(4.0 * pi * xCenter);
term2 = 8.0 * pi * (xCenter - 1.0) .* sin(4.0 * pi * xCenter);
term3 = gamma * (4.0 * xCenter.^2 - 8.0 * xCenter + 6.0) .* cos(4.0 * pi * xCenter);
qFunc = -8.0 * (term1 - term2) + term3;

% f_eps = h(1-phi) + q*phi
fEps = hFunc .* (1.0 - phiCenter) + qFunc .* phiCenter;

% Total RHS Vector F = f_eps - lambda * phi'
rhsVector = fEps - lambda * phiPrimeCenter;

%% 6. Multigrid Configuration

MGParam.nL     = numCells;
MGParam.xLower = xLower;
MGParam.xUpper = xUpper;
MGParam.L      = numLevels;
MGParam.pCycle = 1;       % V-cycle
MGParam.m1     = 2;       % Pre-smoothing
MGParam.m2     = 2;       % Post-smoothing
MGParam.omega  = 2/3;     % Damped Jacobi weight
MGParam.kMax   = 50;      % Max iterations
MGParam.tol    = 1.0e-11; 
MGParam.C      = reacCoefCenter; % Pass computed C array

uInit = zeros(numCells + 2, 1);

%% 7. Run Solver

fprintf('Starting Multigrid Solver with N = %d, epsilon = %.5f...\n', numCells, epsilon);
tic;
% Passing uRef as the "Exact" solution for error tracking column 3
[uSol, errorHistory, iterStop] = multiGridSolver(uInit, rhsVector, diffCoefEdge, ...
                                                 reacCoefCenter, hMesh, ...
                                                 MGParam, uRef);
solveTime = toc;

%% 8. Analysis and Plotting

uInterior = uSol(2:numCells + 1);

% Calculate L2 Difference against Sharp Interface Limit
diffVec = uInterior - uRefCenter;
l2Diff = sqrt(sum(diffVec.^2) / numCells);

fprintf('\nSolver finished in %.4f seconds.\n', solveTime);
fprintf('Final Residual Norm:   %.6e\n', errorHistory(iterStop, 1));
fprintf('Final Correction Norm: %.6e\n', errorHistory(iterStop, 2));
fprintf('L2 Diff (u_eps vs u_0): %.6e\n', l2Diff);

% --- Log-Linear Fit for Convergence Rate ---
% Using the Correction Norm (Column 2)
if iterStop >= 4
    startIndex = iterStop - 3;
    endIndex   = iterStop;
    
    yData = log(errorHistory(startIndex:endIndex, 2));
    xData = (startIndex:endIndex)';
    
    % Polyfit degree 1
    pFit = polyfit(xData, yData, 1);
    convergenceRate = exp(pFit(1));
    
    fprintf('Estimated Convergence Rate (rho): %.4f\n', convergenceRate);
else
    fprintf('Insufficient iterations for convergence fit.\n');
end

% --- Plotting ---

% Figure 1: Solution Comparison
figure(1)
clf
plot(xCenter, uRefCenter, 'r-', 'LineWidth', 3); hold on;
plot(xCenter, uInterior, 'b--', 'LineWidth', 2);
legend('Sharp Interface (u_0)', ['DDM (u_\epsilon), \epsilon=' num2str(epsilon)], 'Location', 'NorthWest');
title(['DDM Solution vs Sharp Interface (\alpha=' num2str(alpha) ')']);
xlabel('x'); ylabel('u');
grid on;
xlim([-1 1]);
ylim([-20 10]); 

% Figure 2: Convergence History
figure(2)
clf
semilogy(1:iterStop, errorHistory(1:iterStop, 1), 'k-o', 'LineWidth', 1.2, 'MarkerSize', 4); hold on;
semilogy(1:iterStop, errorHistory(1:iterStop, 2), 'r-s', 'LineWidth', 1.2, 'MarkerSize', 4); 
semilogy(1:iterStop, errorHistory(1:iterStop, 3), 'b--', 'LineWidth', 1.5);

% Add fit line
if iterStop >= 4
    xFit = 1:iterStop;
    yFit = exp(polyval(pFit, xFit));
    semilogy(xFit, yFit, 'g-', 'LineWidth', 1.0);
    legend('Algebraic Residual', 'Correction Norm', 'Physical Diff ||u_\epsilon - u_0||', ...
           ['Fit (\gamma_{comp} \approx ' num2str(convergenceRate, '%.2f') ')']);
else
    legend('Algebraic Residual', 'Correction Norm', 'Physical Diff ||u_\epsilon - u_0||');
end

title('Multigrid Convergence History');
xlabel('Cycle'); ylabel('Error / Residual');
grid on;