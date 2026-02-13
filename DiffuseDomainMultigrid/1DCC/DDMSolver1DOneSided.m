%
%   DDMSolver1D_OneSided.m
%
%   Reproduces the "One-Sided Problem" experiment shown in Figure 5.1 
%   of the PDF.
%
%   Problem:
%   Left Domain (-1, 0): u_L = 6  (Constant)
%   Right Domain (0, 1): u_R = (4x^2 - 8x + 6)cos(4pi x)
%
%   Parameters from Eqs (5.2)-(5.4) and Figure 5.1 Caption:
%   alpha = 0.01, beta = 0, gamma = 1, kappa = 1.6
%   epsilon = 0.05 (as per Figure 5.1 caption)
%

clear; clc;

%% 1. Simulation Parameters

% N = 2^15 = 32768 (As described in the text for this experiment)
numLevels = 15; 
numCells  = 2^numLevels;

% Domain
xLower = -1.0;
xUpper =  1.0;
hMesh  = (xUpper - xLower) / numCells;

% Problem Parameters
epsilon = 0.0005;   % Matching Figure 5.1 caption
alpha   = 0.01;   % Small alpha for one-sided approximation
beta    = 0.0;    % No reaction in left domain
gamma   = 1.0;
kappa   = 1.6; 

%% 2. Grid and Coordinate Setup

xCenter = linspace(xLower + hMesh/2, xUpper - hMesh/2, numCells)';
xEdge   = linspace(xLower, xUpper, numCells+1)';

%% 3. Construct Sharp Interface Solution (uRef)
% Exact solution u_0 defined in Eq (5.1)

uRefCenter = zeros(numCells, 1);

for i = 1:numCells
    xVal = xCenter(i);
    if xVal <= 0
        % Left domain: u_L(x) = 6
        uRefCenter(i) = 6.0;
    else
        % Right domain: u_R(x) = (4x^2 - 8x + 6)cos(4pi x)
        uRefCenter(i) = (4 * xVal^2 - 8 * xVal + 6) * cos(4 * pi * xVal);
    end
end

% Fill ghost layers for Neumann BCs
uRef = zeros(numCells + 2, 1);
uRef(2:numCells + 1) = uRefCenter;
uRef = applyBCs(uRef);

%% 4. Calculate Lambda Explicitly
% The solver automatically derives lambda from the jump condition of uRef.
% Condition: u'_R(0) - alpha*u'_L(0) = kappa*u(0) + lambda

uZero       = 6.0;
duLeftZero  = 0.0;   % Derivative of constant 6 is 0
duRightZero = -8.0;  % Derivative of u_R at x=0

% lambda = (Jump in Flux) - kappa * u(0)
lambda = (duRightZero - alpha * duLeftZero) - kappa * uZero;
% Expected: (-8 - 0.01*0) - 1.6*6 = -8 - 9.6 = -17.6

fprintf('Calculated Lambda: %.4f (Expected -17.6)\n', lambda);

%% 5. Pre-Compute Coefficient Arrays

% --- Phase Field Computations ---
argEdge = 3.0 * xEdge / epsilon;
phiEdge = 0.5 * (1.0 + tanh(argEdge));

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

% h(x) = 0 (since beta=0 and u_L is constant)
hFunc = zeros(numCells, 1);

% q(x) defined in Eq (5.3)
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
MGParam.pCycle = 1;       
MGParam.m1     = 3;       
MGParam.m2     = 3;       
MGParam.omega  = 2/3;     
MGParam.kMax   = 100;      
MGParam.tol    = 1.0e-11; 
MGParam.C      = reacCoefCenter; 

uInit = zeros(numCells + 2, 1);

%% 7. Run Solver

fprintf('Starting Multigrid Solver with N = %d, epsilon = %.5f...\n', numCells, epsilon);
tic;
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

% --- Plotting Results (reproducing Figure 5.1) ---

% Figure 1: Solution Comparison
figure(1)
clf
% Plot Exact Limit u_0
plot(xCenter, uRefCenter, 'r--', 'LineWidth', 2.5); hold on;
% Plot Numerical Solution u_eps
plot(xCenter, uInterior, 'k-', 'LineWidth', 1.5);
legend('True Limit (u_0)', ['Numerical (u_\epsilon), \epsilon=' num2str(epsilon)], 'Location', 'SouthWest');
title(['One-Sided Problem (\alpha=' num2str(alpha) ', \beta=' num2str(beta) ')']);
xlabel('x'); ylabel('y');
grid on;
xlim([-1 1]);
ylim([-6 8]); % Adjusted to match Figure 5.1 y-axis range roughly

% Figure 2: Convergence History
figure(2)
clf
semilogy(1:iterStop, errorHistory(1:iterStop, 1), 'k-o', 'LineWidth', 1.2, 'MarkerSize', 4); hold on;
semilogy(1:iterStop, errorHistory(1:iterStop, 2), 'r-s', 'LineWidth', 1.2, 'MarkerSize', 4); 
semilogy(1:iterStop, errorHistory(1:iterStop, 3), 'b--', 'LineWidth', 1.5);
legend('Algebraic Residual', 'Correction Norm', 'Physical Diff ||u_\epsilon - u_0||');
title('Solver Convergence History');
xlabel('Cycle'); ylabel('Error / Residual');
grid on;