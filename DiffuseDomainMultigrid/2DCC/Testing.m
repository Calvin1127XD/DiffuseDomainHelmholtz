% RunGlobal_Save.m
% Runs the original Global Solver with the specific test case parameters
% and saves the result for consistency checking.
%
% COMPATIBILITY NOTE: This script is designed for the ORIGINAL library
% (8 arguments, uExact required).

clear; clc;

%% 1. Simulation Parameters
L = 12; 
nL = [2^L, 2^L]; 
N = nL(1);

xyLower = [-2.0, -2.0];
xyUpper = [ 2.0,  2.0];
hMesh = (xyUpper(1) - xyLower(1)) / nL(1);

% Problem Parameters (Specific Test Case)
epsilon = 0.05; 
alpha   = 3.0;
beta    = 2.0;
gamma   = 1.0;
kappa   = 0.01;

% Multigrid Parameters
MGParam.nL     = nL;
MGParam.xLower = xyLower;
MGParam.xUpper = xyUpper;
MGParam.L      = L;
MGParam.pCycle = 1; MGParam.m1 = 3; MGParam.m2 = 3;      
MGParam.omega  = 0.66; MGParam.kMax = 100; MGParam.tol = 1.0e-12;

fprintf('Running Global Solver (N=%d)...\n', N);

%% 2. Grid Setup
xCenter = linspace(xyLower(1) + hMesh/2, xyUpper(1) - hMesh/2, nL(1));
yCenter = linspace(xyLower(2) + hMesh/2, xyUpper(2) - hMesh/2, nL(2));
[XCC, YCC] = ndgrid(xCenter, yCenter);

xEdgePts = linspace(xyLower(1), xyUpper(1), nL(1)+1);
yEdgePts = linspace(xyLower(2), xyUpper(2), nL(2)+1);

[XEW, YEW] = ndgrid(xEdgePts, yCenter);
[XNS, YNS] = ndgrid(xCenter, yEdgePts);

%% 3. Geometry
[polyX, polyY] = generateBoundary(2000);
distCenter = computeSignedDistance(XCC, YCC, polyX, polyY, hMesh);
distEW = computeSignedDistance(XEW, YEW, polyX, polyY, hMesh);
distNS = computeSignedDistance(XNS, YNS, polyX, polyY, hMesh);

%% 4. Compute Coefficients (Exact Match to Test Case)
% Phase Field
argEW = distEW ./ epsilon; phiEW = 0.5 .* (1.0 + tanh(argEW));
diffCoefEW = alpha + (1.0 - alpha) .* phiEW;

argNS = distNS ./ epsilon; phiNS = 0.5 .* (1.0 + tanh(argNS));
diffCoefNS = alpha + (1.0 - alpha) .* phiNS;

argCC = distCenter ./ epsilon;
phiCC = 0.5 .* (1.0 + tanh(argCC));
absGradPhi = (1.0 / (2.0 * epsilon)) .* (sech(argCC).^2);

% Reaction
cEps = beta .* (1.0 - phiCC) + gamma .* phiCC;
reacCoef = cEps + kappa .* absGradPhi;

% Forcing Terms (Test Case Specifics)
qFunc = -XCC.^2 + 15.0;
hFunc = 2.5*sin(XCC)+exp(cos(YCC));
fEps = hFunc .* (1.0 - phiCC) + qFunc .* phiCC;

gFunc = XCC - XCC + 4; % Effectively 4.0
rhsVector = fEps - gFunc .* absGradPhi;

%% 5. Run Solver
uInit = zeros(nL(1)+2, nL(2)+2);

% DUMMY UREF: The old library requires the 8th argument to calculate error.
% We pass a zero matrix so it doesn't crash (error calculation will be junk but we ignore it).
uRefDummy = zeros(nL(1)+2, nL(2)+2);

% Call with 8 arguments (Old Signature)
[uSol, ~, ~] = multiGridSolver(uInit, rhsVector, diffCoefEW, ...
                               diffCoefNS, reacCoef, hMesh, MGParam, uRefDummy);

uInterior = uSol(2:end-1, 2:end-1);

%% 6. Save Data
save('GlobalSolution.mat', 'uInterior', 'xCenter', 'yCenter', 'hMesh', 'N');
fprintf('Saved Global solution to GlobalSolution.mat\n');