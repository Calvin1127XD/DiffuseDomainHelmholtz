function uInterior = computeDiffuseDomainSolution2D(gridInfo, physParam, precompData, solverParam)
% computeDiffuseDomainSolution2D
%
% Solves the DDM Starfish problem for a specific set of parameters.
%
% OPTIMIZATION:
% All spatially-dependent but parameter-independent fields (SDF, q, h, g)
% are passed in 'precompData' to avoid redundant calculation.
%
% Inputs:
%   gridInfo    : Struct (nL, hMesh)
%   physParam   : Struct (epsilon, alpha, beta, gamma, kappa)
%   precompData : Struct (distCenter, distEW, distNS, qFunc, hFunc, gFunc)
%   solverParam : Struct (MG settings)

    % Unpack Grid Params
    nL = gridInfo.nL;
    hMesh = gridInfo.hMesh;

    % Unpack Physics Params
    epsilon = physParam.epsilon;
    alpha   = physParam.alpha;
    beta    = physParam.beta;
    gamma   = physParam.gamma;
    kappa   = physParam.kappa;
    
    % ---------------------------------------------------------
    % 1. Compute Coefficients (Depend on Epsilon/Alpha)
    % ---------------------------------------------------------

    % --- Diffusion Coefficient D (at Edges) ---
    % D = alpha + (1-alpha)*phi
    % precompData.distEW/NS are passed in directly
    
    phiEW = 0.5 .* (1.0 + tanh(precompData.distEW ./ epsilon));
    diffCoefEW = alpha + (1.0 - alpha) .* phiEW;

    phiNS = 0.5 .* (1.0 + tanh(precompData.distNS ./ epsilon));
    diffCoefNS = alpha + (1.0 - alpha) .* phiNS;

    % --- Reaction Coefficient C (at Centers) ---
    argCC = precompData.distCenter ./ epsilon;
    
    % Optimization: Compute common terms once
    tanhPhi = tanh(argCC);
    phiCC = 0.5 .* (1.0 + tanhPhi);

    % Gradient Magnitude: (1/(2eps)) * (1 - tanh^2)
    % sech^2(x) = 1 - tanh^2(x)
    absGradPhi = (1.0 / (2.0 * epsilon)) .* (1.0 - tanhPhi.^2);

    % c_eps = beta*(1-phi) + gamma*phi
    cEps = beta .* (1.0 - phiCC) + gamma .* phiCC;

    % Total Operator C = c_eps + kappa * |grad phi|
    reacCoef = cEps + kappa .* absGradPhi;

    % ---------------------------------------------------------
    % 2. Compute RHS Forcing
    % ---------------------------------------------------------
    
    % f_eps = h*(1-phi) + q*phi
    % hFunc and qFunc are pre-computed!
    fEps = precompData.hFunc .* (1.0 - phiCC) + precompData.qFunc .* phiCC;

    % Total RHS = f_eps - g * |grad phi|
    % gFunc is pre-computed
    rhsVector = fEps - precompData.gFunc .* absGradPhi;

    % ---------------------------------------------------------
    % 3. Run Multigrid Solver
    % ---------------------------------------------------------
    
    MGParam = solverParam;
    MGParam.nL = nL;
    MGParam.C = reacCoef; 
    
    uInit = zeros(nL(1)+2, nL(2)+2);
    uRef  = zeros(nL(1)+2, nL(2)+2); 

    [uSol, ~, ~] = multiGridSolver(uInit, rhsVector, diffCoefEW, ...
                                   diffCoefNS, reacCoef, hMesh, ...
                                   MGParam, uRef);

    uInterior = uSol(2:end-1, 2:end-1);
end