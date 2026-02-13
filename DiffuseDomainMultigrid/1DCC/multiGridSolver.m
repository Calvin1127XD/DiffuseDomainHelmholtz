function [u,errVals,kStop] =...
              multiGridSolver(u0,f,DeW,C,hf,MGParam,uExact)

nf = length(f);
k = 0;
errorVal = MGParam.tol+1.0;
errVals = zeros(MGParam.kMax,3);

uo = u0;

while ((k < MGParam.kMax) && (errorVal > MGParam.tol))

  % Apply one multigrid cycle.
  u = MGOperator(MGParam.L,f,DeW,C,uo,hf,MGParam);

  % Evaluate various norms.
  residErr = normScaledL2(getResidual(u,f,DeW,C,hf));
  correErr = normScaledL2(uo(2:nf+1)-u(2:nf+1));
  exactErr = normScaledL2(uExact(2:nf+1)-u(2:nf+1));

  k = k+1;
  errVals(k,1) = residErr;
  errVals(k,2) = correErr;
  errVals(k,3) = exactErr;

  fprintf('Multigrid iteration %d:\n', k);
  fprintf('  residual   = %.6e\n', residErr);
  fprintf('  correction = %.6e\n', correErr);
  fprintf('  true error = %.6e\n', exactErr);

  % Stopping criterion (artificial, for demonstration).
  errorVal = correErr;  
  uo = u;
end

kStop = k;

end