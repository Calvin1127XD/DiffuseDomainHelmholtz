function u = MGOperator(level,f,DeW,C,u,hf,MGParam)

% Pre-smoothing
% u = smoothQJacDamped(f,DeW,C,u,hf,MGParam.m1,...
%                      MGParam.omega);

% Alternative pre-smoothing options:
% u = smoothRichardson(f,DeW,C,u,hf,MGParam.m1,...
%                      MGParam.omega);
u = smoothQGSDamped(f,DeW,C,u,hf,MGParam.m1,...
                    MGParam.omega,'fwd');

if level > 0
  hc = 2*hf;
  nfPlusGhost = length(u);
  nf = nfPlusGhost-2;

  % Coarse grid interior dimension is nc.
  nc = nf/2;

  % Compute the residual and restric to the coarse grid.
  cGr = restriction(getResidual(u,f,DeW,C,hf));
    
  % Restrict the diffusion coefficients to the coarse
  % grid.
  cDeW = restrictDiff(DeW);
    
  % Restrict the reaction coefficient to the coarse grid.
  cC = restriction(C);
  % Injection based restriction.
  % cC = C(1:2:end);
  
  % Approximate on the coarse grid using recursive MG.
  cGc = zeros(nc+2,1);
  for s = 1:MGParam.pCycle
    cGc = MGOperator(level-1,cGr,cDeW,cC,cGc,hc,MGParam);
  end

  % Prolongate the coarse grid correction and update the
  % fine grid approximation.
  u(2:nf+1) = u(2:nf+1)+prolongation(cGc(2:nc+1));
    
  % Post-smoothing.
  % u = smoothQJacDamped(f,DeW,C,u,hf,MGParam.m2,...
  %                      MGParam.omega);

  % Alternative post-smoothing options:
  % u = smoothRichardson(f,DeW,C,u,hf,MGParam.m2,...
  %                      MGParam.omega);
  u = smoothQGSDamped(f,DeW,C,u,hf,MGParam.m2,...
                      MGParam.omega,'bwd');
end

end