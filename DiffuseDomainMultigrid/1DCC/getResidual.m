function residualResult = getResidual(u,f,DeW,C,hf)

% Apply the boundary conditions before residual 
% computation.
u = applyBCs(u);
nf = length(f);
residualResult = zeros(nf,1);

residualResult = f-FDOperator(u,DeW,C,hf);

end