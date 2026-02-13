function residualResult = getResidual(u,f,DeW,DnS,C,hf)

% Apply the boundary conditions
u = applyBCs(u);

residualResult = f-FDOperator(u,DeW,DnS,C,hf);

end