function u = smoothRichardson(f,DeW,C,u,hf,m,omega)

nfPlusGhostLayers = length(u);
nf = nfPlusGhostLayers-2;
hf2 = hf*hf;
r = zeros(nf,1);

u = applyBCs(u);

for sweep = 1:m
  r(1:nf) = getResidual(u,f,DeW,C,hf);
  u(2:nf+1) = u(2:nf+1)+omega*(hf2/2.0)*r(1:nf);
  u = applyBCs(u);
end

end