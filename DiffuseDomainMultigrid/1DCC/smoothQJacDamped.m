function u = smoothQJacDamped(f,DeW,C,u,hf,m,omega)

nfPlusGhostLayers = length(u);
nf = nfPlusGhostLayers-2;

z = zeros(nf,1);

hf2 = hf*hf;
omegaPrime = 1.0-omega;

u = applyBCs(u);

for sweep = 1:m
  for i = 1:nf

    z(i) = (hf2*f(i)+u(i+2)*DeW(i+1)+u(i)*DeW(i))...
           /(DeW(i+1)+DeW(i)+C(i)*hf2);

  end

  u(2:nf+1) = omega*z(1:nf)+omegaPrime*u(2:nf+1);
  u = applyBCs(u);
  
end

end