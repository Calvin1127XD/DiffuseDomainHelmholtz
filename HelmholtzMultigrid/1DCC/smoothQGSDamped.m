function u = smoothQGSDamped(f,DeW,C,u,hf,m,omega,...
                             sweepDirection)

nfPlusGhostLayers = length(u);
nf = nfPlusGhostLayers-2;

hf2 = hf*hf;
omegaPrime = 1.0-omega;

switch sweepDirection
  case 'fwd'
    sweepIndex = [1:nf];
  otherwise
    sweepIndex = [nf:-1:1];
end

u = applyBCs(u);

for sweep = 1:m
  for i = sweepIndex

    tmp = (hf2*f(i)+u(i+2)*DeW(i+1)+u(i)*DeW(i))...
          /(DeW(i+1)+DeW(i)+C(i)*hf2);
    u(i+1) = omega*tmp+omegaPrime*u(i+1);

  end
  u = applyBCs(u);

end

end