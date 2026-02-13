function operatorResult = FDOperator(u,DeW,C,hf)

nfPlusGhostLayers = length(u);
nf = nfPlusGhostLayers-2;
feW = zeros(nf+1,1);
operatorResult = zeros(nf,1);

hf2 = hf*hf;

% Vectorized code for the operator:
% - (D u')' + C(x) u 

% Compute east-west numerical flux:
feW(1:nf+1) = DeW(1:nf+1).*(u(2:nf+2)-u(1:nf+1));

operatorResult(1:nf) = -(feW(2:nf+1)-feW(1:nf))/hf2...
                       +C.*u(2:nf+1);

end