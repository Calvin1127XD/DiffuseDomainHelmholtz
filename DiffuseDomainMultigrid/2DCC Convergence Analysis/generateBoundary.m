function [x, y] = generateBoundary(n)
% GENERATEBOUNDARY Generates boundary coordinates for the interface.
%
%   Input:
%       n : Number of points to discretize the boundary.
%
%   Output:
%       x, y : Vectors of coordinates defining the closed polygon.
%
%   Shape: 
%       Zero level set of P(x,y) = (x^2+y^2)^2 + 2(x^3 - 3xy^2) - 2.5
%       Polar form: r^4 + 2*r^3*cos(3*theta) - 2.5 = 0
%       Note: The 'tip' of this shape extends to r ~ 2.16.
%

    theta = linspace(0, 2*pi, n);
    r = zeros(size(theta));
    
    % The constant term in the polynomial
    poly_const = 2.5;
    
    % Solve for radius r at each angle theta
    % The root is known to exist in [0.5, 3.0] for this specific shape.
    options = optimset('Display', 'off');
    
    for i = 1:n
        ct = cos(3 * theta(i));
        
        % Function f(r) = r^4 + 2*cos(3t)*r^3 - 2.5
        fun = @(rad) rad^4 + 2*ct*rad^3 - poly_const;
        
        r(i) = fzero(fun, [0.5, 3.0], options);
    end
    
    % Convert to Cartesian
    x = r .* cos(theta);
    y = r .* sin(theta);

end