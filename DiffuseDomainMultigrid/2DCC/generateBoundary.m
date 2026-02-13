function [x, y] = generateBoundary(n)
    % Generates boundary coordinates for the interface.
    
    theta = linspace(0, 2*pi, n);

    % --- Case 1: Ellipse (PDF Section 4.2, First Example) ---
    % Ellipse: x^2 + 4y^2 = 4  =>  (x/2)^2 + (y/1)^2 = 1
    % Semi-major axis a=2, Semi-minor axis b=1
    
    % x = 2.0 * cos(theta);
    % y = 1.0 * sin(theta);
    
    % --- Case 2: Starfish (PDF Section 4.2, Equation 4.9) ---
    % r(theta) = 0.9 * (1.2 + 0.7 * sin(5*theta))
    % This tests the solver on non-convex geometries.
    
    r = 0.9 * (1.2 + 0.7 * sin(5 * theta));
    
    x = r .* cos(theta);
    y = r .* sin(theta);
end