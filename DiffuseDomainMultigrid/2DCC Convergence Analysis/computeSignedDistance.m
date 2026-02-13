function phi = computeSignedDistance(X, Y, polyX, polyY, h)
    
    % Ensure polygon is closed
    if polyX(1) ~= polyX(end) || polyY(1) ~= polyY(end)
        polyX(end+1) = polyX(1);
        polyY(end+1) = polyY(1);
    end

    % Prefer compiled MEX for speed; fall back to pure MATLAB when unavailable.
    if exist('mex_computeSignedDistance', 'file') == 3
        phi = mex_computeSignedDistance(X, Y, polyX, polyY, h);
    elseif exist('computeSignedDistance_MATLAB', 'file') == 2
        warning('computeSignedDistance:NoMex', ...
            ['mex_computeSignedDistance not found. Falling back to ', ...
             'computeSignedDistance_MATLAB (slower).']);
        phi = computeSignedDistance_MATLAB(X, Y, polyX, polyY, h);
    else
        error(['Neither mex_computeSignedDistance nor ', ...
               'computeSignedDistance_MATLAB is available.']);
    end

end
