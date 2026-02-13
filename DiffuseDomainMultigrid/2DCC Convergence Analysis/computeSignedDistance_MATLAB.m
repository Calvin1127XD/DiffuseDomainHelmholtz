function phi = computeSignedDistance_MATLAB(X, Y, polyX, polyY, h)
% computeSignedDistance_MATLAB
% A Pure MATLAB implementation of the Fast Marching Method for benchmarking.
%
% Inputs: X, Y (Meshgrids), polyX, polyY (Boundary points), h (Grid spacing)
% Output: Signed Distance Function phi

    [nx, ny] = size(X);
    
    % Constants
    INF_VAL = 1e10; % Large number acting as infinity
    KNOWN = 2;
    TRIAL = 1;
    FAR   = 0;
    
    % Initialize arrays
    phi = ones(nx, ny) * INF_VAL;
    status = zeros(nx, ny); % 0=Far, 1=Trial, 2=Known
    
    % ---------------------------------------------------------
    % 1. Initialization (Exact Distance to Segments)
    % ---------------------------------------------------------
    % We assume polyX, polyY are dense enough or we check bounding boxes.
    % For MATLAB performance, we vectorize the segment check.
    
    bandWidth = 2.5 * h;
    bandWidthSq = bandWidth^2;
    
    % Ensure closed polygon
    if polyX(1) ~= polyX(end) || polyY(1) ~= polyY(end)
        polyX(end+1) = polyX(1);
        polyY(end+1) = polyY(1);
    end
    
    numSegs = length(polyX)-1;
    
    % We will create a mask of points to update to avoid looping every pixel
    % for every segment.
    
    for k = 1:numSegs
        p1 = [polyX(k), polyY(k)];
        p2 = [polyX(k+1), polyY(k+1)];
        
        % Bounding box of segment + band
        min_x = min(p1(1), p2(1)) - bandWidth;
        max_x = max(p1(1), p2(1)) + bandWidth;
        min_y = min(p1(2), p2(2)) - bandWidth;
        max_y = max(p1(2), p2(2)) + bandWidth;
        
        % Find indices (assuming uniform grid structure for speed)
        % This is faster than find(X > ...)
        j_min = find(Y(1,:) >= min_y, 1, 'first');
        j_max = find(Y(1,:) <= max_y, 1, 'last');
        i_min = find(X(:,1) >= min_x, 1, 'first');
        i_max = find(X(:,1) <= max_x, 1, 'last');
        
        if isempty(i_min) || isempty(j_min), continue; end
        
        % Extract local grid patch
        localX = X(i_min:i_max, j_min:j_max);
        localY = Y(i_min:i_max, j_min:j_max);
        
        % Vectorized Point-Segment Distance
        v = p2 - p1;
        len_sq = sum(v.^2);
        
        if len_sq == 0
            t = zeros(size(localX));
        else
            % Dot product projection
            t = ((localX - p1(1))*v(1) + (localY - p1(2))*v(2)) / len_sq;
            t = max(0, min(1, t));
        end
        
        projX = p1(1) + t*v(1);
        projY = p1(2) + t*v(2);
        
        d2 = (localX - projX).^2 + (localY - projY).^2;
        
        % Update phi
        mask = d2 < bandWidthSq;
        
        % Get current values in this patch
        currentPhiPatch = phi(i_min:i_max, j_min:j_max);
        distPatch = sqrt(d2);
        
        % Update min
        updateMask = mask & (distPatch < currentPhiPatch);
        currentPhiPatch(updateMask) = distPatch(updateMask);
        
        phi(i_min:i_max, j_min:j_max) = currentPhiPatch;
    end
    
    % ---------------------------------------------------------
    % 2. Build Initial Trial Set
    % ---------------------------------------------------------
    % In MATLAB, using a linear array for indices is faster than x,y loops
    
    % Identify points initialized in step 1
    initializedMask = phi < bandWidth;
    
    % Mark these as Known? No, usually boundaries are marked known, 
    % neighbors are trial. Here we adopt the strategy:
    % Fixed values are Known. Neighbors are Trial.
    
    [rows, cols] = find(initializedMask);
    linIdxs = sub2ind([nx, ny], rows, cols);
    
    % Set initial known
    status(linIdxs) = KNOWN;
    
    % Lists for Trial Set (Using vectors instead of Heap for MATLAB simplicity)
    trialVal = [];
    trialIdx = [];
    
    % Add neighbors of Known points to Trial
    % Offsets: N, S, E, W
    offsets = [-1, 1, -nx, nx];
    
    % Loop over known points to find initial trials
    % (Vectorizing this logic is messy, doing a quick pass)
    for k = 1:length(linIdxs)
        curr = linIdxs(k);
        [r, c] = ind2sub([nx, ny], curr);
        
        neighbors = [];
        if r > 1,  neighbors(end+1) = curr - 1; end
        if r < nx, neighbors(end+1) = curr + 1; end
        if c > 1,  neighbors(end+1) = curr - nx; end
        if c < ny, neighbors(end+1) = curr + nx; end
        
        for nb = neighbors
            if status(nb) == FAR
                status(nb) = TRIAL;
                val = solveUpwind(nb, phi, status, h, nx);
                phi(nb) = val;
                trialIdx(end+1) = nb; 
                trialVal(end+1) = val;
            elseif status(nb) == TRIAL
                 val = solveUpwind(nb, phi, status, h, nx);
                 if val < phi(nb)
                     phi(nb) = val;
                     % Update value in list
                     mask = (trialIdx == nb);
                     trialVal(mask) = val;
                 end
            end
        end
    end
    
    % ---------------------------------------------------------
    % 3. Fast Marching Loop
    % ---------------------------------------------------------
    % Using standard array operations (O(N) search) which simulates 
    % why pure MATLAB is slow for FMM compared to C++ Heaps.
    
    while ~isempty(trialVal)
        % Find min (Bottleneck in MATLAB without a Heap class)
        [~, minLoc] = min(trialVal);
        currentIdx = trialIdx(minLoc);
        
        % Remove from lists
        trialIdx(minLoc) = [];
        trialVal(minLoc) = [];
        
        if status(currentIdx) == KNOWN
            continue;
        end
        
        status(currentIdx) = KNOWN;
        
        % Update Neighbors
        [r, c] = ind2sub([nx, ny], currentIdx);
        
        neighbors = [];
        if r > 1,  neighbors(end+1) = currentIdx - 1; end
        if r < nx, neighbors(end+1) = currentIdx + 1; end
        if c > 1,  neighbors(end+1) = currentIdx - nx; end
        if c < ny, neighbors(end+1) = currentIdx + nx; end
        
        for nb = neighbors
            if status(nb) ~= KNOWN
                newVal = solveUpwind(nb, phi, status, h, nx);
                
                if status(nb) == FAR
                    status(nb) = TRIAL;
                    phi(nb) = newVal;
                    trialIdx(end+1) = nb;
                    trialVal(end+1) = newVal;
                elseif status(nb) == TRIAL
                    if newVal < phi(nb)
                        phi(nb) = newVal;
                        mask = (trialIdx == nb);
                        trialVal(mask) = newVal;
                    end
                end
            end
        end
    end
    
    % ---------------------------------------------------------
    % 4. Apply Sign
    % ---------------------------------------------------------
    inMask = inpolygon(X, Y, polyX, polyY);
    phi(~inMask) = -phi(~inMask);

end

function val = solveUpwind(idx, phi, status, h, nx)
    % Convert linear index to subscript
    % We assume column-major (standard MATLAB)
    r = mod(idx-1, nx) + 1;
    c = floor((idx-1)/nx) + 1;
    
    INF_VAL = 1e10;
    
    % Get neighbors (Up, Down, Left, Right)
    % Values
    a = INF_VAL; % Min of Horizontal
    b = INF_VAL; % Min of Vertical
    
    % Check Row Neighbors
    if r > 1 && status(idx-1) == 2
        a = min(a, phi(idx-1));
    end
    if r < nx && status(idx+1) == 2
        a = min(a, min(a, phi(idx+1)));
    end
    
    % Check Col Neighbors
    if c > 1 && status(idx-nx) == 2
        b = min(b, phi(idx-nx));
    end
    if (idx+nx) <= numel(phi) && status(idx+nx) == 2 % Bound check for safety
        b = min(b, min(b, phi(idx+nx)));
    end
    
    % Solve Quadratic
    if a == INF_VAL && b == INF_VAL
        val = INF_VAL;
    elseif a == INF_VAL
        val = b + h;
    elseif b == INF_VAL
        val = a + h;
    else
        d = abs(a - b);
        if d >= h
            val = min(a, b) + h;
        else
            val = 0.5 * (a + b + sqrt(2*h^2 - d^2));
        end
    end
end