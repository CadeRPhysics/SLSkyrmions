% Matlab code created by Nilo Mata Cervera on February 2025
% Image processing algorithm for finding the region of minimum S3
% The code starts with drawing a circle centered at
% [cx,cy] and radius R (all in pixel units). Number of
% points of the initial circle "Npoints", number of
% iterations "Nit", "alpha" weights the elastic forces and "beta" the S3 forces
% The interpolation "method" must be specified: "linear", "spline", etc
% "navg" number of points to smooth the contour to avoid unstabilities

function [mask,min_region] = boundary(S3,cx,cy,R,Npoints,Nit,alpha,beta,gamma,method,navg)

    num_points = Npoints; % Number of points for the initial circle

    % Generate initial circular boundary points
    theta = linspace(0, 2*pi, num_points);
    x_circle = round(cx + R * cos(theta));
    y_circle = round(cy + R * sin(theta));

    % Ensure points are within image bounds
    [rows, cols] = size(S3);
    x_circle = max(min(x_circle, cols), 1);
    y_circle = max(min(y_circle, rows), 1);

    % Initialize the refined boundary points
    x_refined = x_circle;
    y_refined = y_circle;

    % Iterative deformation process
    num_iterations = Nit; % Controls how much the shape evolves
    % alpha = 0.0;  % Weight for internal smoothness (elastic force)
    % beta = 50.0;   % Weight for external force (S3 gradient)

    for iter = 1:num_iterations
        % Compute external forces (gradient of S3)
        [grad_x, grad_y] = gradient(S3);
        
        % Interpolate gradient at contour points
        fx = interp2(grad_x, x_refined, y_refined, method, 0);
        fy = interp2(grad_y, x_refined, y_refined, method, 0);
        
        % Compute internal forces (continuity and smoothness)
        x_smooth = circshift(x_refined, -1) + circshift(x_refined, 1) - 2*x_refined;
        y_smooth = circshift(y_refined, -1) + circshift(y_refined, 1) - 2*y_refined;

        % Compute curvature forces
        x_K = -circshift(x_refined, 2) + 4*circshift(x_refined, 1) - ...
            x_refined + 4*circshift(x_refined, -1)-circshift(x_refined, -2); 
        y_K = -circshift(y_refined, 2) + 4*circshift(y_refined, 1) - ...
            y_refined + 4*circshift(y_refined, -1)-circshift(y_refined, -2);
        
        % Update points using a weighted sum of forces, skyrmion force 
        x_refined = x_refined + alpha * x_smooth - beta * fx + gamma * x_K;
        y_refined = y_refined + alpha * y_smooth - beta * fy + gamma * y_K;

        
        % Update points using a weighted sum of forces
        x_refined = x_refined + alpha * x_smooth - beta * fx;
        y_refined = y_refined + alpha * y_smooth - beta * fy;
        
        % Keep within image bounds
        x_refined = max(min(round(x_refined), cols), 1);
        y_refined = max(min(round(y_refined), rows), 1);
        
        % Optional: Apply smoothing to avoid jagged movements
        x_refined = smooth(x_refined, navg);
        y_refined = smooth(y_refined, navg);
    end

    % Store the final refined boundary (array of points)
    min_region = [x_refined'; y_refined'];

    % Convert boundary to binary mask (1 inside and 0 outside)
    mask = poly2mask(x_refined, y_refined, rows, cols);

end