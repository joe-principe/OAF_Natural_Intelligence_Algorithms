% Script file: principe_oaf_week2.m
%
% Purpose:
%   To demonstrate the Particle Swarm Optimization (PSO) algorithm
%
% Define variables:
%   X            -- The x-coordinates from [0, 5]
%   Y            -- The y-coordinates from [0, 5]
%   Z            -- The z-coordinates of the objective function
%
%   x_min        -- The x-coordinate of the global minimum
%   x_max        -- The x-coordinate of the global maximum
%   y_min        -- The y-coordinate of the global minimum
%   y_max        -- The y-coordinate of the global maximum
%
%   g_min        -- The global minimum value
%   g_max        -- The global maximum value
%   min_index    -- The index of the global minimum value of the object function
%   max_index    -- The index of the global maximum value of the object function
%
%   num_ind      -- The number of individuals in the swarm
%   num_iter     -- The number of times to iterate through the algorithm
%   P            -- An array of the positions of each individual [X; Y]
%   V            -- An array of the velocities of each individual [X; Y]
%   new_pos      -- The new position of each individual
%
%   w            -- The inertia weight coefficient
%   c1           -- The cognitive constant
%   c2           -- The social constant
%   r            -- An array of random numbers between [0, 1]
%
%   p_best       -- An array of the (x, y) coordinates of the personal best of each individual
%   p_best_obj   -- An array of the personal best height of each individual
%   g_best       -- An array of the (x, y) coordinates of the personal best of each individual
%   g_best_obj   -- The current global best height
%   g_best_index -- The index of the individual with the global best value

% Clear the workspace
clear;

% Initialize the constants and coefficients
w = 0.8;
c1 = 0.1;
c2 = 0.1;

% Initialize and plot the objective function
[X, Y] = meshgrid(linspace(0, 5, 100), linspace(0, 5, 100)); % Creates a grid of (X, Y) coordinate pairs
Z = obj_func(X, Y);
subplot(2, 2, 1); % Creates a 2x2 area for plots
contourf(X, Y, Z); % Creates a filled contour plot of the objective function
title('Initial positions of individuals');
colorbar;
hold on;

% Find the global minimum and maximum and their coordinates
% min_index and max_index are used to find the coordinates of the global min/max
[g_min, min_index] = min(Z, [], 'all', 'linear');
% [g_max, max_index] = max(Z, [], 'all', 'linear');

x_min = X(min_index);
y_min = Y(min_index);

% x_max = X(max_index);
% y_max = Y(max_index);

% Initialize the individuals and data
num_ind = 24;
P = rand(2, num_ind) * 5; % The first row are the X values, the second are the Y values
V = randn(2, num_ind) * 0.1; % The first row are the X values, the second are the Y values
plot(P(1, :), P(2, :), 'wo', 'MarkerFaceColor', 'k'); % Plots the starting positions of the individuals
plot(x_min, y_min, 'cx'); % Plots the global minimum (cyan cross)
hold off;

p_best = P; % Sets the the starting location as the initial personal best
p_best_obj = obj_func(P(1, :), P(2, :)); % Gets the personal best of each individual
[~, g_best_index] = min(p_best_obj); % Gets the index of the individual with the global best position
g_best = p_best(:, g_best_index); % Gets the location of the global best
g_best_obj = p_best_obj(:, g_best_index); % Gets the global best value

% Iterates through the algorithm and plots the positions every 10 iterations
num_iter = 30;
for ii = 1:num_iter
    r = rand(1, 2); % Gets random values for r1 and r2
    V = w * V + c1 * r(1) * (p_best - P) + c2 * r(2) * (reshape(g_best, 2, 1) - P);
    P = P + V;
    new_height = obj_func(P(1, :), P(2, :)); % Gets the new z-coordinate of each individual
    p_best(repmat(p_best_obj >= new_height, 2, 1)) = P(repmat(p_best_obj >= new_height, 2, 1)); % Updates the personal best location for individuals with a new personal best
    p_best_obj = min(p_best_obj, new_height); % Sets this iterations personal best value to the minimum between the personal best value and current height
    [~, g_best_index] = min(p_best_obj);
    g_best = p_best(:, g_best_index);
    g_best_obj = min(p_best_obj, [], 'all');
    
    % Plots the iterations of the algorithm
    if (mod(ii, 10) == 0)
        subplot(2, 2, ii/10 + 1);
        hold on;
        contourf(X, Y, Z);
        title(['Individuals at iteration ', num2str(ii)]);
        colorbar;
        plot(P(1, :), P(2, :), 'wo', 'MarkerFaceColor', 'k'); % Plots the starting positions of the individuals
        plot(x_min, y_min, 'cx'); % Plots the global minimum (white cross) and maximum (red circle)
        hold off;
    end
    
end

% Calculates the objective function
function height = obj_func(x, y)
narginchk(2, 2);

height = (x - 3.14).^2 + (y - 2.72).^2 + sin(3 * x + 1.41) + sin(4 * y - 1.73);
% height = sin(x).^2 + sin(y).^2 + sin(x) .* sin(y);
% height = x.^2 + (y + 1).^2 - 5 * cos(1.5 * x + 1.5) - 3 * cos(2 * x - 1.5);
end