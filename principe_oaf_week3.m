% Script file: principe_oaf_week3.m
%
% Purpose:
%   To demonstrate the Particle Swarm Optimization (PSO) algorithm
%
% Define variables:
%   X -- The x-axis search space boundaries [-3, 3], comprised of 100
%   points
%   Y -- The y-axis search space boundaries [-3, 3], comprised of 100
%   points
%   Z -- The objective function f(x, y) = sin(x)^2 + sin(y)^2 + sin(x)sin(y)
%   P -- The position of each particle. Rows are particles, columns are X
%   and Y coordinates
%   V -- The velocity of each particle. Rows are particles, columns are X
%   and Y velocities
%   p_best -- The position of the personal best of each particle. Rows are
%   particles, columns are X and Y coordinates
%   fit_p_best -- The personal best height of each particle
%   g_best -- The global best position of the swarm
%   fit_g_best -- The global best height of the swarm
%   num_iter -- Number of iterations for the PSO algorithm
%   num_particles -- The number of particles in the swarm
%   v_max -- The maximum velocity of each particle

% Clear the workspace
clc;
clear;

% Initialize the constants and coefficients
w = 0.75;
c1 = 0.2;
c2 = 0.2;
num_iter = 100;
num_particles = 30;
v_max = 1;

% Initialize the objective function
[X, Y] = meshgrid(linspace(-3, 3, 100), linspace(-3, 3, 100));
Z = obj_func(X, Y); 
subplot(3, 2, 1); % Creates a 3x2 area for plots
contourf(X, Y, Z); % Creates a filled contour plot of the objective function
title('Initial positions of particles');
colorbar;
hold on;

%%%%% For plotting only %%%%%
% Locate the global minimum and its coordinates
[g_min, min_index] = min(Z, [], 'all', 'linear'); % Returns only the smallest value in Z and its index
x_min = X(min_index); % Gets the x-coordinate of the global minimum
y_min = Y(min_index); % Gets the y-coordinate of the global minimum
%%%%% For plotting only %%%%%

% Initialize the positions of particles
P = -3 + (6 - 1) .* rand(num_particles, 2); % Randomly distributes the particles between [-3, 3]
V = 0.8 * randn(num_particles, 2); % Randomly sets the velocity between [0, 1]
V(V > v_max) = v_max; % Sets velocities > v_max to v_max
V(V < -v_max) = -v_max; % Sets velocities < v_min (aka -v_max) to v_min
plot(P(:, 1), P(:, 2), 'wo', 'MarkerFaceColor', 'k'); % Plots the particles' starting points
plot(x_min, y_min, 'cx'); % Plots the global minimum;
hold off;

p_best = P; % Sets the starting location as the initial personal best since this is the only place the particle has been
fit_p_best = obj_func(P(:, 1), P(:, 2)); % Sets the starting height as the initial personal best

% Find the global best
[~, g_best_index] = min(fit_p_best); % Gets the index of the minimum personal best
g_best = p_best(g_best_index, :); % Sets the position of the global best to be the position of the lowest personal best
fit_g_best = fit_p_best(g_best_index, :); % Sets the global best height to be the lowest personal best height

for ii = 1:num_iter
    w = 1.2 - ((1.1 / num_iter) * (ii + 1)); % Reduces the inertia over time
    
    V = w * V + c1 * rand * (p_best - P) + c2 * rand * (g_best - P); % Updates the velocity
    V(V > v_max) = v_max; % Sets velocities > v_max to v_max
    V(V < -v_max) = -v_max; % Sets velocities < v_min (aka -v_max) to v_min
    P = P + V; % Updates the particles' positions
    
    for jj = 1:num_particles
        if P(jj, 1) > 3 && P(jj, 1) < -3 && P(jj, 2) > 3 && P(jj, 2) < -3
            new_height = obj_func(P(:, 1), P(:, 2)); % Gets the new height of each particle
            p_best(repmat(fit_p_best >= new_height, 2, 1)) = P(repmat(fit_p_best >= new_height, 2, 1)); % Updates the personal best location for individuals with a new personal best
            % Checks for if the new height is lower than the previous best, if
            % it is, then we set the personal best to the new height
            fit_p_best = min(fit_p_best, new_height);

            [~, g_best_index] = min(fit_p_best);
            g_best = p_best(g_best_index, :);
            fit_g_best = fit_p_best(g_best_index, :);
        end
    end
    
    G_best(ii, :) = g_best;
    Fitness(ii) = fit_g_best;
    
    if (mod(ii, 25) == 0)
        subplot(3, 2, ii/25 + 1);
        hold on;
        contourf(X, Y, Z);
        title(['Particles at iteration ', num2str(ii)]);
        colorbar;
        plot(P(:, 1), P(:, 2), 'wo', 'MarkerFaceColor', 'k'); % Plots the starting positions of the individuals
        plot(x_min, y_min, 'cx'); % Plots the global minimum (white cross) and maximum (red circle)
        hold off;
    end
end

G_best
Fitness'
minimum = min(Fitness)
maximum = max(Fitness)
average = mean(Fitness)
std_dev = std(Fitness)
median = median(Fitness)

function result = obj_func(x, y)
% Checks for the correct number of inputs
narginchk(2, 2);

result = sin(x.^2) + sin(y.^2) + sin(x) .* sin(y);
end