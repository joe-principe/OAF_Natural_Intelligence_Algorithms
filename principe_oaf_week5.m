% Script file: principe_oaf_week5.m
%
% Define variables:
%   w                 -- The inertia weight coefficient
%   c1                -- The cognitive constant
%   c2                -- The social constant
%
%   num_iter          -- The number of iterations for the PSO algorithm
%   num_particles     -- The number of particles for the PSO algorithm
%   num_entries       -- The number of entries per iris type
%   num_entries_total -- The total number of entries in the data file
%   num_samples       -- The number of samples for pattern recognition
%   num_classes       -- The number of classes in the iris data set
%   v_max             -- The maximum velocity of the particles
%   Fitness           -- The fitness value of each centroid
%
%   p_best            -- The personal best position of each particle
%   fit_p_best        -- The personal best fit of each particle
%   g_best            -- The global best position of each particle
%   fit_g_best        -- The global best fit of each particle
%
%   centroids         -- The positions of each centroid. First row is
%                        setosa, second versicolor, third virginica
%
%   data_file         -- The iris data file
%   data              -- The iris data as a cell array
%   sepal_length      -- The iris sepal lengths (cm)
%   sepal_width       -- The iris sepal widths (cm)
%   petal_length      -- The iris petal lengths (cm)
%   petal_width       -- The iris petal widths (cm)
%   measurements      -- An array containing all the measurements
%   classes           -- The iris classes (setosa, versicolor, and
%                        virginica)
%
%   P                 -- The positions of the particles. Rows are
%                        particles, columns are coordinates (sl,sw,pl,pw)
%   P_sl              -- The position on the sepal length axis
%   P_sw              -- The position on the sepal width axis
%   P_pl              -- The position on the petal length axis
%   P_pw              -- The position on the petal width axis
%
%   V                 -- The velocities of the particles. Rows are
%                        particles, columns are velocities (sl,sw,pl,pw)
%   V_sl              -- The velocity on the sepal length axis
%   V_sw              -- The velocity on the sepal width axis
%   V_pl              -- The velocity on the petal length axis
%   V_pw              -- The velocity on the petal width axis
%
%   S                 -- The positions of the samples. Rows are samples,
%                        columns are coordinates (sl,sw,pl,pw)
%   S_sl              -- The position on the sepal length axis
%   S_sw              -- The position on the sepal width axis
%   S_pl              -- The position on the petal length axis
%   S_pw              -- The position on the petal width axis
%   sample_classes    -- The classes of each sample
%
%   fit               -- New height for the current particle
%   distances         -- An array of the distance from each sample to each
%                        centroid
%   dist_setosa       -- The distance from the current sample to the setosa
%                        centroid
%   dist_versicolor   -- The distance from the current sample to the
%                        versicolor centroid
%   dist_virginica    -- The distance from the current sample to the
%                        virginica centroid

% Clear the workspace
clc;
clear;

% Initialize constants
w = 0.75;
c1 = 0.2;
c2 = 0.2;
num_iter = 100;
num_particles = 30;
num_entries = 50;
num_entries_total = 150;
num_samples = 30;
num_classes = 3;
v_max = 1;
% Fitness = zeros(num_classes, 1);

% Initialize the particle parameters
P = zeros(num_particles, 4);
V = zeros(num_particles, 4);
p_best = zeros(num_particles, 4);
fit_p_best = zeros(num_particles, 1);
g_best = zeros(1, 4);
fit_g_best = 0;

% Initialize centroid positions at (0, 0, 0, 0)
centroids = zeros(3, 4);

% Load the data file
data_file = fopen('iris.data');
data = textscan(data_file, '%f,%f,%f,%f,%s');
fclose(data_file);

% Get the attributes and classes from the data file
sepal_length = cell2mat(data(1));
sepal_width = cell2mat(data(2));
petal_length = cell2mat(data(3));
petal_width = cell2mat(data(4));
measurements = [sepal_length, sepal_width, petal_length, petal_width];
classes = unique(string(data{5}));

% Generate the samples for pattern recognition
[S_sl, S_sw, S_pl, S_pw] = set_starting_pos(sepal_length,...
    sepal_width, petal_length, petal_width, num_samples);
S = [S_sl, S_sw, S_pl, S_pw];

% Creates an empty string array for the samples' classes
sample_classes = strings(num_samples, 1);

% Run the PSO algorithm for each cluster
% First setosa, then versicolor, then virginica
for ii = 1:num_classes
    % Initialize the starting positions and velocities of the particles
    % Done in this for loop because new particles are needed
    % for each centroid
    [P_sl, P_sw, P_pl, P_pw] = set_starting_pos(sepal_length,...
        sepal_width, petal_length, petal_width, num_particles);
    P = [P_sl, P_sw, P_pl, P_pw];
    [V_sl, V_sw, V_pl, V_pw] = set_starting_vel(num_particles);
    V = [V_sl, V_sw, V_pl, V_pw];

    % Set the starting personal best position and fitness
    p_best = P;
    fit_p_best = min(fit_func(P, measurements, num_entries, ii), [], 2);

    % Set the starting global best position and fitness
    [~, g_best_index] = min(fit_p_best);
    g_best = p_best(g_best_index, :);
    fit_g_best = fit_p_best(g_best_index, :);

    % PSO algorithm (running 100 times)
    for jj = 1:num_iter
        for kk = 1:num_particles
            % Adjust the inertia each iteration
            w = 1.2 - ((1.1 / num_iter) * (jj + 1));

            % Set the new velocity for each axis
            % Sepal length axis vel
            V(kk, 1) = w * V(kk, 1) + ...
                       c1 * rand * (p_best(kk, 1) - P(kk, 1)) + ...
                       c2 * rand * (g_best(1) - P(kk, 1));
            % Sepal width axis vel
            V(kk, 2) = w * V(kk, 2) + ...
                       c1 * rand * (p_best(kk, 2) - P(kk, 2)) + ...
                       c2 * rand * (g_best(2) - P(kk, 2));
            % Petal length axis vel
            V(kk, 3) = w * V(kk, 3) + ...
                       c1 * rand * (p_best(kk, 3) - P(kk, 3)) + ...
                       c2 * rand * (g_best(3) - P(kk, 3));
            % Petal width axis vel
            V(kk, 4) = w * V(kk, 4) + ...
                       c1 * rand * (p_best(kk, 4) - P(kk, 4)) + ...
                       c2 * rand * (g_best(4) - P(kk, 4));

            % Set the new position for each axis
            P(kk, 1) = P(kk, 1) + V(kk, 1);
            P(kk, 2) = P(kk, 2) + V(kk, 2);
            P(kk, 3) = P(kk, 3) + V(kk, 3);
            P(kk, 4) = P(kk, 4) + V(kk, 4);

            % Calculate the new fitness if the particle is in bounds
            if (P(kk, 1) > 0 && P(kk, 1) < max(sepal_length))
                if (P(kk, 2) > 0 && P(kk, 2) < max(sepal_width))
                    if (P(kk, 3) > 0 && P(kk, 3) < max(petal_length))
                        if (P(kk, 4) > 0 && P(kk, 4) < max(petal_width))
                            fit = fit_func(P(kk, :), measurements, num_entries, ii);

                            if (min(fit) < fit_p_best(kk))
                                fit_p_best(kk) = min(fit);
                                p_best(kk, :) = P(kk, :);
                            end
                        end
                    end
                end
            end
        end
        % Calculate the new g_best and fit_g_best
        [~, g_best_index] = min(fit_p_best);
        g_best = p_best(g_best_index, :);
        fit_g_best = fit_p_best(g_best_index, :);
    end
    centroids(ii, :) = g_best;
%     Fitness(ii) = fit_g_best;
end

% Print the locations of the centroids
for ii = 1:num_classes
    fprintf('The %-15s centroid is located at: (%.4f, %.4f, %.4f, %.4f)\n',...
        classes(ii), centroids(ii, 1), centroids(ii, 2),...
                     centroids(ii, 3), centroids(ii, 4));
%     fprintf(' with a fitness value of %.4f\n', Fitness(ii));
end

distances = fit_func(S, centroids, 3, 1);

% Run the pattern recognition algorithm
for ii = 1:num_samples
    % Determine distance to all centroids
    dist_setosa     = fit_func(S(ii, :), centroids, 1, 1);
    dist_versicolor = fit_func(S(ii, :), centroids, 2, 1);
    distances = [dist_setosa; dist_versicolor; dist_virginica];
    dist_virginica  = fit_func(S(ii, :), centroids, 3, 1);

    % Centroid closest to current sample is sample's class
    [~, index] = min(distances(ii,:));
    sample_classes(ii) = classes(index);

    % Print the location and class of the current sample
    fprintf('The sample at (%.2f, %.2f, %.2f, %.2f) is in the %-15s class\n',...
            S(ii, 1), S(ii, 2), S(ii, 3), S(ii, 4), sample_classes(ii));
end

% Calculates the distance of each particle from each data point in 4D space
function dist = fit_func(particles, datapoints, num_datapoints, val)
    % Initialize the distance matrixa
    dist = zeros(size(particles, 1), num_datapoints);

    % For each particle
    for ii = 1:size(particles, 1)
        % For each datapoint (1-50: setosa, 51-100: versicolor, 101-150: virginica)
        % eg, virginica is 50*(2-1)+1:50*2, giving us the range 51-100
        for jj = num_datapoints * (val - 1) + 1 : num_datapoints * val
            % Calculate the euclidean distance from particle to datapoint
            dist(ii, jj) = sqrt((particles(ii,1) - datapoints(jj,1)).^2 ...
                              + (particles(ii,2) - datapoints(jj,2)).^2 ...
                              + (particles(ii,3) - datapoints(jj,3)).^2 ...
                              + (particles(ii,4) - datapoints(jj,4)).^2 );
        end
    end
end

% Sets the starting position for the particles
function [w, x, y, z] = set_starting_pos(sepal_length, sepal_width,...
                        petal_length, petal_width, num_particles)
w = rand(num_particles, 1) .* max(sepal_length); % Sepal length axis
x = rand(num_particles, 1) .* max(sepal_width);  % Sepal width axis
y = rand(num_particles, 1) .* max(petal_length); % Petal length axis
z = rand(num_particles, 1) .* max(petal_width);  % Petal width axis
end

% Sets the starting velocity for the particles
function [w, x, y, z] = set_starting_vel(num_particles)
w = randn(num_particles, 1) .* 0.5; % Sepal length axis
x = randn(num_particles, 1) .* 0.5; % Sepal width axis
y = randn(num_particles, 1) .* 0.5; % Petal length axis
z = randn(num_particles, 1) .* 0.5; % Petal width axis
end
