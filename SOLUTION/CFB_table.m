clc;
clear;
close all;

%% Parameters

l0 = 112.95;
l1 = 60;
l2 = 45;
l3 = 85.47477864;

a0 = -1.28905643  * pi / 180; % zero position for a
b0 =  57.28905643 * pi / 180; % zero position for b
c0 =  156         * pi / 180; % zero position for c

% Singularities
a_sing = [-12.64575024 94.10537442] * pi / 180;
b_sing = [-6.16480061 150.05312056] * pi / 180;

% Allowed range of motion for b: 0 to 130 deg

% Calculated table is wider including some margin around normal operating range
num_points = 4096;
b_table = linspace(-5, 135, num_points)' * pi / 180;
fprintf("Table step: %e deg\n",(b_table(2) - b_table(1)) * 180 / pi)

% Hard-coded Ground Truth
GT_a = [0 10 20 30 40 50 60 70 80 89.13041976]' * pi/180;
GT_b = [0 14.03243845 32.38612781 50.50474484 66.47558509 80.26137339 92.48673122 103.98822317 115.96641313 130]' * pi/180;

% draw plots at the end
make_plots = true;

%% Symbolic derivative setup

% Use implicit differentiation of F(a,b) = 0 for derivative tables.

syms a_sym b_sym real

F = (l0 + l2*cos(b_sym + b0) - l1*cos(a_sym + a0))^2 + ...
    (l1*sin(a_sym + a0) + l2*sin(b_sym + b0))^2 - l3^2;

Fa  = diff(F, a_sym);
Fb  = diff(F, b_sym);

Faa = diff(F, a_sym, 2);
Fab = diff(diff(F, a_sym), b_sym);
Fbb = diff(F, b_sym, 2);

% Derivatives: da/db and d2a/db2
da_db_sym = simplify(-Fb / Fa);
dda_db2_sym = simplify( ...
    -(Fbb + 2*Fab*da_db_sym + Faa*da_db_sym^2) / Fa ...
);

da_fun  = matlabFunction(da_db_sym,   'Vars', {a_sym, b_sym});
dda_fun = matlabFunction(dda_db2_sym, 'Vars', {a_sym, b_sym});

fprintf('Symbolic derivatives generated\n\n');

%% Generate a table from prescribed b table using explicit two-branch geometry

% Define absolute angles A = a + a0 and B = b + b0.
%
% The linkage equation can be reduced to:
%
%     P*cos(A) + Q*sin(A) = K
%
% where:
%
%     P = l0 + l2*cos(B)
%     Q = -l2*sin(B)
%
% The two branches are:
%
%     A = atan2(Q,P) +/- acos(K / hypot(P,Q))

B_table = b_table + b0;

P_table = l0 + l2*cos(B_table);
Q_table = -l2*sin(B_table);
R_table = hypot(P_table, Q_table);

K_table = (P_table.^2 + (l2*sin(B_table)).^2 + l1^2 - l3^2) ./ (2*l1);

acos_arg_table = K_table ./ R_table;

% Numerical safety for acos()
tol = 1e-10;

if any(abs(acos_arg_table) > 1 + tol)
    warning('Some a-table acos arguments are outside [-1, 1]. Check geometry/range.');
end

acos_arg_table = max(-1, min(1, acos_arg_table));

a_plus_table  = unwrap(atan2(Q_table, P_table) + acos(acos_arg_table) - a0);
a_minus_table = unwrap(atan2(Q_table, P_table) - acos(acos_arg_table) - a0);

%% Choose the physical a branch

% The desired physical branch has a near 0 when b is near 0.
[~, idx_zero_b] = min(abs(b_table));

if abs(a_plus_table(idx_zero_b)) <= abs(a_minus_table(idx_zero_b))
    a_table = a_plus_table;
    branch_name = 'plus';
else
    a_table = a_minus_table;
    branch_name = 'minus';
end

fprintf('Selected a branch: %s\n\n', branch_name);

%% Test a branch against ground truth

B_GT = GT_b + b0;

P_GT = l0 + l2*cos(B_GT);
Q_GT = -l2*sin(B_GT);
R_GT = hypot(P_GT, Q_GT);

K_GT = (P_GT.^2 + (l2*sin(B_GT)).^2 + l1^2 - l3^2) ./ (2*l1);

acos_arg_GT = K_GT ./ R_GT;

if any(abs(acos_arg_GT) > 1 + tol)
    warning('Some GT acos arguments are outside [-1, 1]. Check geometry/range.');
end

acos_arg_GT = max(-1, min(1, acos_arg_GT));

a_plus_GT  = unwrap(atan2(Q_GT, P_GT) + acos(acos_arg_GT) - a0);
a_minus_GT = unwrap(atan2(Q_GT, P_GT) - acos(acos_arg_GT) - a0);

if strcmp(branch_name, 'plus')
    a_test = a_plus_GT;
else
    a_test = a_minus_GT;
end

vals = [
    GT_b(:).'    * 180/pi
    a_test(:).'  * 180/pi
    (a_test(:).' - GT_a(:).') * 180/pi
];

fprintf("b = %.8f deg --> a = %.8f deg (error: %.3e deg)\n", vals);
fprintf("\nAverage error with GT: %.5e deg\n\n", mean(abs(a_test - GT_a)) * 180/pi);

%% Generate c table

A_table = a_table + a0;

acos_arg_c = (l1 .* sin(A_table) + l2 .* sin(B_table)) ./ l3;

if any(abs(acos_arg_c) > 1 + tol)
    warning('Some c-table acos arguments are outside [-1, 1]. Check geometry/range.');
end

acos_arg_c = max(-1, min(1, acos_arg_c));

c_table = c0 - (pi/2 - A_table + acos(acos_arg_c));

%% Generate derivative tables

da_table = da_fun(a_table, b_table);
da_table = da_table(:);

dda_table = dda_fun(a_table, b_table);
dda_table = dda_table(:);

% Optional finite-difference checks
da_table_fd = gradient(a_table, b_table);
dda_table_fd = gradient(da_table_fd, b_table);

% Restrict error search to allowed range of motion for b
b_min = 0 * pi/180;
b_max = 130 * pi/180;
allowed_idx = find(b_table >= b_min & b_table <= b_max);

[da_error_max, rel_idx_da_max] = max(abs(da_table(allowed_idx) - da_table_fd(allowed_idx)));
[dda_error_max, rel_idx_dda_max] = max(abs(dda_table(allowed_idx) - dda_table_fd(allowed_idx)));
b_da_max_deg = b_table(allowed_idx(rel_idx_da_max)) * 180/pi;
b_dda_max_deg = b_table(allowed_idx(rel_idx_dda_max)) * 180/pi;

fprintf('Derivative Error check (analytic vs finite diff, allowed b range):\n');
fprintf('  da/db    Maximum Error: %.5e at b = %.8f deg\n', da_error_max  , b_da_max_deg);
fprintf('  dda/db2  Maximum Error: %.5e at b = %.8f deg\n', dda_error_max , b_dda_max_deg);

fprintf('Table calculated\n\n');

%% Save results
params = struct();

params.l0 = l0;
params.l1 = l1;
params.l2 = l2;
params.l3 = l3;

params.a0 = a0;
params.b0 = b0;
params.c0 = c0;

params.a0_deg = a0 * 180/pi;
params.b0_deg = b0 * 180/pi;
params.c0_deg = c0 * 180/pi;

params.a_sing = a_sing;

params.num_points = num_points;
params.branch_name = branch_name;

results_dir = 'results';

if ~exist(results_dir, 'dir')
    mkdir(results_dir);
end

save(fullfile(results_dir, 'cfb_table.mat'), ...
    'params', 'a_table', 'b_table', 'c_table', 'da_table', 'dda_table');

fprintf('Saved table to %s\n', fullfile(results_dir, 'cfb_table.mat'));

%% Optional plots
a_table_deg  = a_table  * 180/pi;
b_table_deg  = b_table  * 180/pi;
c_table_deg  = c_table  * 180/pi;

if make_plots
    fig = figure('Name', 'CFB Tables and Derivatives', 'Position', [100, 100, 1200, 900]);

    subplot(2,2,1);
    plot(b_table_deg, a_table_deg, 'LineWidth', 1.5);
    grid on;
    xlabel('b [deg]');
    ylabel('a [deg]');
    title('a table');

    subplot(2,2,2);
    plot(b_table_deg, c_table_deg, 'LineWidth', 1.5);
    grid on;
    xlabel('b [deg]');
    ylabel('c [deg]');
    title('c table');

    subplot(2,2,3);
    plot(b_table_deg, da_table, 'LineWidth', 1.5);
    hold on;
    plot(b_table_deg, da_table_fd, '--');
    grid on;
    xlabel('b [deg]');
    ylabel('da/db');
    legend('analytic', 'finite difference');
    title('First derivative');

    subplot(2,2,4);
    plot(b_table_deg, dda_table, 'LineWidth', 1.5);
    hold on;
    plot(b_table_deg, dda_table_fd, '--');
    grid on;
    xlabel('b [deg]');
    ylabel('d^2a/db^2');
    legend('analytic', 'finite difference');
    title('Second derivative');

    % Save the figure as a PDF with the same size as the figure window
    pdf_path = fullfile(results_dir, 'cfb_solution_plots.pdf');
    exportgraphics(fig, pdf_path, 'ContentType', 'vector');
    fprintf('Saved combined figure to %s\n', pdf_path);
end
