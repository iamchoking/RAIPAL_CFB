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

% Singularities: a = -12.64575024 deg, a = +94.10537442 deg
a_sing = [-12.64575024 * pi / 180 94.10537442 * pi / 180];

% Allowed range of motion for a: 0 to 89.13041976 deg

% Calculated table is wider including some margin around normal operating range
num_points = 4000;
a_table = linspace(-5, 94, num_points)' * pi / 180;

% Hard-coded Ground Truth
GT_a = [0 10 20 30 40 50 60 70 80 89.13041976]' * pi/180;
GT_b = [0 14.03243845 32.38612781 50.50474484 66.47558509 80.26137339 92.48673122 103.98822317 115.96641313 130]' * pi/180;

% draw plots at the end
make_plots = true;

%% Symbolic derivative setup

% We avoid using solve() for b_table because the symbolic inverse-trig
% expression can jump branches near a = 0.
%
% However, implicit symbolic differentiation is still useful and robust.

syms a_sym b_sym real

F = (l0 + l2*cos(b_sym + b0) - l1*cos(a_sym + a0))^2 + ...
    (l1*sin(a_sym + a0) + l2*sin(b_sym + b0))^2 - l3^2;

Fa  = diff(F, a_sym);
Fb  = diff(F, b_sym);

Faa = diff(F, a_sym, 2);
Fab = diff(diff(F, a_sym), b_sym);
Fbb = diff(F, b_sym, 2);

db_da_sym = simplify(-Fa / Fb);

ddb_da2_sym = simplify( ...
    -(Faa + 2*Fab*db_da_sym + Fbb*db_da_sym^2) / Fb ...
);

db_fun  = matlabFunction(db_da_sym,   'Vars', {a_sym, b_sym});
ddb_fun = matlabFunction(ddb_da2_sym, 'Vars', {a_sym, b_sym});

fprintf('Symbolic derivatives generated\n\n');

%% Generate b table using explicit two-branch geometry

% Define absolute angle A = a + a0 and B = b + b0.
%
% The linkage equation can be reduced to:
%
%     P*cos(B) + Q*sin(B) = K
%
% where:
%
%     P = l0 - l1*cos(A)
%     Q = l1*sin(A)
%
% The two branches are:
%
%     B = atan2(Q,P) +/- acos(K / hypot(P,Q))

A_table = a_table + a0;

P_table = l0 - l1*cos(A_table);
Q_table = l1*sin(A_table);
R_table = hypot(P_table, Q_table);

K_table = (l3^2 - l2^2 - P_table.^2 - Q_table.^2) ./ (2*l2);

acos_arg_table = K_table ./ R_table;

% Numerical safety for acos()
tol = 1e-10;

if any(abs(acos_arg_table) > 1 + tol)
    warning('Some b-table acos arguments are outside [-1, 1]. Check geometry/range.');
end

acos_arg_table = max(-1, min(1, acos_arg_table));

b_plus_table  = unwrap(atan2(Q_table, P_table) + acos(acos_arg_table) - b0);
b_minus_table = unwrap(atan2(Q_table, P_table) - acos(acos_arg_table) - b0);

%% Choose the physical b branch

% The desired physical branch has b near 0 when a is near 0.
[~, idx_zero_a] = min(abs(a_table));

if abs(b_plus_table(idx_zero_a)) <= abs(b_minus_table(idx_zero_a))
    b_table = b_plus_table;
    branch_name = 'plus';
else
    b_table = b_minus_table;
    branch_name = 'minus';
end

fprintf('Selected b branch: %s\n\n', branch_name);

%% Test b branch against ground truth

A_GT = GT_a + a0;

P_GT = l0 - l1*cos(A_GT);
Q_GT = l1*sin(A_GT);
R_GT = hypot(P_GT, Q_GT);

K_GT = (l3^2 - l2^2 - P_GT.^2 - Q_GT.^2) ./ (2*l2);

acos_arg_GT = K_GT ./ R_GT;

if any(abs(acos_arg_GT) > 1 + tol)
    warning('Some GT acos arguments are outside [-1, 1]. Check geometry/range.');
end

acos_arg_GT = max(-1, min(1, acos_arg_GT));

b_plus_GT  = unwrap(atan2(Q_GT, P_GT) + acos(acos_arg_GT) - b0);
b_minus_GT = unwrap(atan2(Q_GT, P_GT) - acos(acos_arg_GT) - b0);



if strcmp(branch_name, 'plus')
    b_test = b_plus_GT;
else
    b_test = b_minus_GT;
end

vals = [
    GT_a(:).'    * 180/pi
    b_test(:).'  * 180/pi
    (b_test(:).' - GT_b(:).') * 180/pi
];

fprintf("a = %.8f deg --> b = %.8f deg (error: %.5f deg)\n", vals);
fprintf("\nAverage error with GT: %.5e rad\n\n", mean(abs(b_test - GT_b)));

%% Generate c table

B_table = b_table + b0;

acos_arg_c = (l1 .* sin(A_table) + l2 .* sin(B_table)) ./ l3;

if any(abs(acos_arg_c) > 1 + tol)
    warning('Some c-table acos arguments are outside [-1, 1]. Check geometry/range.');
end

acos_arg_c = max(-1, min(1, acos_arg_c));

c_table = c0 - (pi/2 - A_table + acos(acos_arg_c));

%% Generate derivative tables

db_table = db_fun(a_table, b_table);
db_table = db_table(:);

ddb_table = ddb_fun(a_table, b_table);
ddb_table = ddb_table(:);

% Optional finite-difference checks
db_table_fd = gradient(b_table, a_table);
ddb_table_fd = gradient(db_table_fd, a_table);

% db_error_max = max(abs(db_table - db_table_fd));
% ddb_error_max = max(abs(ddb_table - ddb_table_fd));

% Find a value (in degrees) where maximum db and ddb error occur

% Restrict error search to allowed range of motion for a
a_min = 0 * pi/180;
a_max = 89.13041976 * pi/180;
allowed_idx = find(a_table >= a_min & a_table <= a_max);

[db_error_max, rel_idx_db_max] = max(abs(db_table(allowed_idx) - db_table_fd(allowed_idx)));
[ddb_error_max, rel_idx_ddb_max] = max(abs(ddb_table(allowed_idx) - ddb_table_fd(allowed_idx)));
a_db_max_deg = a_table(allowed_idx(rel_idx_db_max)) * 180/pi;
a_ddb_max_deg = a_table(allowed_idx(rel_idx_ddb_max)) * 180/pi;

fprintf('Derivative Error check (analytic vs finite diff, allowed a range):\n');
fprintf('  db/da    Maximum Error: %.5e at a = %.8f deg\n', db_error_max  , a_db_max_deg);
fprintf('  ddb/da2  Maximum Error: %.5e at a = %.8f deg\n', ddb_error_max , a_ddb_max_deg);

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
    'params', 'a_table', 'b_table', 'c_table', 'db_table', 'ddb_table');

fprintf('Saved table to %s\n', fullfile(results_dir, 'cfb_table.mat'));

%% Optional plots
a_table_deg  = a_table  * 180/pi;
b_table_deg  = b_table  * 180/pi;
c_table_deg  = c_table  * 180/pi;

if make_plots
    fig = figure('Name', 'CFB Tables and Derivatives', 'Position', [100, 100, 1200, 900]);

    subplot(2,2,1);
    plot(a_table_deg, b_table_deg, 'LineWidth', 1.5);
    grid on;
    xlabel('a [deg]');
    ylabel('b [deg]');
    title('b table');

    subplot(2,2,2);
    plot(a_table_deg, c_table_deg, 'LineWidth', 1.5);
    grid on;
    xlabel('a [deg]');
    ylabel('c [deg]');
    title('c table');

    subplot(2,2,3);
    plot(a_table_deg, db_table, 'LineWidth', 1.5);
    hold on;
    plot(a_table_deg, db_table_fd, '--');
    grid on;
    xlabel('a [deg]');
    ylabel('db/da');
    legend('analytic', 'finite difference');
    title('First derivative');

    subplot(2,2,4);
    plot(a_table_deg, ddb_table, 'LineWidth', 1.5);
    hold on;
    plot(a_table_deg, ddb_table_fd, '--');
    grid on;
    xlabel('a [deg]');
    ylabel('d^2b/da^2');
    legend('analytic', 'finite difference');
    title('Second derivative');

    % Save the figure as a PDF with the same size as the figure window
    pdf_path = fullfile(results_dir, 'cfb_table_plots.pdf');
    exportgraphics(fig, pdf_path, 'ContentType', 'vector');
    fprintf('Saved combined figure to %s\n', pdf_path);
end