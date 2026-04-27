% Closed-form symbolic path for raipal CFB methods
% abandoned due to jumps in b_table near a = 0

clc ;
clear;
close all;
%%
syms a a1 a2 b b1
% a: relevant ange
% a1: angle from alignment with a0
% a2: angle from lower singularity

l0 = 112.95;
l1 = 60;
l2 = 45;
l3 = 85.47477864;

a0 = -1.28905643  * pi / 180; % zero position for a
b0 =  57.28905643 * pi / 180; % zero position for b
c0 =  156 * pi / 180; % zero position for b (generally not needed, but calculated for now)

% allowed range of motion for a: 0 - 89.13041976
% sigularity: a = -12.64575024 / a = +94.10537442

% table parameters
% a_table = linspace(-10,92,4000)' * pi / 180;
a_table = linspace(-5,94,4000)' * pi / 180;

% Hard-coded Ground Truth
GT_a = [0 10 20 30 40 50 60 70 80 89.13041976]' * pi/180;
GT_b = [0 14.03243845 32.38612781 50.50474484 66.47558509 80.26137339 92.48673122 103.98822317 115.96641313 130]' * pi/180;

%% Symbolic solution

% full analytic
equation_b = (l0 + l2*cos(b + b0) - l1*cos(a + a0))^2 + (l1*sin(a + a0) + l2*sin(b + b0))^2 == l3^2;
solution_b = solve(equation_b,b);
solution_b = solution_b(2);

% equation_b1 = (l0 + l2*cos(b1) - l1*cos(a1))^2 + (l1*sin(a1) + l2*sin(b1))^2 == l3^2;
% solution_b1 = solve(equation_b1,b1);
% solution_b = solution_b1(2) - b0;

solution_db = simplify(diff(solution_b,a1));

disp(solution_b)

%% test solution with GT
b_test = unwrap(real(eval(subs(solution_b,a,GT_a))));
% b_test = unwrap(real(eval(subs(solution_b,a1,GT_a + a0))));

vals = double([
    GT_a(:).' * 180/pi
    b_test(:).' * 180/pi
    (b_test(:).' - GT_b(:).') * 180/pi
]);

fprintf("a = %.8f deg --> b = %.8f deg (error: %.5f deg)\n", vals);

fprintf("Average error with GT: %.5e rad\n\n",mean(abs(b_test - GT_b)))

%% generate table

b_table = unwrap(real(eval(subs(solution_b,a,a_table))));
% b_table = unwrap(real(eval(subs(solution_b,a1,a_table + a0))));
c_table = c0 - (pi/2 - (a_table + a0) + acos( ...
    (l1 .* sin(a_table + a0) + l2 .* sin(b_table + b0)) ./ l3 ...
));

% db_table = unwrap(real(eval(subs(solution_db,a,a_table))));
% db_table_legacy = diff(b_table) ./ diff(a_table);

fprintf('table calculated\n\n')

%%
save('results/cfb_table.mat', 'a_table', 'b_table', 'c_table')