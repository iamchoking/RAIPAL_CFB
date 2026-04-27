clc ;
clear;
close all;

syms alpha1 beta1 
syms a b1 cosb1 sinb1 cosa1 sina1
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
num_points = 4000;
a_table = linspace(-5,94,num_points)' * pi / 180;

%% Establish Ground Truth

GT_a = [0 10 20 30 40 50 60 70 80 89.13041976]' * pi/180;
GT_b = [0 14.03243845 32.38612781 50.50474484 66.47558509 80.26137339 92.48673122 103.98822317 115.96641313 130]' * pi/180;

%% generate symbolic solution

% full analytic
% equation_0 = (l0 + l2*cos(b1) - l1*cos(a + a0))^2 + (l1*sin(a + a0) + l2*sin(b1))^2 == l3^2;
% solution = solve(equation_0,b1);
% 
% disp(solution)

% analytic without sin / cos calculation
equation_1 = (l0 + l2*cosb1 - l1*cosa1)^2 + (l1*sina1 + l2*sinb1)^2 == l3^2;
equation_2 = sinb1^2 + cosb1^2 == 1;
solution = solve([equation_1,equation_2],[sinb1,cosb1]);

%% test solution with GT
cosa1_test = cos(GT_a + a0);
sina1_test = sin(GT_a + a0);

b_test = zeros(length(GT_a),1);
for i = 1:length(GT_a)

    sinb1_sol = subs(solution.sinb1(1),[sina1,cosa1],[sina1_test(i),cosa1_test(i)]);
    cosb1_sol = subs(solution.cosb1(1),[sina1,cosa1],[sina1_test(i),cosa1_test(i)]);
    b_test(i) = atan2(sinb1_sol,cosb1_sol) - b0;

    if b_test(i) < -0.1
        b_test(i) = b_test(i) + 2*pi;
    end
    
    fprintf("a = %.8f deg --> b = %.8f deg (error: %.5e deg)\n",GT_a(i) * 180/pi,b_test(i) * 180/pi,(b_test(i)-GT_b(i)) * 180/pi)

end

fprintf("Average error with GT: %.5e rad\n\n",mean(abs(b_test - GT_b)))

%% generate table
cosa1_table = cos(a_table + a0);
sina1_table = sin(a_table + a0);

b_table = zeros(length(a_table),1);
c_table = zeros(length(a_table),1);

for i = 1:length(a_table)

    sinb1_sol = subs(solution.sinb1(1),[sina1,cosa1],[sina1_table(i),cosa1_table(i)]);
    cosb1_sol = subs(solution.cosb1(1),[sina1,cosa1],[sina1_table(i),cosa1_table(i)]);

    b_table(i) = atan2(sinb1_sol,cosb1_sol) - b0;
    c_table(i) = c0 - (pi/2 - (a_table(i) + a0) + acos( (l1 * sina1_table(i) + l2 * double(sinb1_sol))/l3 ) );

    if b_table(i) < -0.1
        b_table(i) = b_table(i) + 2*pi;
    end

    if mod(i,num_points/10) == 0
        fprintf('Table Calculation: %d/%d complete\n',i,num_points)
    end

end

fprintf('table calculated\n\n')
