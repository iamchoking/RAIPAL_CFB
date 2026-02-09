clc 
clear

syms alpha1 beta1 
syms cosa1 sina1 cosb1 sinb1
l0 = 112.95;
l1 = 60;
l2 = 45;
l3 = 85.47477864;

%%
equation      = (l0 + l2*cos(beta1) - l1*cos(alpha1))^2 + (l1*sin(alpha1) + l2*sin(beta1))^2 == l3^2;


solutions = solve(equation, beta1);
fprintf('\nAnalytical solutions:\n');
for i = 1:length(solutions)
    % Try multiple simplification strategies
    sol = solutions(i);
%     sol = rewrite(sol, 'tan');

    sol = vpa(sol,10);
    sol = simplify(sol, 'Steps', 100);
    sol = vpa(sol,8);

    fprintf('\nbeta1_%d = \n', i);
    fprintf('Compact form: %s\n', char(sol));

    fprintf('Substituting pi/4 (compact): %.10f\n\n', subs(sol         ,alpha1,pi/4))

end

%%
equation_alt1 = (l0 + l2*cosb1       - l1*cosa1)^2 + (l1*sina1 + l2*sinb1      )^2 == l3^2;
equation_alt2 = sinb1^2 + cosb1^2 == 1;

solutions = solve([equation_alt1,equation_alt2], [sinb1,cosb1]);
fprintf('\n%d Analytical solutions:\n',length(solutions));
for i = 1:length(solutions)
    % Try multiple simplification strategies
    sol = solutions(i);

    sol_cos = simplify(sol.cosb1(1), 'Steps', 100);
    sol_cos = vpa(sol_cos,8);

    sol_sin = simplify(sol.sinb1(1), 'Steps', 100);
    sol_sin = vpa(sol_sin,8);

    sol = simplify(atan2(sol_sin,sol_cos), 'Steps',100);
    sol = vpa(sol,8);

    fprintf('cos(beta1): %s\n', char(sol_cos));
    fprintf('sin(beta1): %s\n', char(sol_sin));
    fprintf('beta1: %s\n', char(sol));

    % limit 1: 0 (should yield beta = 0)
    alpha_0 = -1.28905643 * pi / 180;
    beta_0 = 57.28905643 * pi / 180;

    fprintf('Substituting Limit 1: ')
    alpha = 0;
    beta = subs(sol,[sina1,cosa1],[sin(alpha+alpha_0),cos(alpha+alpha_0)]) - beta_0;
    fprintf('%.10f deg\n\n', beta * 180 / pi)

    % limit 2: 1.556 (should yield beta = 130deg = 2.26893 rad)
    fprintf('Substituting Limit 2: ')
    alpha = 89.13041976 * pi / 180;
    beta = subs(sol,[sina1,cosa1],[sin(alpha+alpha_0),cos(alpha+alpha_0)]) - beta_0;
    fprintf('%.10f deg\n\n', beta * 180 / pi)
    
end
%%
alpha = 89.13041976 * pi / 180;
x = 180 / pi * double(subs(-beta_0 + atan2(- (1.3473996e-17*(4.9478023e+16*cosa1^2 - 1.8628476e+17*cosa1 + 4.9478023e+16*sina1^2 + 1.0275987e+17))/sina1 - (1.3473996e-17*(400.0*cosa1 - 753.0)*(9.2771294e+13*sina1*(3.2620769e+11*cosa1 + 3.4269867e+11*cosa1*sina1^2 - 7.3177239e+11*cosa1^2 + 3.4269867e+11*cosa1^3 - 4.5511111e+10*cosa1^4 - 8.664215e+10*sina1^2 - 4.5511111e+10*sina1^4 - 9.1022222e+10*cosa1^2*sina1^2 + 1.6657692e+11)^(1/2) - 1.8137637e+20*cosa1 - 1.9791209e+19*cosa1*sina1^2 + 1.1177085e+20*cosa1^2 - 1.9791209e+19*cosa1^3 + 3.7256952e+19*sina1^2 + 7.7378183e+19))/(sina1*(160000.0*cosa1^2 - 602400.0*cosa1 + 160000.0*sina1^2 + 567009.0)), 0.66666667*cosa1 - (1.0*(156472.23*cosa1 + 0.5*sina1*(3.2620769e+11*cosa1 + 3.4269867e+11*cosa1*sina1^2 - 7.3177239e+11*cosa1^2 + 3.4269867e+11*cosa1^3 - 4.5511111e+10*cosa1^4 - 8.664215e+10*sina1^2 - 4.5511111e+10*sina1^4 - 9.1022222e+10*cosa1^2*sina1^2 + 1.6657692e+11)^(1/2) - 294558.97))/(160000.0*cosa1^2 - 602400.0*cosa1 + 160000.0*sina1^2 + 567009.0) - 1.255), [sina1,cosa1], [sin(alpha+alpha_0),cos(alpha+alpha_0)]));
vpa(x,30)

alpha = 0;
x = 180 / pi * double(subs(-beta_0 + atan2(- (1.3473996e-17*(4.9478023e+16*cosa1^2 - 1.8628476e+17*cosa1 + 4.9478023e+16*sina1^2 + 1.0275987e+17))/sina1 - (1.3473996e-17*(400.0*cosa1 - 753.0)*(9.2771294e+13*sina1*(3.2620769e+11*cosa1 + 3.4269867e+11*cosa1*sina1^2 - 7.3177239e+11*cosa1^2 + 3.4269867e+11*cosa1^3 - 4.5511111e+10*cosa1^4 - 8.664215e+10*sina1^2 - 4.5511111e+10*sina1^4 - 9.1022222e+10*cosa1^2*sina1^2 + 1.6657692e+11)^(1/2) - 1.8137637e+20*cosa1 - 1.9791209e+19*cosa1*sina1^2 + 1.1177085e+20*cosa1^2 - 1.9791209e+19*cosa1^3 + 3.7256952e+19*sina1^2 + 7.7378183e+19))/(sina1*(160000.0*cosa1^2 - 602400.0*cosa1 + 160000.0*sina1^2 + 567009.0)), 0.66666667*cosa1 - (1.0*(156472.23*cosa1 + 0.5*sina1*(3.2620769e+11*cosa1 + 3.4269867e+11*cosa1*sina1^2 - 7.3177239e+11*cosa1^2 + 3.4269867e+11*cosa1^3 - 4.5511111e+10*cosa1^4 - 8.664215e+10*sina1^2 - 4.5511111e+10*sina1^4 - 9.1022222e+10*cosa1^2*sina1^2 + 1.6657692e+11)^(1/2) - 294558.97))/(160000.0*cosa1^2 - 602400.0*cosa1 + 160000.0*sina1^2 + 567009.0) - 1.255), [sina1,cosa1], [sin(alpha+alpha_0),cos(alpha+alpha_0)]));
vpa(x,30)

%%
GT_alpha = [0 10 20 30 40 50 60 70 80 89.13041976] * pi/180;
GT_beta = [0 14.03243845 32.38612781 50.50474484 66.47558509 80.26137339 92.48673122 103.98822317 115.96641313 130] * pi/180;
%%
num_points = 3000;
alpha_range = linspace(0,89.13041976 * pi / 180,num_points) ;
cosa1_range = cos(alpha_range +  alpha_0);
sina1_range = sin(alpha_range +  alpha_0);
beta_range = linspace(0,89.13041976 * pi / 180,num_points);

for i = 1:length(alpha_range)
    beta_range(i) = double(subs(-beta_0 + atan2(- (1.3473996e-17*(4.9478023e+16*cosa1^2 - 1.8628476e+17*cosa1 + 4.9478023e+16*sina1^2 + 1.0275987e+17))/sina1 - (1.3473996e-17*(400.0*cosa1 - 753.0)*(9.2771294e+13*sina1*(3.2620769e+11*cosa1 + 3.4269867e+11*cosa1*sina1^2 - 7.3177239e+11*cosa1^2 + 3.4269867e+11*cosa1^3 - 4.5511111e+10*cosa1^4 - 8.664215e+10*sina1^2 - 4.5511111e+10*sina1^4 - 9.1022222e+10*cosa1^2*sina1^2 + 1.6657692e+11)^(1/2) - 1.8137637e+20*cosa1 - 1.9791209e+19*cosa1*sina1^2 + 1.1177085e+20*cosa1^2 - 1.9791209e+19*cosa1^3 + 3.7256952e+19*sina1^2 + 7.7378183e+19))/(sina1*(160000.0*cosa1^2 - 602400.0*cosa1 + 160000.0*sina1^2 + 567009.0)), 0.66666667*cosa1 - (1.0*(156472.23*cosa1 + 0.5*sina1*(3.2620769e+11*cosa1 + 3.4269867e+11*cosa1*sina1^2 - 7.3177239e+11*cosa1^2 + 3.4269867e+11*cosa1^3 - 4.5511111e+10*cosa1^4 - 8.664215e+10*sina1^2 - 4.5511111e+10*sina1^4 - 9.1022222e+10*cosa1^2*sina1^2 + 1.6657692e+11)^(1/2) - 294558.97))/(160000.0*cosa1^2 - 602400.0*cosa1 + 160000.0*sina1^2 + 567009.0) - 1.255), [sina1,cosa1], [sina1_range(i),cosa1_range(i)]));
    if beta_range(i) < 0
        beta_range(i) = beta_range(i) + pi*2;
    end
    if mod(i,num_points/10) == 0
        fprintf('Table Calculation: %d/%d complete\n',i,num_points)
    end
end

fprintf('table calculated\n\n')

p = polyfit(alpha_range,beta_range,15);
beta_fitted = polyval(p,alpha_range);

%%
fprintf('Polynomial coefficients (highest to lowest degree):\n');
disp(vpa(p,10));

L2_error   = norm(beta_range - beta_fitted)^2 / num_points;
mean_error = mean(abs(beta_range - beta_fitted));
max_error  =  max(abs(beta_range - beta_fitted));

fprintf('Mean Square error: %.6e\n', L2_error);
fprintf('Mean difference: %.6e\n', mean_error);
fprintf('Max  difference: %.6e\n',  max_error);

%%
figure;
plot(alpha_range,beta_range)
hold on;
scatter(GT_alpha,GT_beta)
plot(alpha_range,beta_fitted)
hold off;

