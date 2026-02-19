clc ;
clear;
close all;

syms alpha1 beta1 
syms a b1 cosb1 sinb1 cosa1 sina1
l0 = 112.95;
l1 = 60;
l2 = 45;
l3 = 85.47477864;

a0 = -1.28905643  * pi / 180;
b0 =  57.28905643 * pi / 180;
c0 =  156 * pi / 180;
%% Establish GT
GT_a = [0 10 20 30 40 50 60 70 80 89.13041976]' * pi/180;
GT_b = [0 14.03243845 32.38612781 50.50474484 66.47558509 80.26137339 92.48673122 103.98822317 115.96641313 130]' * pi/180;

%% generate symbolic solution
clc 

equation_1 = (l0 + l2*cosb1 - l1*cosa1)^2 + (l1*sina1 + l2*sinb1)^2 == l3^2;
equation_2 = sinb1^2 + cosb1^2 == 1;
% equation_0 = (l0 + l2*cos(b1) - l1*cos(a + a0))^2 + (l1*sin(a + a0) + l2*sin(b1))^2 == l3^2;
solution = solve([equation_1,equation_2],[sinb1,cosb1]);

%% test solution with GT
cosa1_numerical = cos(GT_a + a0);
sina1_numerical = sin(GT_a + a0);

b_numerical = zeros(length(GT_a),1);
for i = 1:length(GT_a)

    sinb1_sol = subs(solution.sinb1(1),[sina1,cosa1],[sina1_numerical(i),cosa1_numerical(i)]);
    cosb1_sol = subs(solution.cosb1(1),[sina1,cosa1],[sina1_numerical(i),cosa1_numerical(i)]);
    b_numerical(i) = atan2(sinb1_sol,cosb1_sol) - b0;

    if b_numerical(i) < -0.1
        b_numerical(i) = b_numerical(i) + 2*pi;
    end
    
    fprintf("a = %.8f deg --> b = %.8f deg (error: %.5e deg)\n",GT_a(i) * 180/pi,b_numerical(i) * 180/pi,(b_numerical(i)-GT_b(i)) * 180/pi)

end

fprintf("Average error with GT: %.5e rad\n\n",mean(abs(b_numerical - GT_b)))

%% generate table of solutions
num_points = 4000;

a_table = linspace(-5,94,num_points)' * pi / 180;
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

%%
a_domain = linspace(0,89.13041976,num_points)' * pi / 180 ;
cosa1_domain = cos(a_domain + a0);
sina1_domain = sin(a_domain + a0);

b_domain = zeros(length(a_domain),1);
c_domain = zeros(length(a_domain),1);

for i = 1:length(a_domain)

    sinb1_sol = subs(solution.sinb1(1),[sina1,cosa1],[sina1_domain(i),cosa1_domain(i)]);
    cosb1_sol = subs(solution.cosb1(1),[sina1,cosa1],[sina1_domain(i),cosa1_domain(i)]);

    b_domain(i) = atan2(sinb1_sol,cosb1_sol) - b0;
    c_domain(i) = c0 - (pi/2 - (a_domain(i) + a0) + acos( (l1 * sina1_domain(i) + l2 * double(sinb1_sol))/l3 ) );

    if b_domain(i) < -0.1
        b_domain(i) = b_domain(i) + 2*pi;
    end

    if mod(i,num_points/10) == 0
        fprintf('Domain Calculation: %d/%d complete\n',i,num_points)
    end

end

fprintf('domain calculated\n\n')

%% preprocess dataset for fitting
n = 16; %polynomial degree
a = preprocess(a_domain,num_points);

%% Fit polynomial for b

[b_poly,b_fit] = polyfit_analyze(a_domain,b_domain,n,"b",1,1.0e-3/(180/pi));

% add GT plot
figure(1);
hold on;
scatter(GT_a,GT_b,"DisplayName","Ground Truth")
hold off;

% vanilla 16 degree
% Polynomial coefficients (highest to lowest degree):
% [5.267736291, -61.50765323, 325.9971618, -1038.387936, 2219.603936, -3370.197818, 3757.847732, -3141.343675, 1991.046707, -952.8842762, 329.2544421, -69.63848238, 4.694541019, -1.732916861, 2.516434102, 1.022471236, 0.000004678279257]
%  
% Mean Square error : 2.140311e-12 rad^2
% Mean difference   : 6.244870e-05 deg
% Max  difference   : 8.861856e-04 deg

% caveman weighting 16 degree
% Polynomial coefficients (highest to lowest degree):
% [7.39492297, -87.32937574, 468.0311062, -1506.68021, 3251.217616, -4972.496218, 5562.446965, -4634.134708, 2899.71248, -1356.776534, 458.0671211, -98.25606708, 8.922560669, -2.119034782, 2.535726081, 1.022052864, 0.000006869097201]
%  
% Mean Square error : 3.332970e-12 rad^2
% Mean difference   : 8.718492e-05 deg
% Max  difference   : 3.935761e-04 deg

% hand-tuned weighting 16 degree
% Polynomial coefficients (highest to lowest degree):
% [7.743760988519409, -91.21738209662894, 487.0807418051794, -1560.06991422005, 3343.487739917685, -5067.667199180803, 5603.110257558747, -4599.901958118891, 2827.536921428657, -1296.048945332545, 427.0774524818472, -88.03097782563684, 6.759940197424146, -1.840913614666093, 2.516241535884542, 1.022642767854859, 0.000002946273850102822]
%  
% Mean Square error : 3.940240e-12 rad^2
% Mean difference   : 9.952158e-05 deg
% Max  difference   : 2.824994e-04 deg

%% Fit polynomial for c

[c_poly,c_fit] = polyfit_analyze(a_domain,c_domain,n,"c",2,1.0e-3/(180/pi));

% hand-tuned weighting
% Polynomial coefficients (highest to lowest degree):
% [-4.599352352699746, 54.74231689702444, -295.2311950384566, 953.4731854791298, -2052.475077842654, 3099.83081925634, -3363.608943070999, 2636.542385991094, -1477.728449862623, 579.0562233124353, -157.9024440428828, 35.05094543077983, -7.773109367008759, -0.656408504334005, 0.798469056970581, 2.098468226080274, -0.000001063969759856979]
%  
% Mean Square error : 8.267375e-13 rad^2
% Mean difference   : 4.477672e-05 deg
% Max  difference   : 1.460069e-04 deg


%% differential kinematics

% numerically differentiate for angular speed
wb_domain = diff(b_domain) ./ diff(a_domain);
wc_domain = diff(c_domain) ./ diff(a_domain);

% fill in first value assuming constant acceleration
wb_domain = [wb_domain(1) - (wb_domain(2) - wb_domain(1));wb_domain];
wc_domain = [wc_domain(1) - (wc_domain(2) - wc_domain(1));wc_domain];

%% Fit polynomial for wb
[wb_poly,wb_fit] = polyfit_analyze(a_domain,wb_domain,n,"\omega_b / \omega_a",3,1.0e-3);
% [wb_poly,wb_fit] = polyfit_analyze(a_domain,wb_domain,n,"\omega_a / \omega_b",3,1.0e-3);

% [\omega_b / \omega_a] Polynomial coefficients (highest to lowest degree):
% [330.1514673459376, -3939.920472235199, 21319.39656007634, -69182.34049371608, 150034.7654546927, -229354.9209106861, 253935.6748254774, -205870.6560746802, 121976.630364867, -52021.60108263157, 15489.35547273981, -3096.107249114025, 433.0072547668175, -57.21741621150876, -0.9670218724995983, 4.921276349693374, 1.022525952484172]
%  
% [\omega_b / \omega_a] Mean Square error : 5.727483e-09
% [\omega_b / \omega_a] Mean difference   : 6.507945e-05
% [\omega_b / \omega_a] Max  difference   : 2.213536e-04

%% Fit polynomial for wc
[wc_poly,wc_fit] = polyfit_analyze(a_domain,wc_domain,n,"\omega_c / \omega_a",4,1.0e-3);

% [\omega_c / \omega_a] Polynomial coefficients (highest to lowest degree):
% [-160.5479994847285, 1895.920138547342, -10134.44134516569, 32423.08074300719, -69176.96131680353, 103833.9828875482, -112772.7868033985, 89865.6754559566, -52817.85947080839, 22895.36338178862, -7274.853311104497, 1636.804103542675, -218.7685080450584, 6.311919990384402, -3.93273889555925, 1.64496743133799, 2.097821834064398]
%  
% [\omega_c / \omega_a] Mean Square error : 1.713449e-09
% [\omega_c / \omega_a] Mean difference   : 3.575284e-05
% [\omega_c / \omega_a] Max  difference   : 1.188441e-04

%% inverse kinematics

% for inverse b (b->a)
[invb_poly,invb_fit] = polyfit_analyze(b_domain,a_domain,n,"inverse-b", 5,1.0e-3/(180/pi));

% [inverse-b] Polynomial coefficients (highest to lowest degree):
% [-0.01820207462949832, 0.3515207083714699, -3.095862950743284, 16.46669621013979, -59.05557308462212, 151.014385452334, -284.0127408330851, 399.5792368412747, -423.9498339027933, 339.9044621761164, -205.6485389436184, 93.83159046256269, -32.59323031770014, 9.041352427928961, -2.195749936345363, 0.9751610264689757, 0.000008330739145766052]
%  
% [inverse-b] Mean Square error : 4.956606e-12
% [inverse-b] Mean difference   : 1.688944e-06
% [inverse-b] Max  difference   : 8.457033e-06

% for inverse c (b->a)
% this is a bad idea since the solution is not unique!

%% Plotting effective reduction ratio
set(groot, "defaultAxesFontSize",   25);
set(groot, "defaultTextFontSize",   30);
set(groot, "defaultLegendFontSize", 20);
set(groot, "defaultColorbarFontSize", 20);
set(groot, 'DefaultLineLineWidth', 3);

if ishandle(99)
    close(99);
end
fig_reduction = figure(99);
fig_reduction.Position = [100 100 960 500];
inv_wb_fit = 1 ./ wb_fit;
normalized_inv_wb_fit = inv_wb_fit / max(inv_wb_fit(300:num_points,:));

title("Elbow Effective Reduction Ratio")
xlabel("Elbow Angle [deg]")
ylabel("Effective Gear Ratio")
% xlim([0 max(b_fit * 180 / pi)]);
% xlim([15 max(b_fit * 180 / pi)]);
xlim([15 120])
ylim([0.5 0.9]);
hold on;
legend('location','southeast');
% plot(b_human, inv_wb_human, "-o", "DisplayName", "Human Elbow Joint");
% exportgraphics(fig_reduction,"results/FIGURE_reduction_human.png","Resolution",300)

plot(b_fit * 180 / pi, inv_wb_fit,"DisplayName","RAIPAL Crossed 4-Bar");

hold off;

exportgraphics(fig_reduction,"results/FIGURE_reduction_full.png","Resolution",300)


%% function definitions
function [x_poly,x_fit] = polyfit_analyze(a_domain,x_domain,n,name,figure_num,err_factor)

    num_points = length(a_domain);

    a = preprocess(a_domain,num_points);
    x = preprocess(x_domain,num_points);
    
    x_poly = polyfit(a,x,n);
    % p(17) = 0;
    
    fprintf('\n\n[%s] Polynomial coefficients (highest to lowest degree):\n', name);
    disp(vpa(x_poly,16));
    
    x_fit = polyval(x_poly,a_domain);
    
    % Calculate Error
    fit_err    = x_domain - x_fit;
    L2_error   = norm(fit_err)^2 / num_points;
    mean_error = mean(abs(fit_err));
    max_error  = max(abs(fit_err));
    
    fprintf('[%s] Mean Square error : %.6e\n',name, L2_error);
    fprintf('[%s] Mean difference   : %.6e\n',name, mean_error);
    fprintf('[%s] Max  difference   : %.6e\n',name, max_error);
    
    if ishandle(figure_num)
        close(figure_num);
    end
    fig = figure(figure_num);
    title(sprintf("RAIPAL Crossed 4-bar Polynomial Analysis: %s",name))
    xlabel("a(input) [rad]")
    ylabel(name)
    xlim([0 max(a_domain)]);
    % ylim([0 GT_b(10)]);
    hold on;
    legend;
    plot(a_domain,x_domain,"DisplayName","Numerical Solution")
    plot(a_domain,abs(fit_err)/err_factor,"DisplayName",sprintf("Approximation Err. (/%.2e)",err_factor))
    hold off;

    exportgraphics(fig,sprintf("results/FIGURE-%d_%s--a.png",figure_num, erase(strrep(name,"/","-"),"\")),"Resolution",300)

end

function arr = preprocess(arr, num_points)
    % AUGMENT_DATA Augment array by appending repeated slices
    %   arr = augment_data(arr, num_points) appends various slices of the
    %   input array to itself multiple times for data augmentation
    %
    %   Inputs:
    %       arr - Input array (column vector or matrix)
    %       num_points - Original number of data points for indexing
    %
    %   Output:
    %       arr - Augmented array
    
    % Append last 10% of data, 4 times
    arr = [arr; repmat(arr(round(num_points * 0.9):num_points, :), 4, 1)];
    
    % Append first 50 points, 10 times
    arr = [arr; repmat(arr(1:50, :), 10, 1)];
    
    % Append last 0.5% of data, 15 times
    arr = [arr; repmat(arr(round(num_points * 0.995):num_points, :), 15, 1)];
    
    % Append last 11 points, 40 times
    arr = [arr; repmat(arr(num_points - 10:num_points, :), 40, 1)];
end