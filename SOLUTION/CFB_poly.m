%%
clc;
clear;
close all;

%%
load("results/cfb_table.mat")

if ~exist("db_table", "var") || ~exist("ddb_table", "var")
    error("results/cfb_table.mat does not contain db_table/ddb_table. Run CFB_table.m first.");
end

%% Polynomial degree
n = 31; % polynomial degree
fit_mode = "least-squares"; 
% fit_mode = "minimax";

% "least-squares", "minimax", or "auto"

range = dictionary;
range("b") = {[0 130] / 180 * pi};
range("a") = {[0 89.13041976] / 180 * pi};

%% Fit bidirectional polynomials
coeffs_b = struct();
[coeffs_b.to_a   , ~] = poly_analyze([b_table a_table],   ["b" "a"]          , ["deg" "deg"     ], n, fit_mode, 1, range);
[coeffs_b.to_da  , ~] = poly_analyze([b_table da_table],  ["b" "da / db"]    , ["deg" "deg/deg" ], n, fit_mode, 2, range);
[coeffs_b.to_dda , ~] = poly_analyze([b_table dda_table], ["b" "d^2a / db^2"], ["deg" "deg^{-1}"], n, fit_mode, 3, range);

coeffs_a = struct();
[coeffs_a.to_b   , ~] = poly_analyze([a_table b_table],   ["a" "b"]          , ["deg" "deg"     ], n, fit_mode, 4, range);
[coeffs_a.to_db  , ~] = poly_analyze([a_table db_table],  ["a" "db / da"]    , ["deg" "deg/deg" ], n, fit_mode, 5, range);
[coeffs_a.to_ddb , ~] = poly_analyze([a_table ddb_table], ["a" "d^2b / da^2"], ["deg" "deg^{-1}"], n, fit_mode, 6, range);

save('results/cfb_coeffs.mat',"coeffs_b","coeffs_a")

%% function definitions
function [coeffs,output_fit] = poly_analyze(table,name,units,num_degrees,fit_mode,fig_num,range_map)

    input_table = table(:,1);
    output_table = table(:,2);
    input_name = name(1);
    output_name = name(2);
    input_units = units(1);
    output_units = units(2);

    func_name = sprintf('%s-%s',output_name,input_name);
    num_points = length(input_table);

    factors = [1 1];
    factors(1) = unit_factor(input_units);
    factors(2) = unit_factor(output_units);

    x_min = min(input_table);
    x_max = max(input_table);
    x_norm = normalize_input(input_table, x_min, x_max);

    [p, ~, fit_method] = normalized_polyfit(x_norm, output_table, num_degrees, fit_mode);

    coeffs = struct();
    coeffs.method = fit_method;
    coeffs.min = x_min;
    coeffs.max = x_max;
    coeffs.degree = num_degrees;
    coeffs.coeffs = p;
    
%     fprintf('[%s] Normalization: x = 2 * (input - %.16g) / %.16g - 1\n', ...
%         func_name, x_min, x_max - x_min);
%     fprintf('\n\n[%s] Normalized %s polynomial coefficients (highest to lowest degree):\n', func_name, fit_method);
%     disp(p);

    output_fit = polyval(p,x_norm);
    input_plot = input_table * factors(1);
    output_fit_plot = output_fit * factors(2);
    output_table_plot = output_table * factors(2);

    % Calculate Error
    fit_err    = output_table_plot - output_fit_plot;
    L2_error   = norm(fit_err)^2 / num_points;
    mean_error = mean(abs(fit_err));
    max_error  = max(abs(fit_err));

    L2_error_suffix = "";
    mean_error_suffix = "";
    max_error_suffix = "";

    if has_range(range_map, input_name)
        input_range = range_map(string(input_name));
        input_range = input_range{1};
        in_range_idx = input_table >= input_range(1) & input_table <= input_range(2);

        if any(in_range_idx)
            fit_err_in_range = fit_err(in_range_idx);
            L2_error_in_range = norm(fit_err_in_range)^2 / nnz(in_range_idx);
            mean_error_in_range = mean(abs(fit_err_in_range));
            max_error_in_range = max(abs(fit_err_in_range));

            L2_error_suffix = sprintf(" (%.3e %s in range)", L2_error_in_range, output_units);
            mean_error_suffix = sprintf(" (%.3e %s in range)", mean_error_in_range, output_units);
            max_error_suffix = sprintf(" (%.3e %s in range)", max_error_in_range, output_units);
        else
            L2_error_suffix = " (no samples in range)";
            mean_error_suffix = " (no samples in range)";
            max_error_suffix = " (no samples in range)";
        end
    end

    fprintf('[%s] Mean Square error : %.3e %s%s\n',func_name, L2_error, output_units, ...
        L2_error_suffix);
    fprintf('[%s] Mean difference   : %.3e %s%s\n',func_name, mean_error, output_units, ...
        mean_error_suffix);
    fprintf('[%s] Max  difference   : %.3e %s%s\n',func_name, max_error, output_units, ...
        max_error_suffix);
    
    if ishandle(fig_num)
        close(fig_num);
    end

    fig = figure(fig_num);
    title(sprintf("RAIPAL Crossed 4-bar Polynomial Analysis: %s",func_name))
    xlabel(sprintf("%s [%s]",input_name,input_units))
    xlim([min(input_plot) max(input_plot)]);
    ylim([min(output_fit_plot) max(output_fit_plot)]);
    yyaxis left
    plot(input_plot,output_table_plot,"DisplayName","Numerical Solution")
    hold on;
    plot(input_plot,output_fit_plot,"--","DisplayName","Polynomial Fit")
    if has_range(range_map, output_name)
        output_range = plot_range(range_map, output_name, output_units);
        yline(output_range(1),"HandleVisibility","off");
        yline(output_range(2),"HandleVisibility","off");
    end
    ylabel(sprintf("%s [%s]",output_name,output_units))

    if has_range(range_map, input_name)
        input_range = plot_range(range_map, input_name, input_units);
        xline(input_range(1),"HandleVisibility","off");
        xline(input_range(2),"HandleVisibility","off");
    end

    yyaxis right
    plot(input_plot,abs(fit_err),"DisplayName",sprintf("Approx. Err."))
    ylabel(sprintf("Approx. Err. [%s]",output_units))

    grid on;
    legend("Location","northwest");
    hold off;

    safe_func_name = erase(strrep(func_name,"/","-"),"\");
    exportgraphics(fig,sprintf("results/FIGURE-%d_%s.png",fig_num, safe_func_name),"Resolution",300)

end

function [p, max_abs_error, fit_method] = normalized_polyfit(x_norm, y, degree, fit_mode)
    fit_mode = string(fit_mode);

    if fit_mode == "least-squares"
        p = polyfit(x_norm, y, degree);
        max_abs_error = max(abs(y - polyval(p, x_norm)));
        fit_method = "least-squares";
        return;
    end

    if fit_mode ~= "minimax" && fit_mode ~= "auto"
        error('Unknown fit_mode "%s". Use "least-squares", "minimax", or "auto".', fit_mode);
    end

    if exist('linprog','file') ~= 2
        if fit_mode == "minimax"
            warning('linprog is unavailable. Falling back to least-squares polyfit.');
        end
        p = polyfit(x_norm, y, degree);
        max_abs_error = max(abs(y - polyval(p, x_norm)));
        fit_method = "least-squares fallback";
        return;
    end

    A = zeros(length(x_norm), degree + 1);
    for k = 0:degree
        A(:,degree + 1 - k) = x_norm.^k;
    end

    % Decision variables are [p(:); t], minimizing t subject to
    % -t <= y - A*p <= t over all sampled table points.
    f = [zeros(degree + 1, 1); 1];
    Aineq = [ A, -ones(length(x_norm), 1);
             -A, -ones(length(x_norm), 1)];
    bineq = [ y;
             -y];
    lb = [-inf(degree + 1, 1); 0];

    [sol, exitflag, output] = solve_minimax_lp(f, Aineq, bineq, lb);

    if isempty(sol) || exitflag <= 0
        warning('linprog minimax fit failed: %s Falling back to least-squares polyfit.', output.message);
        p = polyfit(x_norm, y, degree);
        max_abs_error = max(abs(y - polyval(p, x_norm)));
        fit_method = "least-squares fallback";
        return;
    end

    p = sol(1:end-1).';
    max_abs_error = sol(end);
    fit_method = "minimax";
end

function [sol, exitflag, output] = solve_minimax_lp(f, Aineq, bineq, lb)
    algorithms = ["dual-simplex", "interior-point"];
    sol = [];
    exitflag = -1;
    output = struct("message", "No linprog algorithm was attempted.");

    for i = 1:numel(algorithms)
        try
            options = optimoptions('linprog', ...
                'Algorithm', algorithms(i), ...
                'Display', 'none');
            [candidate_sol, ~, candidate_exitflag, candidate_output] = ...
                linprog(f, Aineq, bineq, [], [], lb, [], options);
        catch err
            candidate_sol = [];
            candidate_exitflag = -1;
            candidate_output = struct("message", err.message);
        end

        if ~isempty(candidate_sol) && candidate_exitflag > 0
            sol = candidate_sol;
            exitflag = candidate_exitflag;
            output = candidate_output;
            return;
        end

        sol = candidate_sol;
        exitflag = candidate_exitflag;
        output = candidate_output;
    end
end

function x_norm = normalize_input(input, x_min, x_max)
    x_norm = 2 * (input - x_min) / (x_max - x_min) - 1;
end

function factor = unit_factor(units)
    units = string(units);
    factor = 1;

    if units == "deg"
        factor = 180/pi;
    elseif units == "deg/deg"
        factor = 1;
    elseif units == "deg^{-1}"
        factor = pi/180;
    elseif units == "rad" || units == "rad/rad" || units == "rad^{-1}"
        factor = 1;
    end
end

function tf = has_range(range_map, range_name)
    tf = any(string(keys(range_map)) == string(range_name));
end

function values = plot_range(range_map, range_name, units)
    values = range_map(string(range_name));
    values = values{1};
    values = values * unit_factor(units);
end
