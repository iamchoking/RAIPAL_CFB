%%
clc;
clear;
close all;

%%
load("results/cfb_table.mat")
%% Polynomial degree
n = 15; %polynomial degree

range = dictionary;
range("joint position") = {[0 130] / 180 * pi};
range("actuator position") = {[0 89.13041976] / 180 * pi};

%% Fit polynomial for a
coeffs = struct();
[coeffs.jointpos_to_actuatorpos,~] = polyfit_analyze(b_table,"joint position",a_table,"actuator position",n,'deg',1,range);

save('results/cfb_coeffs.mat',"coeffs")

%% function definitions
function [coeffs,output_fit] = polyfit_analyze(input_table,input_name,output_table,output_name,num_degrees,fig_units,fig_num,range_map)

    func_name = sprintf('%s -> %s',input_name,output_name);
    num_points = length(input_table);

    coeffs = polyfit(input_table,output_table,num_degrees);
    % p(17) = 0;
    
    fprintf('\n\n[%s] Polynomial coefficients (highest to lowest degree):\n', func_name);
    disp(vpa(coeffs,16));
    
    output_fit = polyval(coeffs,input_table);
    input_plot = input_table;
    
    if(contains(fig_units,'deg'))
        input_plot = input_plot * 180/pi;
        output_fit = output_fit * 180/pi;
        output_table = output_table * 180/pi;
    end
    % Calculate Error
    fit_err    = output_table - output_fit;
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

            L2_error_suffix = sprintf(" (%.3e %s in range)", L2_error_in_range, fig_units);
            mean_error_suffix = sprintf(" (%.3e %s in range)", mean_error_in_range, fig_units);
            max_error_suffix = sprintf(" (%.3e %s in range)", max_error_in_range, fig_units);
        else
            L2_error_suffix = " (no samples in range)";
            mean_error_suffix = " (no samples in range)";
            max_error_suffix = " (no samples in range)";
        end
    end

    fprintf('[%s] Mean Square error : %.3e %s%s\n',func_name, L2_error, fig_units, ...
        L2_error_suffix);
    fprintf('[%s] Mean difference   : %.3e %s%s\n',func_name, mean_error, fig_units, ...
        mean_error_suffix);
    fprintf('[%s] Max  difference   : %.3e %s%s\n',func_name, max_error, fig_units, ...
        max_error_suffix);
    
    if ishandle(fig_num)
        close(fig_num);
    end

    fig = figure(fig_num);
    title(sprintf("RAIPAL Crossed 4-bar Polynomial Analysis: %s",func_name))
    xlabel(sprintf("%s [%s]",input_name,fig_units(1:3)))
    xlim([min(input_plot) max(input_plot)]);
    ylim([min(output_fit) max(output_fit)]);
    yyaxis left
    plot(input_plot,output_table,"DisplayName","Numerical Solution")
    hold on;
    plot(input_plot,output_fit,"--","DisplayName","Polynomial Fit")
    if has_range(range_map, output_name)
        output_range = plot_range(range_map, output_name, fig_units);
        yline(output_range(1),"HandleVisibility","off");
        yline(output_range(2),"HandleVisibility","off");
    end
    ylabel(sprintf("%s [%s]",output_name,fig_units))

    if has_range(range_map, input_name)
        input_range = plot_range(range_map, input_name, fig_units);
        xline(input_range(1),"HandleVisibility","off");
        xline(input_range(2),"HandleVisibility","off");
    end

    yyaxis right
    plot(input_plot,abs(fit_err),"DisplayName",sprintf("Approx. Err."))
    ylabel(sprintf("Approx. Err. [%s]",fig_units))

    grid on;
    legend("Location","northwest");
    hold off;

    safe_func_name = erase(strrep(func_name,"/","-"),"\");
    exportgraphics(fig,sprintf("results/FIGURE-%d_%s.png",fig_num, safe_func_name),"Resolution",300)

end

function tf = has_range(range_map, range_name)
    tf = any(string(keys(range_map)) == string(range_name));
end

function values = plot_range(range_map, range_name, fig_units)
    values = range_map(string(range_name));
    values = values{1};
    if contains(fig_units,'deg')
        values = values * 180/pi;
    end
end
