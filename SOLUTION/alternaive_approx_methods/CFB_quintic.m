%%
clc;
clear;
close all;

%%
load("results/cfb_table.mat")

%% Quintic Hermite settings
num_segments = 28; % 5 deg segments over the full -5 to 135 deg table range

range = dictionary;
range("joint position") = {[0 130] / 180 * pi};
range("actuator position") = {[0 89.13041976] / 180 * pi};

%% Build quintic Hermite approximation for a(b)
coeffs = struct();
[coeffs.jointpos_to_actuatorpos, ~] = quintic_hermite_analyze( ...
    b_table, "joint position", ...
    a_table, da_table, dda_table, "actuator position", ...
    num_segments, 'deg', 1, range);

save('results/cfb_quintic_coeffs.mat',"coeffs")

%% function definitions
function [fit_data,fit] = quintic_hermite_analyze( ...
    input_table,input_name, ...
    output_table,doutput_table,ddoutput_table,output_name, ...
    num_segments,fig_units,fig_num,range_map)

    func_name = sprintf('%s -> %s',input_name,output_name);
    fit_range = [min(input_table) max(input_table)];

    edges = linspace(fit_range(1), fit_range(2), num_segments + 1);
    segment_widths = diff(edges);

    y_edge   = interp1(input_table, output_table,   edges, 'pchip');
    dy_edge  = interp1(input_table, doutput_table,  edges, 'pchip');
    ddy_edge = interp1(input_table, ddoutput_table, edges, 'pchip');

    % Coefficients are ascending powers of local s:
    %   y = c0 + c1*s + c2*s^2 + c3*s^3 + c4*s^4 + c5*s^5
    % where s = (input - edge_left) / segment_width.
    segment_coeffs = zeros(num_segments, 6);

    for i = 1:num_segments
        h = segment_widths(i);

        y0 = y_edge(i);
        y1 = y_edge(i + 1);
        v0 = h * dy_edge(i);
        v1 = h * dy_edge(i + 1);
        a0 = h^2 * ddy_edge(i);
        a1 = h^2 * ddy_edge(i + 1);

        c0 = y0;
        c1 = v0;
        c2 = 0.5 * a0;

        d0 = y1 - (c0 + c1 + c2);
        d1 = v1 - (c1 + 2*c2);
        d2 = a1 - 2*c2;

        c3 = 10*d0 - 4*d1 + 0.5*d2;
        c4 = -15*d0 + 7*d1 - d2;
        c5 = 6*d0 - 3*d1 + 0.5*d2;

        segment_coeffs(i,:) = [c0 c1 c2 c3 c4 c5];
    end

    [fit.output, fit.doutput, fit.ddoutput] = eval_quintic_segments( ...
        input_table, edges, segment_widths, segment_coeffs);

    fit_data = struct();
    fit_data.degree = 5;
    fit_data.num_segments = num_segments;
    fit_data.edges = edges;
    fit_data.segment_widths = segment_widths;
    fit_data.coeffs = segment_coeffs;
    fit_data.coeff_order = "ascending powers of normalized local s";
    fit_data.input_name = input_name;
    fit_data.output_name = output_name;
    fit_data.input_range = fit_range;

    fprintf('\n\n[%s] Quintic Hermite coefficients:\n', func_name);
    fprintf('Rows are segments; columns are c0..c5 for normalized local s in [0,1].\n');
    disp(segment_coeffs);

    print_error_report(func_name, "position", input_table, output_table, fit.output, fig_units, range_map, input_name);
    print_error_report(func_name, "velocity", input_table, doutput_table, fit.doutput, "ratio", range_map, input_name);
    print_error_report(func_name, "acceleration", input_table, ddoutput_table, fit.ddoutput, "1/rad", range_map, input_name);

    coeff_abs = abs(segment_coeffs(segment_coeffs ~= 0));
    coeff_ratio = max(coeff_abs) / min(coeff_abs);
    fprintf('[%s] Coefficient magnitude ratio: %.3e\n',func_name, coeff_ratio);

    plot_quintic_report( ...
        input_table,input_name, ...
        output_table,doutput_table,ddoutput_table, ...
        fit,output_name, ...
        edges,num_segments,fig_units,fig_num,range_map,func_name);
end

function [output,doutput,ddoutput] = eval_quintic_segments(input,edges,segment_widths,segment_coeffs)
    num_segments = size(segment_coeffs, 1);
    output = nan(size(input));
    doutput = nan(size(input));
    ddoutput = nan(size(input));

    for i = 1:num_segments
        if i < num_segments
            segment_idx = input >= edges(i) & input < edges(i + 1);
        else
            segment_idx = input >= edges(i) & input <= edges(i + 1);
        end

        h = segment_widths(i);
        s = (input(segment_idx) - edges(i)) ./ h;
        c = segment_coeffs(i,:);

        output(segment_idx) = c(1) + s .* (c(2) + s .* (c(3) + s .* (c(4) + s .* (c(5) + s .* c(6)))));
        doutput(segment_idx) = (c(2) + s .* (2*c(3) + s .* (3*c(4) + s .* (4*c(5) + s .* 5*c(6))))) ./ h;
        ddoutput(segment_idx) = (2*c(3) + s .* (6*c(4) + s .* (12*c(5) + s .* 20*c(6)))) ./ h^2;
    end
end

function print_error_report(func_name,quantity_name,input_table,true_table,fit_table,units,range_map,input_name)
    valid_idx = ~isnan(fit_table);
    err = true_table - fit_table;

    if contains(units,'deg')
        err = err * 180/pi;
    end

    L2_error = norm(err(valid_idx))^2 / nnz(valid_idx);
    mean_error = mean(abs(err(valid_idx)));
    max_error = max(abs(err(valid_idx)));

    L2_error_suffix = "";
    mean_error_suffix = "";
    max_error_suffix = "";

    if has_range(range_map, input_name)
        input_range = raw_range(range_map, input_name);
        in_range_idx = input_table >= input_range(1) & input_table <= input_range(2) & valid_idx;

        if any(in_range_idx)
            err_in_range = err(in_range_idx);
            L2_error_in_range = norm(err_in_range)^2 / nnz(in_range_idx);
            mean_error_in_range = mean(abs(err_in_range));
            max_error_in_range = max(abs(err_in_range));

            L2_error_suffix = sprintf(" (%.3e %s in range)", L2_error_in_range, units);
            mean_error_suffix = sprintf(" (%.3e %s in range)", mean_error_in_range, units);
            max_error_suffix = sprintf(" (%.3e %s in range)", max_error_in_range, units);
        else
            L2_error_suffix = " (no samples in range)";
            mean_error_suffix = " (no samples in range)";
            max_error_suffix = " (no samples in range)";
        end
    end

    fprintf('[%s] %s Mean Square error : %.3e %s%s\n',func_name,quantity_name,L2_error,units,L2_error_suffix);
    fprintf('[%s] %s Mean difference   : %.3e %s%s\n',func_name,quantity_name,mean_error,units,mean_error_suffix);
    fprintf('[%s] %s Max  difference   : %.3e %s%s\n',func_name,quantity_name,max_error,units,max_error_suffix);
end

function plot_quintic_report( ...
    input_table,input_name, ...
    output_table,doutput_table,ddoutput_table, ...
    fit,output_name, ...
    edges,num_segments,fig_units,fig_num,range_map,func_name)

    input_plot = input_table;
    output_plot = output_table;
    output_fit_plot = fit.output;

    if contains(fig_units,'deg')
        input_plot = input_plot * 180/pi;
        output_plot = output_plot * 180/pi;
        output_fit_plot = output_fit_plot * 180/pi;
    end

    edge_plot = edges;
    if contains(fig_units,'deg')
        edge_plot = edge_plot * 180/pi;
    end

    if ishandle(fig_num)
        close(fig_num);
    end

    fig = figure(fig_num);
    fig.Position = [100, 100, 1200, 900];

    subplot(3,1,1)
    yyaxis left
    plot(input_plot,output_plot,"DisplayName","Numerical Solution")
    hold on;
    plot(input_plot,output_fit_plot,"--","DisplayName","Quintic Hermite Fit")
    ylabel(sprintf("%s [%s]",output_name,fig_units))
    add_range_lines(input_name, output_name, range_map, fig_units, true);

    yyaxis right
    plot(input_plot,abs(output_plot - output_fit_plot),"DisplayName","Approx. Err.")
    ylabel(sprintf("Approx. Err. [%s]",fig_units))
    title(sprintf("RAIPAL Crossed 4-bar Quintic Hermite: %s",func_name))
    grid on;
    legend("Location","northwest");

    subplot(3,1,2)
    yyaxis left
    plot(input_plot,doutput_table,"DisplayName","Numerical da/db")
    hold on;
    plot(input_plot,fit.doutput,"--","DisplayName","Hermite da/db")
    ylabel("da/db")
    add_range_lines(input_name, output_name, range_map, fig_units, false);

    yyaxis right
    plot(input_plot,abs(doutput_table - fit.doutput),"DisplayName","Derivative Err.")
    ylabel("Derivative Err.")
    grid on;
    legend("Location","northwest");

    subplot(3,1,3)
    yyaxis left
    plot(input_plot,ddoutput_table,"DisplayName","Numerical d^2a/db^2")
    hold on;
    plot(input_plot,fit.ddoutput,"--","DisplayName","Hermite d^2a/db^2")
    ylabel("d^2a/db^2")
    xlabel(sprintf("%s [%s]",input_name,fig_units(1:3)))
    add_range_lines(input_name, output_name, range_map, fig_units, false);

    yyaxis right
    plot(input_plot,abs(ddoutput_table - fit.ddoutput),"DisplayName","Second Derivative Err.")
    ylabel("Second Derivative Err.")
    grid on;
    legend("Location","northwest");

    for ax = findall(fig,'Type','axes')'
        axes(ax);
        xlim([min(input_plot) max(input_plot)]);
        for i = 2:num_segments
            xline(edge_plot(i),":","HandleVisibility","off");
        end
    end

    safe_func_name = erase(strrep(func_name,"/","-"),"\");
    exportgraphics(fig,sprintf("results/FIGURE-%d_quintic_%s.png",fig_num, safe_func_name),"Resolution",300)
end

function add_range_lines(input_name,output_name,range_map,fig_units,show_output_range)
    if has_range(range_map, input_name)
        input_range = plot_range(range_map, input_name, fig_units);
        xline(input_range(1),"HandleVisibility","off");
        xline(input_range(2),"HandleVisibility","off");
    end

    if show_output_range && has_range(range_map, output_name)
        output_range = plot_range(range_map, output_name, fig_units);
        yline(output_range(1),"HandleVisibility","off");
        yline(output_range(2),"HandleVisibility","off");
    end
end

function tf = has_range(range_map, range_name)
    tf = any(string(keys(range_map)) == string(range_name));
end

function values = raw_range(range_map, range_name)
    values = range_map(string(range_name));
    values = values{1};
end

function values = plot_range(range_map, range_name, fig_units)
    values = raw_range(range_map, range_name);
    if contains(fig_units,'deg')
        values = values * 180/pi;
    end
end
