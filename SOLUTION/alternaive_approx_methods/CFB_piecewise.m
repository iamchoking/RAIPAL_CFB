%%
clc;
clear;
close all;

%%
load("results/cfb_table.mat")

%% Piecewise polynomial settings
degree = 5;
num_segments = 28; % 10 deg segments over the full -5 to 135 deg table range

range = dictionary;
range("joint position") = {[0 130] / 180 * pi};
range("actuator position") = {[0 89.13041976] / 180 * pi};

%% Fit piecewise polynomial for a
coeffs = struct();
[coeffs.jointpos_to_actuatorpos, ~] = piecewise_polyfit_analyze( ...
    b_table, "joint position", ...
    a_table, "actuator position", ...
    degree, num_segments, 'deg', 1, range);

save('results/cfb_piecewise_coeffs.mat',"coeffs")

%% function definitions
function [fit_data,output_fit] = piecewise_polyfit_analyze( ...
    input_table,input_name,output_table,output_name,degree,num_segments,fig_units,fig_num,range_map)

    func_name = sprintf('%s -> %s',input_name,output_name);
    num_points = length(input_table);

    fit_range = [min(input_table) max(input_table)];

    edges = linspace(fit_range(1), fit_range(2), num_segments + 1);
    centers = 0.5 * (edges(1:end-1) + edges(2:end));
    half_widths = 0.5 * diff(edges);
    segment_coeffs = zeros(num_segments, degree + 1);
    output_fit = nan(size(output_table));

    for i = 1:num_segments
        if i < num_segments
            segment_idx = input_table >= edges(i) & input_table < edges(i + 1);
        else
            segment_idx = input_table >= edges(i) & input_table <= edges(i + 1);
        end

        if nnz(segment_idx) < degree + 1
            error('Segment %d has only %d points. Need at least %d.', ...
                i, nnz(segment_idx), degree + 1);
        end

        x_segment = (input_table(segment_idx) - centers(i)) ./ half_widths(i);
        segment_coeffs(i,:) = polyfit(x_segment, output_table(segment_idx), degree);
        output_fit(segment_idx) = polyval(segment_coeffs(i,:), x_segment);
    end

    fit_data = struct();
    fit_data.degree = degree;
    fit_data.num_segments = num_segments;
    fit_data.edges = edges;
    fit_data.centers = centers;
    fit_data.half_widths = half_widths;
    fit_data.coeffs = segment_coeffs;
    fit_data.input_name = input_name;
    fit_data.output_name = output_name;
    fit_data.input_range = fit_range;

    fprintf('\n\n[%s] Piecewise degree-%d polynomial coefficients:\n', func_name, degree);
    fprintf('Rows are segments; columns are highest to lowest powers of normalized local x.\n');
    disp(vpa(segment_coeffs,16));

    input_plot = input_table;
    output_plot = output_table;
    output_fit_plot = output_fit;

    if contains(fig_units,'deg')
        input_plot = input_plot * 180/pi;
        output_plot = output_plot * 180/pi;
        output_fit_plot = output_fit_plot * 180/pi;
    end

    fit_err = output_plot - output_fit_plot;
    valid_idx = ~isnan(output_fit);

    L2_error = norm(fit_err(valid_idx))^2 / nnz(valid_idx);
    mean_error = mean(abs(fit_err(valid_idx)));
    max_error = max(abs(fit_err(valid_idx)));

    L2_error_suffix = "";
    mean_error_suffix = "";
    max_error_suffix = "";

    if has_range(range_map, input_name)
        input_range = raw_range(range_map, input_name);
        in_range_idx = input_table >= input_range(1) & input_table <= input_range(2) & valid_idx;

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

    coeff_abs = abs(segment_coeffs(segment_coeffs ~= 0));
    coeff_ratio = max(coeff_abs) / min(coeff_abs);
    fprintf('[%s] Coefficient magnitude ratio: %.3e\n',func_name, coeff_ratio);

    if ishandle(fig_num)
        close(fig_num);
    end

    fig = figure(fig_num);
    title(sprintf("RAIPAL Crossed 4-bar Piecewise Polynomial Analysis: %s",func_name))
    xlabel(sprintf("%s [%s]",input_name,fig_units(1:3)))
    xlim([min(input_plot) max(input_plot)]);

    yyaxis left
    plot(input_plot,output_plot,"DisplayName","Numerical Solution")
    hold on;
    plot(input_plot,output_fit_plot,"--","DisplayName","Piecewise Polynomial Fit")
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

    edge_plot = edges;
    if contains(fig_units,'deg')
        edge_plot = edge_plot * 180/pi;
    end
    for i = 2:num_segments
        xline(edge_plot(i),":","HandleVisibility","off");
    end

    yyaxis right
    plot(input_plot,abs(fit_err),"DisplayName",sprintf("Approx. Err."))
    ylabel(sprintf("Approx. Err. [%s]",fig_units))

    grid on;
    legend("Location","northwest");
    hold off;

    safe_func_name = erase(strrep(func_name,"/","-"),"\");
    exportgraphics(fig,sprintf("results/FIGURE-%d_piecewise_%s.png",fig_num, safe_func_name),"Resolution",300)
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
