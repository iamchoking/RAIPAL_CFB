%%
clc;
clear;
close all;

%%
load("../results/cfb_table.mat")

%% Comparison settings
simple_poly_degree = 15;
piecewise_degree = 5;
num_segments = 28; % 5 deg segments over the full -5 to 135 deg table range
num_timing_repeats = 200;

input = b_table;
truth.a = a_table;
truth.da = da_table;
truth.dda = dda_table;

fit_range = [min(input) max(input)];
edges = linspace(fit_range(1), fit_range(2), num_segments + 1);

%% Build approximations
models = {};

models{end + 1} = build_linear_model(input, truth.a, edges);
models{end + 1} = build_cubic_model(input, truth.a, edges);
models{end + 1} = build_quintic_model(input, truth.a, truth.da, truth.dda, edges);
models{end + 1} = build_simple_poly_model(input, truth, simple_poly_degree);
models{end + 1} = build_piecewise_poly_model(input, truth, edges, piecewise_degree);

%% Evaluate and time
names = strings(numel(models), 1);
coeff_counts = zeros(numel(models), 1);
rms_a_deg = zeros(numel(models), 1);
max_a_deg = zeros(numel(models), 1);
rms_da = zeros(numel(models), 1);
max_da = zeros(numel(models), 1);
rms_dda = zeros(numel(models), 1);
max_dda = zeros(numel(models), 1);
combined_rms = zeros(numel(models), 1);
time_per_call_us = zeros(numel(models), 1);
time_per_sample_ns = zeros(numel(models), 1);

for i = 1:numel(models)
    model = models{i};
    fit = model.eval(input);

    err_a_deg = (truth.a - fit.a) * 180/pi;
    err_da = truth.da - fit.da;
    err_dda = truth.dda - fit.dda;

    names(i) = model.name;
    coeff_counts(i) = model.num_coeffs;

    rms_a_deg(i) = sqrt(mean(err_a_deg.^2));
    max_a_deg(i) = max(abs(err_a_deg));
    rms_da(i) = sqrt(mean(err_da.^2));
    max_da(i) = max(abs(err_da));
    rms_dda(i) = sqrt(mean(err_dda.^2));
    max_dda(i) = max(abs(err_dda));
    combined_rms(i) = sqrt(mean([err_a_deg(:); err_da(:); err_dda(:)].^2));

    elapsed = time_inference(model.eval, input, num_timing_repeats);
    time_per_call_us(i) = elapsed / num_timing_repeats * 1e6;
    time_per_sample_ns(i) = elapsed / num_timing_repeats / numel(input) * 1e9;
end

results = table(names, coeff_counts, ...
    rms_a_deg, max_a_deg, ...
    rms_da, max_da, ...
    rms_dda, max_dda, ...
    combined_rms, ...
    time_per_call_us, time_per_sample_ns, ...
    'VariableNames', {'method','num_coeffs', ...
    'rms_a_deg','max_a_deg', ...
    'rms_da','max_da', ...
    'rms_dda','max_dda', ...
    'combined_rms', ...
    'time_per_call_us','time_per_sample_ns'});

disp(results)

if ~exist("results", "dir")
    mkdir("results");
end

save("results/cfb_approx_comparison.mat", "results")

%% Plot comparison
fig = figure(1);
fig.Position = [100, 100, 1400, 900];

subplot(2,3,1)
bar(categorical(names), coeff_counts)
ylabel("Stored scalar coefficients")
title("Total Coefficient Count")
grid on

subplot(2,3,2)
bar(categorical(names), rms_a_deg)
ylabel("RMS a error [deg]")
title("Position RMS Error")
grid on

subplot(2,3,3)
bar(categorical(names), max_a_deg)
ylabel("Max a error [deg]")
title("Position Max Error")
grid on

subplot(2,3,4)
bar(categorical(names), rms_da)
ylabel("RMS da/db error")
title("Velocity RMS Error")
grid on

subplot(2,3,5)
bar(categorical(names), rms_dda)
ylabel("RMS d^2a/db^2 error")
title("Acceleration RMS Error")
grid on

subplot(2,3,6)
bar(categorical(names), time_per_sample_ns)
ylabel("Time per sample [ns]")
title("MATLAB Inference Time")
grid on

exportgraphics(fig, "results/FIGURE_approx_comparison.png", "Resolution", 300)

%% Local functions
function model = build_linear_model(input, output, edges)
    y_edge = interp1(input, output, edges, 'pchip');
    widths = diff(edges);
    coeffs = zeros(numel(widths), 2);

    for i = 1:numel(widths)
        coeffs(i,:) = [y_edge(i), y_edge(i + 1) - y_edge(i)];
    end

    model = struct();
    model.name = "linear interpolation";
    model.num_coeffs = numel(coeffs);
    model.eval = @(x) eval_piecewise_ascending_with_derivatives(x, edges, coeffs);
end

function model = build_cubic_model(input, output, edges)
    y_edge = interp1(input, output, edges, 'pchip');
    pp = pchip(edges, y_edge);

    model = struct();
    model.name = "cubic interpolation";
    model.num_coeffs = numel(pp.coefs);
    model.eval = @(x) eval_pp_with_derivatives(pp, x);
end

function model = build_quintic_model(input, output, doutput, ddoutput, edges)
    widths = diff(edges);

    y_edge = interp1(input, output, edges, 'pchip');
    dy_edge = interp1(input, doutput, edges, 'pchip');
    ddy_edge = interp1(input, ddoutput, edges, 'pchip');

    coeffs = zeros(numel(widths), 6);

    for i = 1:numel(widths)
        h = widths(i);

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

        coeffs(i,:) = [c0 c1 c2 c3 c4 c5];
    end

    model = struct();
    model.name = "quintic Hermite";
    model.num_coeffs = numel(coeffs);
    model.eval = @(x) eval_piecewise_ascending_with_derivatives(x, edges, coeffs);
end

function model = build_simple_poly_model(input, truth, degree)
    coeffs_a = polyfit(input, truth.a, degree);
    coeffs_da = polyfit(input, truth.da, degree);
    coeffs_dda = polyfit(input, truth.dda, degree);

    model = struct();
    model.name = sprintf("simple polynomial d%d", degree);
    model.num_coeffs = numel(coeffs_a) + numel(coeffs_da) + numel(coeffs_dda);
    model.eval = @(x) eval_three_poly(coeffs_a, coeffs_da, coeffs_dda, x);
end

function model = build_piecewise_poly_model(input, truth, edges, degree)
    centers = 0.5 * (edges(1:end-1) + edges(2:end));
    half_widths = 0.5 * diff(edges);

    coeffs_a = fit_piecewise_centered_poly(input, truth.a, edges, centers, half_widths, degree);
    coeffs_da = fit_piecewise_centered_poly(input, truth.da, edges, centers, half_widths, degree);
    coeffs_dda = fit_piecewise_centered_poly(input, truth.dda, edges, centers, half_widths, degree);

    model = struct();
    model.name = sprintf("piecewise polynomial d%d", degree);
    model.num_coeffs = numel(coeffs_a) + numel(coeffs_da) + numel(coeffs_dda);
    model.eval = @(x) eval_three_piecewise_centered_poly( ...
        x, edges, centers, half_widths, coeffs_a, coeffs_da, coeffs_dda);
end

function coeffs = fit_piecewise_centered_poly(input, output, edges, centers, half_widths, degree)
    coeffs = zeros(numel(half_widths), degree + 1);

    for i = 1:numel(half_widths)
        if i < numel(half_widths)
            idx = input >= edges(i) & input < edges(i + 1);
        else
            idx = input >= edges(i) & input <= edges(i + 1);
        end

        x_local = (input(idx) - centers(i)) ./ half_widths(i);
        coeffs(i,:) = polyfit(x_local, output(idx), degree);
    end
end

function fit = eval_piecewise_ascending_with_derivatives(x, edges, coeffs)
    fit.a = nan(size(x));
    fit.da = nan(size(x));
    fit.dda = nan(size(x));

    num_segments = size(coeffs, 1);
    widths = diff(edges);

    for i = 1:num_segments
        if i < num_segments
            idx = x >= edges(i) & x < edges(i + 1);
        else
            idx = x >= edges(i) & x <= edges(i + 1);
        end

        h = widths(i);
        s = (x(idx) - edges(i)) ./ h;
        c = coeffs(i,:);

        fit.a(idx) = polyval(fliplr(c), s);

        if numel(c) >= 2
            dc = (1:numel(c)-1) .* c(2:end);
            fit.da(idx) = polyval(fliplr(dc), s) ./ h;
        else
            fit.da(idx) = 0;
        end

        if numel(c) >= 3
            ddc = (2:numel(c)-1) .* (1:numel(c)-2) .* c(3:end);
            fit.dda(idx) = polyval(fliplr(ddc), s) ./ h^2;
        else
            fit.dda(idx) = 0;
        end
    end
end

function fit = eval_pp_with_derivatives(pp, x)
    fit.a = ppval(pp, x);
    fit.da = ppval(fnder_local(pp, 1), x);
    fit.dda = ppval(fnder_local(pp, 2), x);
end

function ppd = fnder_local(pp, derivative_order)
    ppd = pp;
    for k = 1:derivative_order
        order = size(ppd.coefs, 2);
        powers = order-1:-1:1;
        ppd.coefs = ppd.coefs(:,1:end-1) .* powers;
        ppd.order = ppd.order - 1;
    end
end

function fit = eval_three_poly(coeffs_a, coeffs_da, coeffs_dda, x)
    fit.a = polyval(coeffs_a, x);
    fit.da = polyval(coeffs_da, x);
    fit.dda = polyval(coeffs_dda, x);
end

function fit = eval_three_piecewise_centered_poly(x, edges, centers, half_widths, coeffs_a, coeffs_da, coeffs_dda)
    fit.a = eval_piecewise_centered_poly(x, edges, centers, half_widths, coeffs_a);
    fit.da = eval_piecewise_centered_poly(x, edges, centers, half_widths, coeffs_da);
    fit.dda = eval_piecewise_centered_poly(x, edges, centers, half_widths, coeffs_dda);
end

function y = eval_piecewise_centered_poly(x, edges, centers, half_widths, coeffs)
    y = nan(size(x));
    num_segments = size(coeffs, 1);

    for i = 1:num_segments
        if i < num_segments
            idx = x >= edges(i) & x < edges(i + 1);
        else
            idx = x >= edges(i) & x <= edges(i + 1);
        end

        x_local = (x(idx) - centers(i)) ./ half_widths(i);
        y(idx) = polyval(coeffs(i,:), x_local);
    end
end

function elapsed = time_inference(eval_fun, input, num_repeats)
    eval_fun(input); % warm up

    if exist('timeit', 'file') == 2
        elapsed = timeit(@() repeated_eval(eval_fun, input, num_repeats));
    else
        tic;
        repeated_eval(eval_fun, input, num_repeats);
        elapsed = toc;
    end
end

function repeated_eval(eval_fun, input, num_repeats)
    for i = 1:num_repeats
        eval_fun(input);
    end
end
