%% Paper CFB figures
% Creates:
%   figures/PAPER_reduction_with-markers.png
%   figures/PAPER_operating-regions_with-markers.png

clear; clc; close all;

script_dir = fileparts(mfilename("fullpath"));
out_dir = fullfile(script_dir, "figures");
if ~exist(out_dir, "dir")
    mkdir(out_dir);
end

%% Global font scaling (for paper figures)
set(groot, "defaultAxesFontSize",   12);
set(groot, "defaultTextFontSize",   15);
set(groot, "defaultLegendFontSize", 12);
set(groot, "defaultColorbarFontSize", 12);
set(groot, 'DefaultLineLineWidth', 2);

blue   = [0.0000 0.4470 0.7410];
orange = [0.8500 0.3250 0.0980];
yellow = [0.9290 0.6940 0.1250];
purple = [0.4940 0.1840 0.5560];
colors = [blue; orange; yellow; purple];

%% Effective reduction ratio data
% Polynomial coefficients are from the fitted CFB solution used in
% VISUALIZATION/proposal_log_analysis.m. They are ordered highest degree to
% lowest degree for polyval().
coeff_b_of_a = [
     7.487348861572753, -88.17268027058002, 470.7676589474978, ...
    -1507.949571146181, 3232.924175366478, -4903.587201465159, ...
     5428.248858581333, -4464.681332539945, 2751.745273804874, ...
    -1265.589676288956, 418.4687974259759, -86.36058832217348, ...
     6.540591054554072, -1.820824957906633, 2.514881046035373, ...
     1.022698355939166, 0.000002629580961635315];

coeff_wb_of_a = [
     330.1514673459376, -3939.920472235199, 21319.39656007634, ...
    -69182.34049371608, 150034.7654546927, -229354.9209106861, ...
     253935.6748254774, -205870.6560746802, 121976.630364867, ...
    -52021.60108263157, 15489.35547273981, -3096.107249114025, ...
     433.0072547668175, -57.21741621150876, -0.9670218724995983, ...
     4.921276349693374, 1.022525952484172];

a_domain = linspace(0, 89.13041976, 4000)' * pi / 180;
elbow_deg = polyval(coeff_b_of_a, a_domain) * 180 / pi;
ratio = 1 ./ polyval(coeff_wb_of_a, a_domain);

extrema_idx = elbow_deg >= 15 & elbow_deg <= 120;
elbow_extrema = elbow_deg(extrema_idx);
ratio_extrema = ratio(extrema_idx);
[~, local_min_idx] = min(ratio_extrema);
[~, local_max_idx] = max(ratio_extrema);

marker_elbow_deg = [
    elbow_extrema(local_min_idx), ...
    50, ...
    75, ...
    elbow_extrema(local_max_idx)];
marker_ratio = interp1(elbow_deg, ratio, marker_elbow_deg, "pchip");

%% Figure 1: effective reduction ratio
fig = figure("Color", "w", "Units", "centimeters", "Position", [2 2 12.0 9.0]);
ax = axes(fig);
hold(ax, "on");
grid(ax, "on");
box(ax, "on");

plot(ax, elbow_deg, ratio, "Color", blue);
for i = 1:numel(marker_elbow_deg)
    plot(ax, marker_elbow_deg(i), marker_ratio(i), "o", ...
        "MarkerSize", 4.2, ...
        "MarkerFaceColor", "w", ...
        "MarkerEdgeColor", colors(i, :), ...
        "LineWidth", 1.4);

    text(ax, marker_elbow_deg(i), marker_ratio(i) - 0.012, sprintf("%d", i), ...
        "HorizontalAlignment", "center", ...
        "VerticalAlignment", "top", ...
        "FontWeight", "bold", ...
        "Color", colors(i, :));
end

title(ax, "Elbow Effective Reduction Ratio", "FontWeight", "bold");
xlabel(ax, "Elbow Angle [deg]");
ylabel(ax, "Reduction Ratio");
xlim(ax, [15 120]);
ylim(ax, [0.5 0.9]);
xticks(ax, 20:20:120);
yticks(ax, 0.5:0.1:0.9);

exportgraphics(fig, fullfile(out_dir, "PAPER_reduction_with-markers.png"), ...
    "Resolution", 300);

%% Figure 2: elbow joint operating region
load(fullfile(script_dir, "MOR", "MOR_MHR.mat"), "MOR_keypoints_MHR_peak");

fig = figure("Color", "w", "Units", "centimeters", "Position", [2 2 12.0 9.0]);
ax = axes(fig);
hold(ax, "on");
grid(ax, "on");
box(ax, "on");

xlim(ax, [-30 30]);
ylim(ax, [-220 220]);
xticks(ax, -30:10:30);
yticks(ax, -200:100:200);

xline(ax, 0, "-", "Color", [0.45 0.45 0.45], "LineWidth", 0.9, ...
    "HandleVisibility", "off");
yline(ax, 0, "-", "Color", [0.45 0.45 0.45], "LineWidth", 0.9, ...
    "HandleVisibility", "off");

for i = 1:numel(marker_ratio)
    r = marker_ratio(i);
    x_or = MOR_keypoints_MHR_peak(:, 1) ./ r;
    y_or = MOR_keypoints_MHR_peak(:, 2) .* r;

    plot(ax, x_or, y_or, "-", ...
        "Color", colors(i, :), ...
        "LineWidth", 2.0, ...
        "HandleVisibility", "off");

    marker_idx = 4; % positive-speed, negative-torque vertex
    x_pt = x_or(marker_idx);
    y_pt = y_or(marker_idx);
    plot(ax, x_pt, y_pt, "o", ...
        "MarkerSize", 4.6, ...
        "MarkerFaceColor", "w", ...
        "MarkerEdgeColor", colors(i, :), ...
        "LineWidth", 1.4);

    text(ax, x_pt + 1.0, y_pt - 3.0, sprintf("%d", i), ...
        "Color", colors(i, :), ...
        "FontWeight", "bold", ...
        "HorizontalAlignment", "left", ...
        "VerticalAlignment", "top");
end

plot(ax, MOR_keypoints_MHR_peak(:, 1), MOR_keypoints_MHR_peak(:, 2), "--", ...
    "Color", [0 0 1], ...
    "LineWidth", 2.0, ...
    "HandleVisibility", "off");

text(ax, 8.0, 190, "Actuator", ...
    "Color", [0 0 1], ...
    "FontWeight", "bold", ...
    "HorizontalAlignment", "left", ...
    "VerticalAlignment", "top");

title(ax, "Elbow Joint Operating Region", "FontWeight", "bold");
xlabel(ax, "Speed [rad/s]");
ylabel(ax, "Torque [Nm]");

exportgraphics(fig, fullfile(out_dir, "PAPER_operating-regions_with-markers.png"), ...
    "Resolution", 300);

fprintf("Saved %s\n", fullfile(out_dir, "PAPER_reduction_with-markers.png"));
fprintf("Saved %s\n", fullfile(out_dir, "PAPER_operating-regions_with-markers.png"));
