%% 
clear; clc;

%% Global font scaling (for paper figures)
set(groot, "defaultAxesFontSize",   14);
set(groot, "defaultTextFontSize",   12);
set(groot, "defaultLegendFontSize", 12);
set(groot, "defaultColorbarFontSize", 12);
set(groot, 'DefaultLineLineWidth', 2);

fig_size = [2 2 9 12];

%% Load data
num_points = 4000;
a_domain = linspace(0,89.13041976,num_points)' * pi / 180;
da = a_domain(2) - a_domain(1);

load('results/v1-5_result_2026-02-09.mat')
load('MOR/MOR_MHR.mat')

filename = "raipal_2026-02-16-05-00-10_pitching-demo_x1.25_20.6_elbow";
load("log_data/elbow_/" + filename)
load("MOR/MOR_MHR")

markers = [39.2, 39.45, 39.55, 39.61, 39.6415, 39.76];

% Impact-time offset
impact_time_offset = 39.6415;

% Only keep marker #3, #4, #5
keep_idx = [3 4 5];
markers_kept = markers(keep_idx);

%% Compute output torque from input torque and ratio (robust ratio index lookup)
log_gf_output = zeros(length(log_gf_input),1);
for i = 1:length(log_gf_input)
    [~, idx_ratio] = min(abs(a_domain - log_gc_input(i)));
    idx_ratio = min(max(idx_ratio,1), numel(inv_wb_fit));
    ratio_current = inv_wb_fit(idx_ratio);
    log_gf_output(i) = log_gf_input(i) * ratio_current;
end

%% Find nearest sample indices for prescribed marker times (with impact-time offset)
marker_idx = zeros(size(markers_kept));
for k = 1:numel(markers_kept)
    [~, marker_idx(k)] = min(abs(log_time - (markers_kept(k) - impact_time_offset)));
end

%% Single figure: output trace + selected highlighted points + instantaneous ORs
fig = figure('Color','w');
set(fig, 'Units', 'centimeters', 'Position', fig_size);

ax = axes(fig);
hold(ax, 'on'); grid(ax, 'on');
title(ax, "Torque-Speed (Elbow Joint)");
xlabel(ax, "Speed [rad/s]");
ylabel(ax, "Torque [Nm]");

% Axis limits
x_limits = [-30 30];
y_limits = [-190 190];
xlim(ax, x_limits);
ylim(ax, y_limits);

% Force square grid cells:
% same pixel size for 10 rad/s in x and 50 Nm in y
% daspect(ax, [10 50 1]);

% Major ticks aligned to desired grid increments
xticks(ax, x_limits(1):10:x_limits(2));
yticks(ax, -150:50:150);

colors = lines(numel(markers_kept));

% Full output torque-speed trace
plot(ax, log_gv_output, log_gf_output, 'Color', [0.2 0.2 0.2]);

% Label offsets
pt_dx = 0.0;    pt_dy = -8.0;    % directly below points
or_dx = 28.0;   or_dy = 0.0;     % near upper-left region area

for k = 1:numel(markers_kept)
    idx = marker_idx(k);

    % Ratio at this marker
    [~, idx_ratio] = min(abs(a_domain - log_gc_input(idx)));
    idx_ratio = min(max(idx_ratio,1), numel(inv_wb_fit));
    ratio_current = inv_wb_fit(idx_ratio);

    % Elbow joint angle (deg) for label
    elbow_deg = log_gc_input(idx) * 180/pi;

    % fprintf('kept-k=%d, marker#=%d, t=%.4f s, idx=%d, gc=%.6f rad (%.2f deg), ratio=%.6f\n', ...
    %     k, keep_idx(k), markers_kept(k), idx, log_gc_input(idx), elbow_deg, ratio_current);

    % Instantaneous operating region
    x_or = MOR_keypoints_MHR_peak(:,1) / ratio_current;
    y_or = MOR_keypoints_MHR_peak(:,2) * ratio_current;
    plot(ax, x_or, y_or, '-', 'Color', colors(k,:), 'LineWidth', 1.6);

    % Highlighted trace point
    x_pt = log_gv_output(idx);
    y_pt = log_gf_output(idx);
    scatter(ax, x_pt, y_pt, 55, ...
        'MarkerEdgeColor', colors(k,:), ...
        'MarkerFaceColor', 'none', ...
        'LineWidth', 1.9);

    % Point label: just 3 / 4 / 5, directly below point
    text(ax, x_pt + pt_dx, y_pt + pt_dy, sprintf('%d', keep_idx(k)), ...
        'Color', colors(k,:), 'FontWeight', 'bold', ...
        'HorizontalAlignment', 'center', 'VerticalAlignment', 'top', ...
        'Clipping', 'on');

    % OR label: angle + instantaneous reduction ratio (2 decimals)
    [~, iUL] = max(y_or - 0.02*x_or);  % bias toward upper-left boundary
    text(ax, x_or(iUL) + or_dx, y_or(iUL) + or_dy, ...
        sprintf('OR%d (%.1f°, r=%.2f)', keep_idx(k), elbow_deg, ratio_current), ...
        'Color', colors(k,:), 'FontWeight', 'bold', ...
        'HorizontalAlignment', 'right', 'VerticalAlignment', 'bottom', ...
        'Clipping', 'on');
end

hold(ax, 'off');

exportgraphics(fig, "C:\Users\railab-chh\Downloads\DEMO_pitching-OR.pdf", "ContentType", "vector");


%% Configs (kept separate because they live in different CSV logs)
start_idx = 31200;
end_idx = 31750;

cfg_speed = struct();
cfg_speed.log.name = "raipal_2026-02-16-05-00-10_pitching-demo_x1.25_20.6_TIP_";
cfg_speed.log.startIdx = start_idx;
cfg_speed.log.endIdx = end_idx;
cfg_speed.plot.title   = "Tip Speed & Wrist Velocity";
cfg_speed.plot.x       = {'time'};
cfg_speed.plot.xLabel  = 'Time [s]';
cfg_speed.plot.y       = {'speed'};
cfg_speed.plot.yLabel  = 'Translational Speed [m/s]';
cfg_speed.plot.legends = {'Tip Spd.'};

cfg_pitch = struct();
cfg_pitch.log.name = "raipal_2026-02-16-05-00-10_pitching-demo_x1.25_20.6_";
cfg_pitch.log.startIdx = start_idx;
cfg_pitch.log.endIdx = end_idx;
cfg_pitch.plot.x       = {'time'};
cfg_pitch.plot.y       = {'genVelocity_16'};
cfg_pitch.plot.yLabel  = 'Angular Velocity [rad/s]';
cfg_pitch.plot.legends = {'Wrist Ang. Vel.'};

cfg_plot = struct();
cfg_plot.markers = [39.2, 39.45, 39.55, 39.61, 39.6415, 39.76];

%% Load both logs
log_speed = csvToDictionary("log_data\raw\" + cfg_speed.log.name + ".csv", cfg_speed.log.startIdx, cfg_speed.log.endIdx);
log_pitch = csvToDictionary("log_data\raw\" + cfg_pitch.log.name + ".csv", cfg_pitch.log.startIdx, cfg_pitch.log.endIdx);

% X (time) for each log
t_speed = log_speed('time')/1e6;
t_pitch = log_pitch('time')/1e6;

% Y series
v_tip   = log_speed('speed');
w_pitch = log_pitch('genVelocity_16');

%% Dual-axis plot
fig = figure('Color','w');
set(fig, 'Units', 'centimeters', 'Position', fig_size);
ax = axes(fig); %#ok<LAXES>
hold(ax, 'on');
grid(ax, 'on');

% Plot left axis
yyaxis(ax, 'left');
p1 = plot(ax, t_speed, v_tip, 'DisplayName', cfg_speed.plot.legends{1});
ylabel(ax, cfg_speed.plot.yLabel);

% Plot right axis
yyaxis(ax, 'right');
p2 = plot(ax, t_pitch, w_pitch, 'DisplayName', cfg_pitch.plot.legends{1});
ylabel(ax, cfg_pitch.plot.yLabel);

% Common x-axis + title
xlabel(ax, cfg_speed.plot.xLabel);
title(ax, cfg_speed.plot.title);

% x-limits that cover both logs
xlim(ax, [min([t_speed; t_pitch]) max([t_speed; t_pitch])]);

% --- Align y=0 on both y-axes ---
yyaxis(ax,'left');  ylL = ylim(ax);
yyaxis(ax,'right'); ylR = ylim(ax);

if ylL(1) > 0, ylL(1) = 0; end
if ylL(2) < 0, ylL(2) = 0; end
if ylR(1) > 0, ylR(1) = 0; end
if ylR(2) < 0, ylR(2) = 0; end

pL = (0 - ylL(1)) / (ylL(2) - ylL(1) + eps);
pR = (0 - ylR(1)) / (ylR(2) - ylR(1) + eps);
p = min(max(max(pL, pR), 0.1), 0.9);

yyaxis(ax,'left');
topL = max(ylL(2), 0);
botL = min(ylL(1), 0);
scaleL = max( topL/(1-p+eps), (-botL)/(p+eps) );
ylim(ax, [-p*scaleL, (1-p)*scaleL]);

yyaxis(ax,'right');
topR = max(ylR(2), 0);
botR = min(ylR(1), 0);
scaleR = max( topR/(1-p+eps), (-botR)/(p+eps) );
ylim(ax, [-p*scaleR, (1-p)*scaleR]);

y_tick_step_L = 5;
y_tick_step_R = 10;

% --- y-ticks: 1 decimal place ---
yyaxis(ax,'left');
y_limits = ylim(ax);
y_ticks_base = floor(y_limits(1)/y_tick_step_L)*y_tick_step_L : y_tick_step_L : floor(y_limits(2)/y_tick_step_L)*y_tick_step_L;
y_ticks = sort(unique([y_ticks_base, max(v_tip)]));
y_ticks = round(y_ticks, 1);
if numel(y_ticks) >= 2
  y_ticks(end-1) = [];
end
yticks(ax, y_ticks);
yticklabels(ax, compose('%.1f', y_ticks));   % <-- 1 decimal

yyaxis(ax,'right');
y_limits = ylim(ax);
y_ticks_base = floor(y_limits(1)/y_tick_step_R)*y_tick_step_R : y_tick_step_R : floor(y_limits(2)/y_tick_step_R)*y_tick_step_R;
y_ticks = sort(unique([y_ticks_base, max(w_pitch)]));
y_ticks = round(y_ticks, 1);
if numel(y_ticks) >= 2
  y_ticks(end-1) = [];
end
yticks(ax, y_ticks);
yticklabels(ax, compose('%.1f', y_ticks));   % <-- 1 decimal

% --- Add a right-axis horizontal line at max(w_pitch) ---
yyaxis(ax,'right');
yline(ax, max(w_pitch), '-', 'Color', [0.5 0.5 0.5], 'LineWidth', 0.75, 'HandleVisibility','off');

% --- Vertical markers ---
if isfield(cfg_plot,'markers') && ~isempty(cfg_plot.markers)
  m = sort(cfg_plot.markers(:)', 'ascend');
  xl = xlim(ax);

  yyaxis(ax,'left');
  yl = ylim(ax);
  yText = yl(2) - 0.02*(yl(2)-yl(1));
  xTextOffset = 0.002 * (xl(2)-xl(1));

  for k = 1:numel(m)
    xmk = m(k);
    xline(ax, xmk, 'k--', 'LineWidth', 1.5, 'HandleVisibility','off');
    text(ax, xmk + xTextOffset, yText, sprintf('%d', k), ...
      'Color','k', 'FontWeight','bold', ...
      'HorizontalAlignment','left', 'VerticalAlignment','top', ...
      'Clipping','on', 'HandleVisibility','off');
  end
end

% One combined legend
lgd = legend(ax, [p1 p2], {cfg_speed.plot.legends{1}, cfg_pitch.plot.legends{1}}, 'Location', 'southwest');
lgd.BackgroundAlpha = 0.85;

hold(ax, 'off');

exportgraphics(fig, "C:\Users\railab-chh\Downloads\DEMO_pitching-speed.pdf", "Resolution", 300);