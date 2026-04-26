%%
cfg = struct();
cfg.log.name = "raipal_2026-04-01-14-48-18_20kg-LEFT";
cfg.log.startIdx = 28000;
cfg.log.endIdx = 36000;

cfg.plot.title = "Upper Arm Torque";
% cfg.plot.x = {'time','time','time','time'};
cfg.plot.x = {'index','index','index','index'};
cfg.plot.xLabel = 'Time [s]';
cfg.plot.y = {'actualTorque_7','actualTorque_8','actualTorque_9','actualTorque_10'};
cfg.plot.yLabel = 'Torque [Nm]';
cfg.plot.legends = {'\tau_{0L} (Flexion)','\tau_{1L} (Adduction)','\tau_{2L} (Rotation)','\tau_{3L} (Elbow)'};

%%
cfg = struct();
cfg.log.name = "raipal_2026-02-20-13-10-21_fore-demo_x1.2_asym";
cfg.log.startIdx = 26800;
cfg.log.endIdx = 27300;

cfg.plot.title = "Forearm Angular Velocity";
cfg.plot.x = {'time','time','time'};
cfg.plot.xLabel = 'Time [s]';
cfg.plot.y = {'genVelocity_16','genVelocity_15','genVelocity_17'};
cfg.plot.yLabel = 'Angular Velocity [rad/s]';
cfg.plot.legends = {'\omega_{15} (Roll)','\omega_{16} (Pitch)','\omega_{17} (Yaw)'};

%% Endpoint Speed
cfg = struct();
cfg.log.name = "raipal_2026-02-16-05-00-10_pitching-demo_x1.25_20.6_TIP_";
cfg.log.startIdx = 31000;
cfg.log.endIdx = 32000;

cfg.plot.title = "Tip Speed";
cfg.plot.x = {'time'};
cfg.plot.xLabel = 'Time [s]';
cfg.plot.y = {'speed'};
cfg.plot.yLabel = 'Translational Speed [m/s]';
cfg.plot.legends = {'v_{L} (Left Tip)'};

%% Pitch velocity
cfg = struct();
cfg.log.name = "raipal_2026-02-16-05-00-10_pitching-demo_x1.25_20.6_";
cfg.log.startIdx = 31000;
cfg.log.endIdx = 32000;

cfg.plot.title = "Pitch Axis Velocity";
cfg.plot.x = {'time'};
cfg.plot.xLabel = 'Time [s]';
cfg.plot.y = {'genVelocity_16'};
cfg.plot.yLabel = 'Angular Velocity [rad/s]';
cfg.plot.legends = {'\omega_{16} (Roll)'};

%%
set(groot, "defaultAxesFontSize",   25);
set(groot, "defaultTextFontSize",   30);
set(groot, "defaultLegendFontSize", 20);
set(groot, "defaultColorbarFontSize", 20);
set(groot, 'DefaultLineLineWidth', 3);
%% Actual plotting part
log_data = csvToDictionary("log_data\raw\" + cfg.log.name + ".csv", cfg.log.startIdx, cfg.log.endIdx);

if numel(cfg.plot.x) ~= numel(cfg.plot.y) || numel(cfg.plot.y) ~= numel(cfg.plot.legends)
	error('plot_x, plot_y, and plot_legends must have the same length.');
end

fig = figure;
hold on;

x_limits = [0 0];
y_limits = [0 0];

for i = 1:numel(cfg.plot.y)
	x_key = cfg.plot.x{i};
	y_key = cfg.plot.y{i};

	if ~isKey(log_data, x_key) && ~strcmp(x_key,'index')
		error('Missing x key in log_data: %s', x_key);
	end
	if ~isKey(log_data, y_key)
		error('Missing y key in log_data: %s', y_key);
	end
    
    if strcmp(x_key, 'index')
        x_vals = (cfg.log.startIdx:1:cfg.log.endIdx)';
    elseif strcmp(x_key, 'time')
        x_vals = log_data(x_key)/1e6;
    else
    	x_vals = log_data(x_key);
    end

    y_vals = log_data(y_key);

    if i == 1
        x_limits = [min(x_vals) max(x_vals)];
        y_limits = [min(y_vals) max(y_vals)];
    else
        x_limits = [min(min(x_vals),x_limits(1)) max(max(x_vals),x_limits(2))];
        y_limits = [min(min(y_vals),y_limits(1)) max(max(y_vals),y_limits(2))];
    end

	plot(x_vals, y_vals, 'DisplayName', cfg.plot.legends{i});
    
end

xlim(x_limits)
ylim(y_limits)

% Ensure grid lines include exact min/max y-limits
y_tick_step = 50;
y_ticks_base = floor(y_limits(1)/y_tick_step)*y_tick_step : y_tick_step : ceil(y_limits(2)/y_tick_step)*y_tick_step;
y_ticks = sort(unique([y_ticks_base, y_limits(1), y_limits(2)]));
yticks(y_ticks);

xlabel(cfg.plot.xLabel);
ylabel(cfg.plot.yLabel);
title(cfg.plot.title);

axis_handle = gca;
axis_handle.YLabel.Units = 'normalized';
axis_handle.YLabel.Position = [-0.15 0.5 0];
axis_handle.YLabel.VerticalAlignment = 'middle';

lgd = legend('Location', 'best');
lgd.BackgroundAlpha = 0.85;
grid on;
hold off;

exportgraphics(fig, "log_data/plot/" + cfg.log.name + "_[" + cfg.plot.title + "].png", "Resolution", 300)