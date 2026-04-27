%% Global font scaling (for paper figures)
set(groot, "defaultAxesFontSize",   12);
set(groot, "defaultTextFontSize",   15);
set(groot, "defaultLegendFontSize", 12);
set(groot, "defaultColorbarFontSize", 12);
set(groot, 'DefaultLineLineWidth', 2);

%% 1) Upper Arm Torque
cfg = struct();
cfg.log.name = "raipal_2026-04-01-14-48-18_20kg-LEFT";
cfg.log.startIdx = 32000;
cfg.log.endIdx = 38000;

cfg.plot.title   = "Upper Arm Torque";
cfg.plot.x       = {'time','time','time','time'};
cfg.plot.xLabel  = 'Time [s]';
cfg.plot.y       = {'actualTorque_7','actualTorque_8','actualTorque_9','actualTorque_10'};
cfg.plot.yLabel  = 'Torque [Nm]';
cfg.plot.legends = {'\tau_{0L} (Flexion)','\tau_{1L} (Adduction)','\tau_{2L} (Rotation)','\tau_{3L} (Elbow)'};
cfg.plot.legendPlacement = 'southwest';

cfg.markers = [76.7,82.4,85.0,86.5];

% plot
fig1 = plotLogSeries(cfg);

% export
exportgraphics(fig1, "C:\Users\railab-chh\Downloads\DEMO_20kg-torque.pdf");

%% 2) Forearm Angular Velocity
cfg = struct();
cfg.log.name = "raipal_2026-02-20-13-10-21_fore-demo_x1.2_asym";
cfg.log.startIdx = 26950;
cfg.log.endIdx = 27250;

cfg.plot.title   = "Forearm Angular Velocity";
cfg.plot.x       = {'index','index','index'};
cfg.plot.xLabel  = 'Time [s]';
cfg.plot.y       = {'genVelocity_16','genVelocity_15','genVelocity_17'};
cfg.plot.yLabel  = 'Angular Velocity [rad/s]';
cfg.plot.legends = {'\omega_{5L} (Roll)','\omega_{6L} (Pitch)','\omega_{7L} (Yaw)'};
cfg.plot.legendPlacement = 'southwest';


cfg.filename = "forearm_angular_velocity";

fig2 = plotLogSeries(cfg);

exportgraphics(fig2, "C:\Users\railab-chh\Downloads\DEMO_forearm-speed.pdf");

