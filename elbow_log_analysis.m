% Coefficients from CFB_approx polynomial fit
coeff_gc_3to5 = {
  7.487348861572753, -88.17268027058002, 470.7676589474978, -1507.949571146181, ...
  3232.924175366478, -4903.587201465159, 5428.248858581333, -4464.681332539945, ...
  2751.745273804874, -1265.589676288956, 418.4687974259759, -86.36058832217348, ...
  6.540591054554072, -1.820824957906633, 2.514881046035373, 1.022698355939166,  ...
  0.000002629580961635315};

coeff_gc_3to4 = {
  -4.599352352699746, 54.74231689702444, -295.2311950384566, 953.4731854791298, ...
  -2052.475077842654, 3099.83081925634, -3363.608943070999, 2636.542385991094,  ...
  -1477.728449862623, 579.0562233124353, -157.9024440428828, 35.05094543077983, ...
  -7.773109367008759, -0.656408504334005, 0.798469056970581, 2.098468226080274, ...
  -0.000001063969759856979};


coeff_gv_3to5 = {
  330.1514673459376, -3939.920472235199, 21319.39656007634, -69182.34049371608, ...
  150034.7654546927, -229354.9209106861, 253935.6748254774, -205870.6560746802, ...
  121976.630364867, -52021.60108263157, 15489.35547273981, -3096.107249114025,  ...
  433.0072547668175, -57.21741621150876, -0.9670218724995983, 4.921276349693374,...
  1.022525952484172};

coeff_gv_3to4 = {
  -160.5479994847285, 1895.920138547342, -10134.44134516569, 32423.08074300719, ...
  -69176.96131680353, 103833.9828875482, -112772.7868033985, 89865.6754559566,  ...
  -52817.85947080839, 22895.36338178862, -7274.853311104497, 1636.804103542675, ...
  -218.7685080450584, 6.311919990384402, -3.93273889555925, 1.64496743133799,   ...
  2.097821834064398};

coeff_gc_5to3= {
  -0.01820207462949832, 0.3515207083714699, -3.095862950743284, 16.46669621013979, ...
  -59.05557308462212, 151.014385452334, -284.0127408330851, 399.5792368412747, ...
  -423.9498339027933, 339.9044621761164, -205.6485389436184, 93.83159046256269,...
  -32.59323031770014, 9.041352427928961, -2.195749936345363, 0.9751610264689757,...
  0.000008330739145766052};


num_points = 4000;
a_domain = linspace(0,89.13041976,num_points)' * pi / 180 ;

%% Global font scaling
set(groot, "defaultAxesFontSize",   25);
set(groot, "defaultTextFontSize",   30);
set(groot, "defaultLegendFontSize", 20);
set(groot, "defaultColorbarFontSize", 20);
set(groot, 'DefaultLineLineWidth', 3);
%%
load('results/v1-5_result_2026-02-09.mat')
load('MOR/MOR_MHR.mat')
%% Plot effective reduction ratio using polynomial coefficients
% b(a) and wb(a) polynomial fits (highest to lowest degree)
b_fit = polyval(cell2mat(coeff_gc_3to5), a_domain);
wb_fit = polyval(cell2mat(coeff_gv_3to5), a_domain);

inv_wb_fit = 1 ./ wb_fit;
normalized_inv_wb_fit = inv_wb_fit / max(inv_wb_fit(300:num_points,:));

if ishandle(99)
  close(99);
end
fig_reduction = figure(99);
title("RAIPAL Crossed 4-bar Effective Reduction Ratio")
xlabel("Elbow Angle [deg]")
ylabel("Normalized Torque Outupt")
xlim([15 max(b_fit * 180 / pi)]);
ylim([0 1]);
hold on;
legend;
plot(b_fit * 180 / pi, inv_wb_fit, "DisplayName", "RAIPAL Crossed 4-Bar");
hold off;

% exportgraphics(fig_reduction, "results/FIGURE_reduction_from_coeffs.fig", "Resolution", 300)
%% 
ratio_min = 0.5311;
ratio_max = 0.875;

%% A simple movie about changing MOR (single plot)
idx_start = length(a_domain);
idx_end = 100;
frame_step = -15;
save_movie = true;
movie_file = "videos/MOR_animation.avi";

if ishandle(1)
  close(1);
end

fig_ts = figure(1);
target_w = 1280;
target_h = 720;
set(fig_ts, "Units", "pixels", "Position", [100 100 target_w target_h]);
set(fig_ts, "Resize", "off");

if save_movie
  v = VideoWriter(movie_file, "Motion JPEG AVI");
  v.FrameRate = 60;
  v.Quality = 100;
  open(v);
end

x_limits = [-30 30];
y_limits = [-190 190];

for idx = idx_start:frame_step:idx_end
  clf(fig_ts);

  hold on;
  title("Operating Region (Output)");
  xlabel("Speed [rad/s]");
  ylabel("Torque [Nm]");
  plot(MOR_keypoints_MHR_peak(:,1) / ratio_min, MOR_keypoints_MHR_peak(:,2) * ratio_min,"r--", "DisplayName", "O.R. (max. speed)");
  plot(MOR_keypoints_MHR_peak(:,1) / ratio_max, MOR_keypoints_MHR_peak(:,2) * ratio_max,"g--", "DisplayName", "O.R. (max. torque)");
  plot(MOR_keypoints_MHR_peak(:,1), MOR_keypoints_MHR_peak(:,2),"b--", "DisplayName", "O.R. (input)");
  plot(MOR_keypoints_MHR_peak(:,1) / inv_wb_fit(idx), MOR_keypoints_MHR_peak(:,2) * inv_wb_fit(idx),"k-", "DisplayName", sprintf("O.R. (gc = %.3f)",b_fit(idx)));
  xlim(x_limits);
  ylim(y_limits);
  legend;
  hold off;

  drawnow;
  if save_movie
    frame = getframe(fig_ts);
    frame_img = frame.cdata;
    % h = size(frame_img,1);
    % w = size(frame_img,2);
    % if h ~= target_h || w ~= target_w
    %   resized = zeros(target_h, target_w, 3, "uint8");
    %   h_min = min(h, target_h);
    %   w_min = min(w, target_w);
    %   resized(1:h_min,1:w_min,:) = frame_img(1:h_min,1:w_min,:);
    %   frame.cdata = resized;
    % end
    writeVideo(v, frame);
  end
end

if save_movie
  close(v);
end

%%
filename="raipal_2026-02-19-00-32-43_strongman-sim_55kg_clamped-clean";
% filename="elbow_/raipal_2026-02-11-14-04-39_payload_20kg";
% filename="elbow_/raipal_2026-02-16-05-00-10_pitching-demo_x1.25_20.6_elbow";
load("log_data/elbow_/" + filename)
load("MOR/MOR_MHR")

da = a_domain(2) - a_domain(1);

%% calculate log_gf_output
log_gf_output = zeros(length(log_gf_input),1);

for i = 1:length(log_gf_input)
    ratio_current = inv_wb_fit(round(log_gc_input(i) / da));
%     log_gc_input(i) - a_domain(round(log_gc_input(i) / da))
    log_gf_output(i) = log_gf_input(i) * ratio_current;
end

%% Animating torque-speed curve with position-dependent MOR
idx_start = 1;
idx_end = length(log_time);
frame_step = 10;
% frame_step = 2;
save_movie = true;
movie_file = "videos/ts_"+filename+".avi";

if ishandle(1)
  close(1);
end

fig_ts = figure(1);
target_w = 1280;
target_h = 720;
set(fig_ts, "Units", "pixels", "Position", [100 100 target_w target_h]);
set(fig_ts, "Resize", "off");

if save_movie
  v = VideoWriter(movie_file, "Motion JPEG AVI");
  v.FrameRate = 60;
  v.Quality = 100;
  open(v);
end

x_limits = [-30 30];
y_limits = [-190 190];

for idx = idx_start:frame_step:idx_end
  clf(fig_ts);

  subplot(1,2,1);
  hold on;
  title("Torque-Speed (Input)");
  xlabel("Speed [rad/s]");
  ylabel("Torque [Nm]");
  plot(log_gv_input, log_gf_input, "DisplayName", "Input");
  plot(MOR_keypoints_MHR_peak(:,1), MOR_keypoints_MHR_peak(:,2), "DisplayName", "O.R. (input)");
  scatter([log_gv_input(idx)], [log_gf_input(idx)],"MarkerEdgeColor","red","LineWidth",2.0, "DisplayName", sprintf("t = %.3f",log_time(idx)))
  xlim(x_limits);
  ylim(y_limits);
  lgd = legend('location','southwest');
  lgd.BackgroundAlpha = 0.95;
  hold off;

  subplot(1,2,2);
  hold on;
  title("Torque-Speed (Output)");
  xlabel("Speed [rad/s]");
  ylabel("Torque [Nm]");
  plot(log_gv_output, log_gf_output, "DisplayName", "Output");
  plot(MOR_keypoints_MHR_peak(:,1) / ratio_min, MOR_keypoints_MHR_peak(:,2) * ratio_min,"r--", "DisplayName", "O.R. (max. speed)");
  plot(MOR_keypoints_MHR_peak(:,1) / ratio_max, MOR_keypoints_MHR_peak(:,2) * ratio_max,"g--", "DisplayName", "O.R. (max. torque)");
  plot(MOR_keypoints_MHR_peak(:,1), MOR_keypoints_MHR_peak(:,2),"b--", "DisplayName", "O.R. (input)");
%   ratio_current = log_gv_input(idx) / log_gv_output(idx);
  ratio_current = inv_wb_fit(round(log_gc_input(idx) / da));
  plot(MOR_keypoints_MHR_peak(:,1) / ratio_current, MOR_keypoints_MHR_peak(:,2) * ratio_current,"k-", "DisplayName", sprintf("O.R. (gc = %.3f)",log_gc_output(idx)));
  scatter([log_gv_output(idx)], [log_gf_output(idx)],"MarkerEdgeColor","red","LineWidth",2.0, "DisplayName", sprintf("t = %.3f",log_time(idx)))
  xlim(x_limits);
  ylim(y_limits);
  lgd = legend('location','southwest');
  lgd.BackgroundAlpha = 0.95;
  hold off;

  drawnow;
  if save_movie
    frame = getframe(fig_ts);
    frame_img = frame.cdata;
    % h = size(frame_img,1);
    % w = size(frame_img,2);
    % if h ~= target_h || w ~= target_w
    %   disp("warning! size descrepancy.")
    %   fprintf("h: %d\n",h)
    %   fprintf("w: %d\n",w)
    %   resized = zeros(target_h, target_w, 3, "uint8");
    %   h_min = min(h, target_h);
    %   w_min = min(w, target_w);
    %   resized(1:h_min,1:w_min,:) = frame_img(1:h_min,1:w_min,:);
    %   frame.cdata = resized;
    % end
    writeVideo(v, frame);
  end
end

if save_movie
  close(v);
end