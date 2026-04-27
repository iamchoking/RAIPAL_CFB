function fig = plotLogSeries(cfg, varargin)
% plotLogSeries  Plot time series from one log CSV using a cfg struct.
%
% Required cfg fields:
%   cfg.log.name, cfg.log.startIdx, cfg.log.endIdx
%   cfg.plot.title, cfg.plot.x, cfg.plot.xLabel, cfg.plot.y, cfg.plot.yLabel, cfg.plot.legends
%
% Optional cfg fields:
%   cfg.width_cm   (default: 16)
%   cfg.height_cm  (default: 9)
%   cfg.markers    numeric vector of x-values where vertical lines are drawn
%
% Optional plot fields:
%   cfg.plot.legendPlacement (default: "best")  % any valid legend Location value
%
% Name-value options:
%   "InDir"           (default: "log_data\raw\")
%   "TightPad"        (default: 0.01)  extra padding added to TightInset (all sides)
%   "LeftInsetScale"  (default: 0.80)  shrink left margin relative to TightInset

opts = struct();
opts.InDir = "log_data\raw\";
opts.CloseFig = false;
opts.TightPad = 0.01;
opts.LeftInsetScale = 0.80;

if ~isempty(varargin)
  if mod(numel(varargin), 2) ~= 0
    error('Name-value arguments must come in pairs.');
  end
  for k = 1:2:numel(varargin)
    opts.(string(varargin{k})) = varargin{k+1};
  end
end

% --- Figure sizing defaults ---
if ~isfield(cfg, "width_cm")  || isempty(cfg.width_cm),  cfg.width_cm  = 16; end
if ~isfield(cfg, "height_cm") || isempty(cfg.height_cm), cfg.height_cm = 9;  end

% --- Legend placement default ---
if ~isfield(cfg, "plot") || ~isfield(cfg.plot, "legendPlacement") || isempty(cfg.plot.legendPlacement)
  cfg.plot.legendPlacement = "best";
end

hasMarkers = isfield(cfg, "markers") && ~isempty(cfg.markers);

% --- Load data ---
csvPath = opts.InDir + cfg.log.name + ".csv";
log_data = csvToDictionary(csvPath, cfg.log.startIdx, cfg.log.endIdx);

% --- Validate plot configuration ---
if numel(cfg.plot.x) ~= numel(cfg.plot.y) || numel(cfg.plot.y) ~= numel(cfg.plot.legends)
  error('plot.x, plot.y, and plot.legends must have the same length.');
end

% --- Create figure with prescribed size ---
fig = figure('Color','w');
set(fig, "Units","centimeters");
pos = get(fig, "Position");     % [x y w h] in cm
pos(3) = cfg.width_cm;
pos(4) = cfg.height_cm;
set(fig, "Position", pos);

ax = axes(fig); %#ok<LAXES>
hold(ax, "on");

x_limits = [ inf, -inf ];
y_limits = [ inf, -inf ];

for i = 1:numel(cfg.plot.y)
  x_key = cfg.plot.x{i};
  y_key = cfg.plot.y{i};

  if ~strcmp(x_key,'index') && ~isKey(log_data, x_key)
    error('Missing x key in log_data: %s', x_key);
  end
  if ~isKey(log_data, y_key)
    error('Missing y key in log_data: %s', y_key);
  end

  if strcmp(x_key, 'index')
    x_vals = (cfg.log.startIdx:1:cfg.log.endIdx)';
  elseif strcmp(x_key, 'time')
    x_vals = log_data(x_key) / 1e6;
  else
    x_vals = log_data(x_key);
  end

  y_vals = log_data(y_key);

  x_limits = [ min(x_limits(1), min(x_vals)), max(x_limits(2), max(x_vals)) ];
  y_limits = [ min(y_limits(1), min(y_vals)), max(y_limits(2), max(y_vals)) ];

  plot(ax, x_vals, y_vals, 'DisplayName', cfg.plot.legends{i});
end

xlim(ax, x_limits);
ylim(ax, y_limits);

% --- Ticks (keeps your behavior) ---
y_tick_step = 50;
y_ticks_base = floor(y_limits(1)/y_tick_step)*y_tick_step : y_tick_step : ceil(y_limits(2)/y_tick_step)*y_tick_step;
y_ticks = sort(unique([y_ticks_base, y_limits(1), y_limits(2)]));
yticks(ax, y_ticks);

% --- Vertical marker lines + numbered labels at TOP ---
if hasMarkers
  m = sort(cfg.markers(:)', 'ascend');

  y0 = y_limits(1);
  y1 = y_limits(2);
  yr = y1 - y0;
  yText = y1 - 0.02 * yr;
  xTextOffset = 0.002 * (x_limits(2) - x_limits(1));

  for k = 1:numel(m)
    xmk = m(k);
    xline(ax, xmk, 'k--', 'LineWidth', 1.5, 'HandleVisibility','off');
    text(ax, xmk + xTextOffset, yText, sprintf('%d', k), ...
      'Color','k', 'FontWeight','bold', ...
      'HorizontalAlignment','left', 'VerticalAlignment','top', ...
      'Clipping','on', 'HandleVisibility','off');
  end
end

% --- Labels / title ---
xlabel(ax, cfg.plot.xLabel);
ylabel(ax, cfg.plot.yLabel);
title(ax, cfg.plot.title);

% Optional label spacing controls (if available)
if isprop(ax, "LabelSpacing")
  ax.LabelSpacing = 0.2;
end
if isprop(ax, "TickLabelSpacing")
  ax.TickLabelSpacing = 0.2;
end

% --- Legend / grid ---
lgd = legend(ax, 'Location', char(cfg.plot.legendPlacement));
lgd.BackgroundAlpha = 0.85;
grid(ax, 'on');
hold(ax, "off");

% --- Tighten outside margins, then reduce LEFT margin ---
drawnow;
ax.Units = 'normalized';
tight = ax.TightInset;   % [left bottom right top]

pad = [opts.TightPad opts.TightPad opts.TightPad opts.TightPad];
li = tight + pad;

% Shrink left inset, but keep a small safety minimum
li(1) = max(0.005, li(1) * opts.LeftInsetScale);

ax.LooseInset = li;

end