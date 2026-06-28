close all; clear; clc;

T = readtable('xf_vs_Pe.csv');

Ma_vals = unique(T.Ma);

cl_colors = {
    [44 105 176]/255,    % dark blue
    [82 141 208]/255,    % mid blue
    [165 195 226]/255,   % light blue
    [255 247 153]/255,   % pale yellow
    [253 174 97]/255,    % light orange
    [234 94 56]/255,     % red-orange
    [178 24 43]/255      % deep red
};

set(groot,'defaultAxesTickLabelInterpreter','latex');

figure; hold on;

for i = 1:length(Ma_vals)
    Ma = Ma_vals(i);

    idx = T.Ma == Ma;
    Pe = T.Pe(idx);
    xf = T.x_f(idx);

    loglog(Pe, xf, 'o', ...
        'MarkerSize', 12, ...
        'MarkerFaceColor', cl_colors{i}, ...
        'MarkerEdgeColor', 'k', ...
        'LineWidth', 2);
end

Pe_ref = [1 2 3 5];
loglog(Pe_ref, 0.15 * Pe_ref.^1, 'k', 'LineWidth', 2);

ax = gca;
ax.XScale = 'log';
ax.YScale = 'log';
ax.FontSize = 20;
ax.TickLength = [0.015, 0.015];
ax.LineWidth = 2;
ax.XColor = [0 0 0];
ax.YColor = [0 0 0];

set(gcf,'Color','w','InvertHardcopy','off','Renderer','painters');