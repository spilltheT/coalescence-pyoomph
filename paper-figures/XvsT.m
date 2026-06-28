close; clear; clc;

tmp = cell(1,8);

for k = 1:8
    file_path = sprintf('h_f_vals/min_h_values_%d.txt', k);
    T = readtable(file_path, 'Delimiter', '\t');
    tmp{k} = T.Min_x;         
end

[min_x, min_x1, min_x2, min_x3, min_x4, min_x5, min_x6, min_x7] = tmp{:};

t1 = linspace(0,1000,10001)';
t2 = linspace(0,250,2501)';

beta  = [0.8 0.8 0.5 0.5  0.8 0.8 0.5 0.5];
Pe    = [1 10 1 10  1 10 1 10];

theta1 = 10*pi/180;
theta2 = 20*pi/180;

theta = [theta1 theta1 theta1 theta1  theta2 theta2 theta2 theta2];

Gamma_0 = 0.8;

h_a = (1 - beta * Gamma_0/2) .* theta.^4;
x_a = (beta .* Pe.^(1/2) .* h_a .* Gamma_0)./theta;

blue = [31 119 180]./255;
orange = [255 127 14]./255;
green = [44 160 44]./255;
red = [214 39 40]./255;
purple = [148 103 189]./255;
brown = [140 86 75]./255;
pink = [227 119 194]./255;
gray = [127 127 127]./255;

cl_colors = {blue, orange, green, red, purple, brown, pink, gray};

set(groot,'defaultAxesTickLabelInterpreter','latex'); 

% Uncompensated
loglog(t1(1:3000), min_x(1:3000),'o','MarkerSize',13,'MarkerFaceColor',cl_colors{1}, 'MarkerEdgeColor','k');
hold on;
loglog(t1(1:3000),min_x1(1:3000),'o','MarkerSize',13,'MarkerFaceColor',cl_colors{2}, 'MarkerEdgeColor','k');
loglog(t1(1:3000),min_x3(1:3000),'o','MarkerSize',13,'MarkerFaceColor',cl_colors{4}, 'MarkerEdgeColor','k');
loglog(t2(1:420),min_x4(1:420),'o','MarkerSize',13,'MarkerFaceColor',cl_colors{5}, 'MarkerEdgeColor','k');
loglog(t2(1:550),min_x5(1:550),'o','MarkerSize',13,'MarkerFaceColor',cl_colors{6}, 'MarkerEdgeColor','k');
loglog(t2(1:300),min_x6(1:300),'o','MarkerSize',13,'MarkerFaceColor',cl_colors{7}, 'MarkerEdgeColor','k');
loglog(t2(5:70),0.03*t2(5:70).^(1.5), 'Color', 'k', LineWidth=5);

% Compensated
% loglog(t1(1:3000), min_x(1:3000)./x_a(1),'o','MarkerSize',13,'MarkerFaceColor',cl_colors{1}, 'MarkerEdgeColor','k');
% hold on;
% loglog(t1(1:3000),min_x1(1:3000)./x_a(2),'o','MarkerSize',13,'MarkerFaceColor',cl_colors{2}, 'MarkerEdgeColor','k');
% loglog(t1(1:3000),min_x3(1:3000)./x_a(4),'o','MarkerSize',13,'MarkerFaceColor',cl_colors{4}, 'MarkerEdgeColor','k');
% loglog(t2(1:420),min_x4(1:420)./x_a(5),'o','MarkerSize',13,'MarkerFaceColor',cl_colors{5}, 'MarkerEdgeColor','k');
% loglog(t2(1:550),min_x5(1:550)./x_a(6),'o','MarkerSize',13,'MarkerFaceColor',cl_colors{6}, 'MarkerEdgeColor','k');
% loglog(t2(1:300),min_x6(1:300)./x_a(7),'o','MarkerSize',13,'MarkerFaceColor',cl_colors{7}, 'MarkerEdgeColor','k');
% loglog(t2(2:110),0.1881*t2(2:110).^(1.5), 'Color', 'b', LineWidth=5);

ax = gca;
ax.FontSize = 20;
ax.TickLength = [0.015, 0.015];
ax.LineWidth = 2;