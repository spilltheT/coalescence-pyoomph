close; clear; clc;

tmp = cell(1,8);

for k = 1:8
    file_path = sprintf('h_f_vals/min_h_values_%d.txt', k);
    T = readtable(file_path, 'Delimiter', '\t');
    tmp{k} = T.Min_h;         
end

[min_h, min_h1, min_h2, min_h3, min_h4, min_h5, min_h6, min_h7] = tmp{:};

t1 = linspace(0,1000,10001)';
t2 = linspace(0,250,2501)';

t1 = +0.39 + t1;   % t_0 correction

beta  = [0.8 0.8 0.5 0.5  0.8 0.8 0.5 0.5];
Pe    = [1 10 1 10  1 10 1 10];

theta1 = 10*pi/180;
theta2 = 20*pi/180;

theta = [theta1 theta1 theta1 theta1  theta2 theta2 theta2 theta2];

Gamma_0 = 0.8;

h_a = (1 - beta * Gamma_0/2) .* theta.^4;

blue1 = [173, 216, 230]/255;  % LightBlue
blue2 = [100, 149, 237]/255;  % CornflowerBlue
blue3 = [65, 105, 225]/255;   % RoyalBlue
blue4 = [0, 0, 139]/255;

green1 = [144, 238, 144]/255;  % LightGreen
green2 = [60, 179, 113]/255;   % MediumSeaGreen
green3 = [34, 139, 34]/255;    % ForestGreen
green4 = [0, 100, 0]/255;

blue_shades = {blue1, blue2, blue3, blue4};
green_shades = {green1, green2, green3, green4};

set(groot,'defaultAxesTickLabelInterpreter','latex'); 

% Uncompensated
loglog(t1(1:end), min_h(1:end),'o','MarkerSize',13,'MarkerFaceColor',blue_shades{1}, 'MarkerEdgeColor','k');
hold on;
loglog(t1(1:end),min_h1(1:end),'o','MarkerSize',13,'MarkerFaceColor',blue_shades{2}, 'MarkerEdgeColor','k');
loglog(t1(1:end),min_h2(1:end),'o','MarkerSize',13,'MarkerFaceColor',blue_shades{3}, 'MarkerEdgeColor','k');
loglog(t1(1:end),min_h3(1:end),'o','MarkerSize',13,'MarkerFaceColor',blue_shades{4}, 'MarkerEdgeColor','k');
loglog(t2(1:end),min_h4(1:end),'o','MarkerSize',13,'MarkerFaceColor',green_shades{1}, 'MarkerEdgeColor','k');
loglog(t2(1:end),min_h5(1:end),'o','MarkerSize',13,'MarkerFaceColor',green_shades{2}, 'MarkerEdgeColor','k');
loglog(t2(1:end),min_h6(1:end),'o','MarkerSize',13,'MarkerFaceColor',green_shades{3}, 'MarkerEdgeColor','k');
loglog(t2(1:end),min_h7(1:end),'o','MarkerSize',13,'MarkerFaceColor',green_shades{4}, 'MarkerEdgeColor','k');
loglog(t1(10:110),0.000872*t1(10:110).^(1), 'k', LineWidth=4)

% Compensated
% loglog(t1(1:end), min_h(1:end)./h_a(1),'o','MarkerSize',13,'MarkerFaceColor',blue_shades{1}, 'MarkerEdgeColor','k');
% hold on;
% loglog(t1(1:end),min_h1(1:end)./h_a(2),'o','MarkerSize',13,'MarkerFaceColor',blue_shades{2}, 'MarkerEdgeColor','k');
% loglog(t1(1:end),min_h2(1:end)./h_a(3),'o','MarkerSize',10,'MarkerFaceColor',blue_shades{3}, 'MarkerEdgeColor','k');
% loglog(t1(1:end),min_h3(1:end)./h_a(4),'o','MarkerSize',13,'MarkerFaceColor',blue_shades{4}, 'MarkerEdgeColor','k');
% loglog(t2(1:end),min_h4(1:end)./h_a(5),'o','MarkerSize',13,'MarkerFaceColor',green_shades{1}, 'MarkerEdgeColor','k');
% loglog(t2(1:end),min_h5(1:end)./h_a(6),'o','MarkerSize',13,'MarkerFaceColor',green_shades{2}, 'MarkerEdgeColor','k');
% loglog(t2(1:end),min_h6(1:end)./h_a(7),'o','MarkerSize',13,'MarkerFaceColor',green_shades{3}, 'MarkerEdgeColor','k');
% loglog(t2(1:end),min_h7(1:end)./h_a(8),'o','MarkerSize',10,'MarkerFaceColor',green_shades{4}, 'MarkerEdgeColor','k');
% loglog(t1(2:700),0.272*t1(2:700).^(1), 'r', LineWidth=4)

ax = gca;
ax.FontSize = 20;
ax.TickLength = [0.015, 0.015];
ax.LineWidth = 2;