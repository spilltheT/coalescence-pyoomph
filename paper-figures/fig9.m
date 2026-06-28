base_path = '9029/';

time_steps = [0, 40, 100]*10;

x_values = cell(length(time_steps), 1);
h_values = cell(length(time_steps), 1);
Gamma_values = cell(length(time_steps), 1);

for i = 1:length(time_steps)
    time_step = time_steps(i);

    tt = time_step*0.1;
    
    file_path = sprintf('%sdomain_%06d.txt', base_path, time_step);
    
    data = readtable(file_path, 'Delimiter', '\t');
    
    x = data{:,1}; 
    h = data{:,2};  
    Gamma = data{:,4};

    x_values{i} = x;
    h_values{i} = h;
    Gamma_values{i} = Gamma;
end

colors = lines(length(time_steps));

fig = figure;
fig.Position = [100, 100, 800, 400];

hold on;
for i = 1:length(time_steps)
%     plot(x_values{i}, h_values{i}, 'o', 'MarkerSize', 8, 'Color', colors(i,:), 'MarkerFaceColor', colors(i,:), ...
%         'DisplayName', sprintf('Time Step %d', time_steps(i)));
%     plot(x_values{i}, Gamma_values{i}, 'o', 'MarkerSize', 8, 'Color', colors(i,:), 'MarkerFaceColor', colors(i,:), ...
%         'DisplayName', sprintf('Time Step %d', time_steps(i)));
    plot(x_values{i}, Gamma_values{i}, 'o', 'MarkerSize', 8, 'Color', 'k', 'MarkerFaceColor', 'k', ...
        'DisplayName', sprintf('Time Step %d', time_steps(i)));
end

hold off;

ax = gca;
ax.FontSize = 24;
ax.TickLength = [0.015, 0.015];
ax.LineWidth = 2;
% ax.XTick = [];  
% ax.YTick = [];
box on;
% daspect([10 1 1]);
% axis([-2.02 2.2 0.0 0.25])
% axis([-2.02 2.2 0.0 1.0])
% axis([-0.1 2 0.0 1.0])
