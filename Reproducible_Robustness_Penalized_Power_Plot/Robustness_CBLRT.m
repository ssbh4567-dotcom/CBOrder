% Small sample size
N_small = [0.0476, 0.1316, 0.1684, 0.2068, 0.2474, 0.2954, 0.3388, 0.3764, 0.4192, 0.4554, 0.4898, 0.5136, 0.5302, 0.5414, 0.5514, 0.5580, 0.5672, 0.5762];
Mix1_small = [0.0438, 0.1316, 0.1648, 0.2018, 0.2480, 0.2838, 0.3286, 0.3628, 0.4006, 0.4372, 0.4676, 0.4880, 0.5032, 0.5164, 0.5290, 0.5408, 0.5492, 0.5500];
Mix2_small = [0.0440, 0.1378, 0.1712, 0.2112, 0.2514, 0.2976, 0.3388, 0.3834, 0.4150, 0.4396, 0.4692, 0.4912, 0.5072, 0.5200, 0.5328, 0.5476, 0.5516, 0.5610];
Mix3_small = [0.0464, 0.1450, 0.1840, 0.2222, 0.2638, 0.3036, 0.3368, 0.3756, 0.4090, 0.4480, 0.4700, 0.4916, 0.5100, 0.5316, 0.5406, 0.5506, 0.5540, 0.5606];
Mix4_small = [0.0440, 0.1488, 0.1898, 0.2320, 0.2750, 0.3152, 0.3508, 0.3840, 0.4212, 0.4578, 0.4806, 0.5018, 0.5186, 0.5320, 0.5458, 0.5540, 0.5590, 0.5672];
laplace_small = [0.0520, 0.1510, 0.2014, 0.2484, 0.2976, 0.3458, 0.3980, 0.4404, 0.4760, 0.5108, 0.5344, 0.5594, 0.5726, 0.5856, 0.5926, 0.5970, 0.6044, 0.6082];

% Moderate sample size
N_Mod = [0.0506, 0.2960, 0.3864, 0.4782, 0.5614, 0.6198, 0.6694, 0.7050, 0.7394, 0.7614, 0.7850, 0.8124, 0.8362, 0.8582, 0.8756, 0.8960, 0.9128, 0.9276];
Mix1_Mod = [0.0478, 0.2824, 0.3784, 0.4636, 0.5426, 0.6056, 0.6544, 0.6912, 0.7218, 0.7518, 0.7770, 0.7960, 0.8200, 0.8346, 0.8594, 0.8784, 0.8976, 0.9148];
Mix2_Mod = [0.0480, 0.2920, 0.3882, 0.4752, 0.5416, 0.6064, 0.6556, 0.6856, 0.7194, 0.7414, 0.7664, 0.7864, 0.8050, 0.8192, 0.8384, 0.8570, 0.8732, 0.8916];
Mix3_Mod = [0.0506, 0.2868, 0.3784, 0.4690, 0.5416, 0.6106, 0.6580, 0.6950, 0.7222, 0.7514, 0.7760, 0.8010, 0.8238, 0.8444, 0.8650, 0.8830, 0.9040, 0.9182];
Mix4_Mod = [0.0516, 0.2990, 0.3878, 0.4744, 0.5442, 0.6118, 0.6564, 0.6912, 0.7174, 0.7430, 0.7664, 0.7904, 0.8124, 0.8328, 0.8524, 0.8708, 0.8888, 0.9070];
laplace_Mod = [0.0470, 0.3036, 0.3996, 0.4904, 0.5730, 0.6264, 0.6740, 0.7032, 0.7332, 0.7624, 0.7856, 0.8148, 0.8382, 0.8538, 0.8748, 0.8900, 0.9058, 0.9192];



N_small = [0.0476, N_small(2:end) ./ sqrt(1 + abs(1 - (0.0476 / 0.05)))];
Mix1_small = [0.0438, Mix1_small(2:end) ./ sqrt(1 + abs(1 - (0.0438 / 0.05)))];
Mix2_small = [0.0440, Mix2_small(2:end) ./ sqrt(1 + abs(1 - (0.0440 / 0.05)))];
Mix3_small = [0.0464, Mix3_small(2:end) ./ sqrt(1 + abs(1 - (0.0464 / 0.05)))];
Mix4_small = [0.0440, Mix4_small(2:end) ./ sqrt(1 + abs(1 - (0.0440 / 0.05)))];
laplace_small = [0.0520, laplace_small(2:end) ./ sqrt(1 + abs(1 - (0.0520 / 0.05)))];


N_Mod = [0.0506, N_Mod(2:end) ./ sqrt(1 + abs(1 - (0.0506 / 0.05)))];
Mix1_Mod = [0.0478, Mix1_Mod(2:end) ./ sqrt(1 + abs(1 - (0.0478 / 0.05)))];
Mix2_Mod = [0.0480, Mix2_Mod(2:end) ./ sqrt(1 + abs(1 - (0.0480 / 0.05)))];
Mix3_Mod = [0.0506, Mix3_Mod(2:end) ./ sqrt(1 + abs(1 - (0.0506 / 0.05)))];
Mix4_Mod = [0.0516, Mix4_Mod(2:end) ./ sqrt(1 + abs(1 - (0.0516 / 0.05)))];
laplace_Mod = [0.0470, laplace_Mod(2:end) ./ sqrt(1 + abs(1 - (0.0470 / 0.05)))];




% Generate sequence 'c'
u = 1:0.25:5;  % Creates 18 elements
c = [0, u];     % Prepend 0 to get 19 elements

% Create a figure with 2x2 subplots
figure;
set(gcf, 'PaperPositionMode', 'auto'); % Keeps figure size same as on-screen
set(gcf, 'Position', [100, 100, 900, 400]); % Adjust as needed
set(gca, 'FontSize', 14, 'FontWeight', 'bold');

subplot(1,2,1);
plot_power(c, N_small, Mix1_small, Mix2_small, Mix3_small, Mix4_small, laplace_small, '(e)', false);

subplot(1,2,2);
plot_power(c, N_Mod, Mix1_Mod, Mix2_Mod, Mix3_Mod, Mix4_Mod, laplace_Mod, '(f)', true);


%print('Robustness_LRT', '-depsc', '-r1200'); % EPS for LaTeX/journal submission


%% Local Function Definition
function plot_power(c, power_CBLRT, power_CBMax, power_CBMaxMin, power_CBMin, CBMinMax, power_Lap, title_str, show_legend)
    % Define fresh, contrasting colors
    color1 = [1 0 0];              
    color2 = [0 0 1]; 
    color3 = [0.4660 0.6740 0.1880]; 
    color4 = [1 0.5 1];     
    color5 = [0 1 0]; 
    color6 = [0.5 0.5 0]; 

    % Smooth line styles, minimal markers
    h1 = plot(c, power_CBLRT, ':>', 'Color', color1, 'LineWidth', 1.5, 'MarkerSize', 3); hold on;
    h2 = plot(c, power_CBMax, '<--', 'Color', color2, 'LineWidth', 1.5, 'MarkerSize', 3);
    h3 = plot(c, power_CBMaxMin, '<-', 'Color', color3, 'LineWidth', 1.5, 'MarkerSize', 3);
    h4 = plot(c, power_CBMin, '-.>', 'Color', color4, 'LineWidth', 1.5, 'MarkerSize', 3);
    h5 = plot(c, CBMinMax, '->', 'Color', color5, 'LineWidth', 1.5, 'MarkerSize', 3);
    h6 = plot(c, power_Lap, '-.>', 'Color', color6, 'LineWidth', 1.5, 'MarkerSize', 3);

    % Configure plot appearance
    xlim([0 5.1]);
    ylim([0 1.05]);
    title(title_str, 'FontWeight', 'bold', 'FontSize', 16);
    xlabel('c', 'FontSize', 14, 'FontAngle', 'italic');
    ylabel('Power', 'FontSize', 14);
    grid on; % light grid

    % Legend only if requested
    if show_legend
        legend([h1, h2, h3, h4, h5, h6], ...
            {'N', 'Mix1', 'Mix2', 'Mix3', 'Mix4', 'Lap'}, ...
            'Location', 'southeast', 'FontSize', 12, 'Box', 'off');
    end

    box on;
    hold off;
end
