
% Small sample size group (N3)
N_small = [0.0460, 0.1578, 0.2046, 0.2564, 0.3072, 0.3714, 0.4358, 0.5078, ...
           0.5732, 0.6360, 0.6998, 0.7514, 0.7968, 0.8320, 0.8636, 0.8932, ...
           0.9166, 0.9366];

Mix1_small = [0.0324, 0.1164, 0.1506, 0.1918, 0.2436, 0.2944, 0.3530, 0.4124, ...
              0.4790, 0.5460, 0.6046, 0.6612, 0.7136, 0.7592, 0.8020, 0.8366, ...
              0.8694, 0.8982];

Mix2_small = [0.0318, 0.1106, 0.1438, 0.1858, 0.2334, 0.2884, 0.3434, 0.4088, ...
              0.4652, 0.5222, 0.5834, 0.6402, 0.6870, 0.7314, 0.7748, 0.8118, ...
              0.8470, 0.8752];

Mix3_small = [0.0338, 0.1236, 0.1578, 0.2046, 0.2556, 0.3154, 0.3786, 0.4428, ...
              0.5094, 0.5736, 0.6280, 0.6868, 0.7360, 0.7832, 0.8252, 0.8580, ...
              0.8882, 0.9148];

Mix4_small = [0.0322, 0.1254, 0.1576, 0.2006, 0.2540, 0.3148, 0.3768, 0.4408, ...
              0.5048, 0.5660, 0.6238, 0.6806, 0.7220, 0.7718, 0.8130, 0.8454, ...
              0.8762, 0.9038];

laplace_small = [0.0480, 0.1680, 0.2182, 0.2736, 0.3340, 0.4010, 0.4716, 0.5426, ...
                 0.6066, 0.6602, 0.7162, 0.7632, 0.8036, 0.8410, 0.8712, 0.8966, ...
                 0.9150, 0.9306];

% Moderate sample size group (N5)
N_Mod = [0.0542, 0.3238, 0.4220, 0.5440, 0.6498, 0.7396, 0.8178, 0.8814, ...
         0.9232, 0.9482, 0.9682, 0.9800, 0.9886, 0.9950, 0.9964, 0.9984, ...
         0.9988, 0.9994];

Mix1_Mod = [0.0380, 0.2818, 0.3868, 0.4944, 0.5986, 0.6922, 0.7742, 0.8380, ...
            0.8962, 0.9310, 0.9540, 0.9704, 0.9806, 0.9874, 0.9936, 0.9974, ...
            0.9988, 0.9996];

Mix2_Mod = [0.0382, 0.2784, 0.3816, 0.4828, 0.5894, 0.6888, 0.7654, 0.8308, ...
            0.8812, 0.9200, 0.9460, 0.9630, 0.9770, 0.9838, 0.9920, 0.9958, ...
            0.9976, 0.9990];

Mix3_Mod = [0.0406, 0.2892, 0.3870, 0.5060, 0.6130, 0.7082, 0.7842, 0.8508, ...
            0.8978, 0.9316, 0.9538, 0.9718, 0.9836, 0.9910, 0.9956, 0.9978, ...
            0.9988, 0.9994];

Mix4_Mod = [0.0406, 0.2884, 0.3858, 0.4974, 0.6068, 0.7072, 0.7788, 0.8410, ...
            0.8934, 0.9272, 0.9512, 0.9676, 0.9816, 0.9892, 0.9944, 0.9972, ...
            0.9986, 0.9992];

laplace_Mod = [0.0476, 0.3184, 0.4290, 0.5422, 0.6508, 0.7394, 0.8150, 0.8760, ...
               0.9120, 0.9426, 0.9656, 0.9808, 0.9896, 0.9944, 0.9972, 0.9982, ...
               0.9988, 1.0000];


N_small = [0.0460, N_small(2:end) ./ sqrt(1 + abs(1 - (0.0460 / 0.05)))];
Mix1_small = [0.0324, Mix1_small(2:end) ./ sqrt(1 + abs(1 - (0.0324 / 0.05)))];
Mix2_small = [0.0318, Mix2_small(2:end) ./ sqrt(1 + abs(1 - (0.0318 / 0.05)))];
Mix3_small = [0.0338, Mix3_small(2:end) ./ sqrt(1 + abs(1 - (0.0338 / 0.05)))];
Mix4_small = [0.0322, Mix4_small(2:end) ./ sqrt(1 + abs(1 - (0.0322 / 0.05)))];
laplace_small = [0.0480, laplace_small(2:end) ./ sqrt(1 + abs(1 - (0.0480 / 0.05)))];


N_Mod = [0.0542, N_Mod(2:end) ./ sqrt(1 + abs(1 - (0.0542 / 0.05)))];
Mix1_Mod = [0.0380, Mix1_Mod(2:end) ./ sqrt(1 + abs(1 - (0.0380 / 0.05)))];
Mix2_Mod = [0.0382, Mix2_Mod(2:end) ./ sqrt(1 + abs(1 - (0.0382 / 0.05)))];
Mix3_Mod = [0.0406, Mix3_Mod(2:end) ./ sqrt(1 + abs(1 - (0.0406 / 0.05)))];
Mix4_Mod = [0.0406, Mix4_Mod(2:end) ./ sqrt(1 + abs(1 - (0.0406 / 0.05)))];
laplace_Mod = [0.0476, laplace_Mod(2:end) ./ sqrt(1 + abs(1 - (0.0476 / 0.05)))];



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


%print('Robustness_MinMax', '-depsc', '-r1200'); % EPS for LaTeX/journal submission


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