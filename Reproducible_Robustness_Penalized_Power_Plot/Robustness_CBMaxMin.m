
N_small = [0.0504, 0.1614, 0.1998, 0.2370, 0.2844, 0.3256, 0.3702, 0.4142, 0.4600, ...
           0.4966, 0.5360, 0.5812, 0.6248, 0.6670, 0.7030, 0.7336, 0.7652, 0.7934];

Mix1_small = [0.0580, 0.1716, 0.2076, 0.2480, 0.2836, 0.3206, 0.3630, 0.4044, 0.4398, ...
              0.4766, 0.5134, 0.5498, 0.5854, 0.6182, 0.6460, 0.6764, 0.7036, 0.7302];

Mix2_small = [0.0632, 0.1780, 0.2146, 0.2512, 0.2902, 0.3308, 0.3722, 0.4126, 0.4512, ...
              0.4852, 0.5194, 0.5520, 0.5846, 0.6164, 0.6446, 0.6726, 0.6956, 0.7190];

Mix3_small = [0.0670, 0.1786, 0.2136, 0.2506, 0.2876, 0.3266, 0.3718, 0.4170, 0.4564, ...
              0.4970, 0.5330, 0.5706, 0.6032, 0.6368, 0.6674, 0.6972, 0.7254, 0.7506];

Mix4_small = [0.0690, 0.1818, 0.2198, 0.2560, 0.2996, 0.3376, 0.3816, 0.4202, 0.4652, ...
              0.4982, 0.5364, 0.5734, 0.6070, 0.6350, 0.6616, 0.6918, 0.7200, 0.7460];

laplace_small = [0.0486, 0.1658, 0.2092, 0.2510, 0.2992, 0.3466, 0.3936, 0.4400, 0.4878, ...
                 0.5288, 0.5754, 0.6162, 0.6536, 0.6886, 0.7226, 0.7532, 0.7818, 0.8068];

N_Mod = [0.0520, 0.2748, 0.3492, 0.4314, 0.5118, 0.5964, 0.6672, 0.7228, 0.7854, ...
         0.8352, 0.8738, 0.9082, 0.9296, 0.9486, 0.9626, 0.9742, 0.9822, 0.9886];

Mix1_Mod = [0.0544, 0.2876, 0.3714, 0.4558, 0.5270, 0.5958, 0.6642, 0.7240, 0.7736, ...
            0.8244, 0.8624, 0.8930, 0.9176, 0.9380, 0.9570, 0.9710, 0.9780, 0.9862];

Mix2_Mod = [0.0578, 0.2980, 0.3776, 0.4568, 0.5304, 0.5988, 0.6592, 0.7138, 0.7732, ...
            0.8186, 0.8546, 0.8818, 0.9100, 0.9312, 0.9498, 0.9622, 0.9754, 0.9828];

Mix3_Mod = [0.0584, 0.2958, 0.3738, 0.4516, 0.5278, 0.6036, 0.6610, 0.7226, 0.7752, ...
            0.8202, 0.8628, 0.8972, 0.9244, 0.9426, 0.9594, 0.9718, 0.9814, 0.9880];

Mix4_Mod = [0.0618, 0.2976, 0.3768, 0.4538, 0.5316, 0.6012, 0.6604, 0.7188, 0.7676, ...
            0.8200, 0.8576, 0.8904, 0.9190, 0.9394, 0.9580, 0.9702, 0.9784, 0.9848];

laplace_Mod = [0.0460, 0.2804, 0.3602, 0.4426, 0.5220, 0.6006, 0.6696, 0.7284, 0.7842, ...
               0.8282, 0.8664, 0.8942, 0.9222, 0.9402, 0.9568, 0.9714, 0.9782, 0.9838];



N_small = [0.0504, N_small(2:end) ./ sqrt(1 + abs(1 - (0.0504 / 0.05)))];
Mix1_small = [0.0580, Mix1_small(2:end) ./ sqrt(1 + abs(1 - (0.0580 / 0.05)))];
Mix2_small = [0.0632, Mix2_small(2:end) ./ sqrt(1 + abs(1 - (0.0632 / 0.05)))];
Mix3_small = [0.0670, Mix3_small(2:end) ./ sqrt(1 + abs(1 - (0.0670 / 0.05)))];
Mix4_small = [0.0690, Mix4_small(2:end) ./ sqrt(1 + abs(1 - (0.0690 / 0.05)))];
laplace_small = [0.0486, laplace_small(2:end) ./ sqrt(1 + abs(1 - (0.0486 / 0.05)))];


N_Mod = [0.0520, N_Mod(2:end) ./ sqrt(1 + abs(1 - (0.0520 / 0.05)))];
Mix1_Mod = [0.0544, Mix1_Mod(2:end) ./ sqrt(1 + abs(1 - (0.0544 / 0.05)))];
Mix2_Mod = [0.0578, Mix2_Mod(2:end) ./ sqrt(1 + abs(1 - (0.0578 / 0.05)))];
Mix3_Mod = [0.0584, Mix3_Mod(2:end) ./ sqrt(1 + abs(1 - (0.0584 / 0.05)))];
Mix4_Mod = [0.0618, Mix4_Mod(2:end) ./ sqrt(1 + abs(1 - (0.0618 / 0.05)))];
laplace_Mod = [0.0460, laplace_Mod(2:end) ./ sqrt(1 + abs(1 - (0.0460 / 0.05)))];


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


%print('Robustness_MaxMin', '-depsc', '-r1200'); % EPS for LaTeX/journal submission


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