
N_small = [0.0468, 0.1454, 0.1826, 0.2176, 0.2548, 0.2868, 0.3198, 0.3494, 0.3802, 0.4070, 0.4276, 0.4460, 0.4680, 0.4846, 0.5004, 0.5144, 0.5284, 0.5394];

Mix1_small = [0.0508, 0.1424, 0.1736, 0.2018, 0.2378, 0.2728, 0.3062, 0.3350, 0.3632, 0.3874, 0.4072, 0.4282, 0.4450, 0.4622, 0.4770, 0.4936, 0.5074, 0.5192];

Mix2_small = [0.0520, 0.1462, 0.1708, 0.2052, 0.2376, 0.2694, 0.3006, 0.3302, 0.3584, 0.3848, 0.4086, 0.4254, 0.4456, 0.4604, 0.4780, 0.4924, 0.5040, 0.5158];

Mix3_small = [0.0558, 0.1496, 0.1806, 0.2174, 0.2494, 0.2794, 0.3100, 0.3360, 0.3638, 0.3894, 0.4132, 0.4352, 0.4552, 0.4736, 0.4920, 0.5074, 0.5220, 0.5304];

Mix4_small = [0.0570, 0.1558, 0.1850, 0.2156, 0.2492, 0.2794, 0.3104, 0.3372, 0.3654, 0.3914, 0.4124, 0.4366, 0.4538, 0.4720, 0.4918, 0.5066, 0.5196, 0.5322];

laplace_small = [0.0520, 0.1634, 0.1978, 0.2338, 0.2686, 0.3034, 0.3350, 0.3630, 0.3952, 0.4206, 0.4430, 0.4650, 0.4822, 0.4990, 0.5174, 0.5326, 0.5474, 0.5590];

N_Mod = [0.0510, 0.2472, 0.2998, 0.3528, 0.4012, 0.4386, 0.4672, 0.4922, 0.5136, 0.5362, 0.5562, 0.5732, 0.5868, 0.6016, 0.6176, 0.6284, 0.6424, 0.6548];

Mix1_Mod = [0.0502, 0.2538, 0.3100, 0.3610, 0.4042, 0.4406, 0.4730, 0.4994, 0.5198, 0.5404, 0.5586, 0.5806, 0.5968, 0.6124, 0.6250, 0.6370, 0.6466, 0.6592];

Mix2_Mod = [0.0512, 0.2570, 0.3102, 0.3626, 0.4060, 0.4412, 0.4686, 0.4960, 0.5190, 0.5382, 0.5572, 0.5774, 0.5942, 0.6098, 0.6242, 0.6338, 0.6452, 0.6574];

Mix3_Mod = [0.0500, 0.2580, 0.3112, 0.3606, 0.4036, 0.4450, 0.4726, 0.4998, 0.5238, 0.5456, 0.5634, 0.5794, 0.5970, 0.6138, 0.6266, 0.6366, 0.6472, 0.6584];

Mix4_Mod = [0.0512, 0.2544, 0.3136, 0.3612, 0.4046, 0.4410, 0.4720, 0.4980, 0.5216, 0.5440, 0.5598, 0.5794, 0.5972, 0.6104, 0.6228, 0.6348, 0.6438, 0.6568];

laplace_Mod = [0.0508, 0.2590, 0.3144, 0.3668, 0.4052, 0.4458, 0.4792, 0.5030, 0.5276, 0.5470, 0.5648, 0.5830, 0.5988, 0.6128, 0.6292, 0.6434, 0.6536, 0.6664];



N_small = [0.0468, N_small(2:end) ./ sqrt(1 + abs(1 - (0.0468 / 0.05)))];
Mix1_small = [0.0508, Mix1_small(2:end) ./ sqrt(1 + abs(1 - (0.0508 / 0.05)))];
Mix2_small = [0.0520, Mix2_small(2:end) ./ sqrt(1 + abs(1 - (0.0520 / 0.05)))];
Mix3_small = [0.0558, Mix3_small(2:end) ./ sqrt(1 + abs(1 - (0.0558 / 0.05)))];
Mix4_small = [0.0570, Mix4_small(2:end) ./ sqrt(1 + abs(1 - (0.0570 / 0.05)))];
laplace_small = [0.0520, laplace_small(2:end) ./ sqrt(1 + abs(1 - (0.0520 / 0.05)))];


N_Mod = [0.0510, N_Mod(2:end) ./ sqrt(1 + abs(1 - (0.0510 / 0.05)))];
Mix1_Mod = [0.0502, Mix1_Mod(2:end) ./ sqrt(1 + abs(1 - (0.0502 / 0.05)))];
Mix2_Mod = [0.0512, Mix2_Mod(2:end) ./ sqrt(1 + abs(1 - (0.0512 / 0.05)))];
Mix3_Mod = [0.0500, Mix3_Mod(2:end) ./ sqrt(1 + abs(1 - (0.0500 / 0.05)))];
Mix4_Mod = [0.0512, Mix4_Mod(2:end) ./ sqrt(1 + abs(1 - (0.0512 / 0.05)))];
laplace_Mod = [0.0508, laplace_Mod(2:end) ./ sqrt(1 + abs(1 - (0.0508 / 0.05)))];


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


%print('Robustness_Min', '-depsc', '-r1200'); % EPS for LaTeX/journal submission


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