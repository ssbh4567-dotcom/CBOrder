
% Small sample size group
N_small = [0.0450 0.1444 0.1936 0.2420 0.3012 0.3702 0.4414 0.5114 0.5850 0.6570 0.7282 0.7900 0.8386 0.8856 0.9214 0.9490 0.9664 0.9782];
Mix1_small = [0.0352 0.1280 0.1674 0.2148 0.2686 0.3290 0.3928 0.4602 0.5268 0.5932 0.6534 0.7124 0.7676 0.8138 0.8564 0.8904 0.9160 0.9386];
Mix2_small = [0.0318 0.1300 0.1770 0.2262 0.2852 0.3462 0.4112 0.4816 0.5456 0.6080 0.6636 0.7202 0.7688 0.8108 0.8518 0.8814 0.9060 0.9282];
Mix3_small = [0.0398 0.1396 0.1844 0.2358 0.2936 0.3558 0.4156 0.4792 0.5488 0.6148 0.6796 0.7420 0.7980 0.8408 0.8774 0.9090 0.9346 0.9522];
Mix4_small = [0.0364 0.1426 0.1906 0.2460 0.3096 0.3726 0.4380 0.4994 0.5644 0.6340 0.6930 0.7512 0.8024 0.8432 0.8768 0.9048 0.9290 0.9484];
laplace_small = [0.0464 0.1700 0.2252 0.2866 0.3556 0.4382 0.5158 0.5944 0.6594 0.7244 0.7854 0.8370 0.8776 0.9098 0.9334 0.9522 0.9680 0.9766];

% Moderate sample size group
N_Mod = [0.0530 0.3430 0.4766 0.6204 0.7498 0.8544 0.9220 0.9668 0.9884 0.9954 0.9990 0.9994 0.9998 1.0000 1.0000 1.0000 1.0000 1.0000];
Mix1_Mod = [0.0410 0.3184 0.4596 0.5950 0.7230 0.8266 0.9062 0.9534 0.9768 0.9910 0.9974 0.9996 0.9998 1.0000 1.0000 1.0000 1.0000 1.0000];
Mix2_Mod = [0.0384 0.3276 0.4616 0.6038 0.7256 0.8278 0.9010 0.9494 0.9728 0.9876 0.9954 0.9990 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000];
Mix3_Mod = [0.0422 0.3302 0.4688 0.6014 0.7300 0.8360 0.9074 0.9538 0.9808 0.9946 0.9980 0.9994 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000];
Mix4_Mod = [0.0464 0.3396 0.4698 0.6052 0.7286 0.8348 0.9050 0.9524 0.9786 0.9926 0.9976 0.9992 0.9998 1.0000 1.0000 1.0000 1.0000 1.0000];
laplace_Mod = [0.0474 0.3632 0.5062 0.6468 0.7638 0.8528 0.9210 0.9662 0.9848 0.9950 0.9988 0.9992 0.9996 0.9998 1.0000 1.0000 1.0000 1.0000];

N_small = [0.0450, N_small(2:end) ./ sqrt(1 + abs(1 - (0.0450 / 0.05)))];
Mix1_small = [0.0352, Mix1_small(2:end) ./ sqrt(1 + abs(1 - (0.0352 / 0.05)))];
Mix2_small = [0.0318, Mix2_small(2:end) ./ sqrt(1 + abs(1 - (0.0318 / 0.05)))];
Mix3_small = [0.0398, Mix3_small(2:end) ./ sqrt(1 + abs(1 - (0.0398 / 0.05)))];
Mix4_small = [0.0364, Mix4_small(2:end) ./ sqrt(1 + abs(1 - (0.0364 / 0.05)))];
laplace_small = [0.0464, laplace_small(2:end) ./ sqrt(1 + abs(1 - (0.0464 / 0.05)))];


N_Mod = [0.0530, N_Mod(2:end) ./ sqrt(1 + abs(1 - (0.0530 / 0.05)))];
Mix1_Mod = [0.0410, Mix1_Mod(2:end) ./ sqrt(1 + abs(1 - (0.0410 / 0.05)))];
Mix2_Mod = [0.0384, Mix2_Mod(2:end) ./ sqrt(1 + abs(1 - (0.0384 / 0.05)))];
Mix3_Mod = [0.0422, Mix3_Mod(2:end) ./ sqrt(1 + abs(1 - (0.0422 / 0.05)))];
Mix4_Mod = [0.0464, Mix4_Mod(2:end) ./ sqrt(1 + abs(1 - (0.0464 / 0.05)))];
laplace_Mod = [0.0474, laplace_Mod(2:end) ./ sqrt(1 + abs(1 - (0.0474 / 0.05)))];



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


%print('Robustness_MAX', '-depsc', '-r1200'); % EPS for LaTeX/journal submission


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
