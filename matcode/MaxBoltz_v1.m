%% =====================================================================
%  多体积分数批量处理与速率分布对比脚本 (对应文献 Fig. 7c)
%  功能: 自动遍历多个原始数据文件 -> 提取速度 -> 清洗数据 -> 拟合 -> 同图对比
% =====================================================================

clear all; close all; clc;

%% ========== 1. 用户配置区 ==========

% --- 1.1 填入你所有需要对比的原始 .mat 数据文件名 ---
% 假设这些文件都在你当前运行脚本的同一目录下
% raw_files = {
%     'Data_FixSpace_R0_04_mode0td_phi0.05_Tag2_Step1_to_3000_20260402_141317.mat', ... % %
%     'Data_FixSpace_R0_04_mode0td_phi0.03_Tag2_Step1_to_3000_20260402_140424.mat', ...  
%     'Data_FixSpace_R0_04_mode0td_phi0.01_Tag2_Step1_to_3000_20260402_134707.mat', ...
%     'Data_FixSpace_R0_04_mode0td_phi0.004_Tag2_Step1_to_3000_20260402_134222.mat'
%     };
% --- 1.2 图例标签 (与上方文件顺序一一对应) ---
% phi_labels = {'\phi = 0.05','\phi = 0.03', ...
%     '\phi = 0.01','\phi = 0.004'
%     };
raw_files = {

    'Data_FixSpace_R0_04_mode0td_phi0.02_Tag2_Step1_to_3000_20260402_140145.mat', ... %
    'Data_FixSpace_R0_04_mode0td_phi0.015_Tag2_Step1_to_3000_20260402_135256.mat', ...
    'Data_FixSpace_R0_04_mode0td_phi0.01_Tag2_Step1_to_3000_20260402_134707.mat' , ...% 替换
    'Data_FixSpace_R0_04_mode0td_phi0.004_Tag2_Step1_to_3000_20260402_134222.mat'
    };
% --- 1.2 图例标签 (与上方文件顺序一一对应) ---
phi_labels = {'\phi = 0.02','\phi = 0.015', ...
    '\phi = 0.01','\phi = 0.004'
    };

% --- 1.3 基础系统参数 ---
data_in_cm = true;           % 原始数据是否为 cm 和 cm/s
tag_name = 'Tag_2';          % 颗粒字段名
delta_step = 1;              % 帧步长 (设为 5 提速，设为 1 用全量数据)
n_bins_vdf = 60;             % 直方图 bin 数量

% 设置缩放因子 (cm -> mm)
if data_in_cm
    scale_factor = 10.0;
else
    scale_factor = 1.0;
end

% 预设统一的速率网格边界 (用于画图对齐)
speed_edges = linspace(0, 5, n_bins_vdf + 1);
speed_bin_centers = (speed_edges(1:end-1) + speed_edges(2:end)) / 2;

% 配色方案
color_map = lines(length(raw_files)); 

%% ========== 2. 准备画图容器 ==========

figure('Name', 'Speed Distribution Variation with Volume Fraction', ...
    'Position', [150, 150, 950, 750]);
hold on;

fprintf('=== 开始批量处理多组数据 ===\n');

% ========== 3. 循环批量处理与画图 ==========

for i = 1:length(raw_files)
    current_file = raw_files{i};
    fprintf('\n[%d/%d] 正在处理: %s ...\n', i, length(raw_files), current_file);
    
    % --- 3.1 加载原始数据 ---
    try
        loaded = load(current_file);
        PD = loaded.ExtractedData.(tag_name);
    catch
        warning('文件加载失败或格式不对: %s，跳过。', current_file);
        continue;
    end
    
    total_frames = length(PD.Pos);
    valid_frames = 1 : delta_step : total_frames;
    num_valid = length(valid_frames);
    
    % --- 3.2 极速提取与清洗速度数据 ---
    vel_cell = cell(num_valid, 1);
    for j = 1:num_valid
        fi = valid_frames(j);
        vel_i = PD.Vel{fi};
        if ~isempty(vel_i)
            vel_cell{j} = vel_i * scale_factor; % 直接转成 mm/s
        end
    end
    
    all_vel = cell2mat(vel_cell);
    
    % 数据清洗: 剔除包含 NaN 的行
    nan_mask = any(isnan(all_vel), 2);
    if sum(nan_mask) > 0
        all_vel = all_vel(~nan_mask, :);
    end
    
    % 计算速率 |v| 并归一化 v / <|v|>
    speed = sqrt(sum(all_vel.^2, 2));
    avg_speed_global = mean(speed);
    speed_normalized = speed / avg_speed_global;
    
    fprintf('      平均速率 <|v|> = %.4f mm/s\n', avg_speed_global);
    
    % --- 3.3 计算概率分布 (PDF) ---
    counts_speed = histcounts(speed_normalized, speed_edges, 'Normalization', 'pdf');
    valid_bins = counts_speed > 0;
    c_data = speed_bin_centers(valid_bins);
    f_data = counts_speed(valid_bins);
    
    % % 画散点 (Measurement)
    % semilogy(speed_bin_centers, counts_speed, 'o', ...
    %     'MarkerSize', 8, 'MarkerFaceColor', color_map(i,:), ...
    %     'MarkerEdgeColor', 'none', 'DisplayName', sprintf('%s Data', phi_labels{i}));
    semilogy(speed_bin_centers, counts_speed, 'o', ...
        'MarkerSize', 8, 'MarkerFaceColor', 'none', ...            % 背景掏空
        'MarkerEdgeColor', color_map(i,:), 'LineWidth', 1.5, ...   % 颜色赋给边框并加粗
        'DisplayName', sprintf('%s Data', phi_labels{i}));
    % --- 3.4 麦克斯韦-玻尔兹曼 (MB) 拟合 ---
    mb_func = @(a, c) 4*pi * (a/pi)^(1.5) .* c.^2 .* exp(-a * c.^2);
    opts = optimset('Display', 'off');
    
    try
        a_fit = lsqcurvefit(mb_func, 0.8, c_data, f_data, 0.1, 5, opts);
    catch
        a_fit = 3 / (2 * mean(speed_normalized.^2));
    end
    
    c_fit = linspace(0, 5, 200);
    mb_fit_curve = mb_func(a_fit, c_fit);
    
    % 画 MB 拟合实线
    semilogy(c_fit, mb_fit_curve, '-', 'Color', color_map(i,:), 'LineWidth', 1.5, ...
        'HandleVisibility', 'off'); % 让图例不显示这根实线，保持整洁
    
    % --- 3.5 高能胖尾拟合 (~exp(-α * c^β)) ---
    tail_mask = c_data > 2;
    if sum(tail_mask) > 5
        c_tail = c_data(tail_mask);
        log_f_tail = log(f_data(tail_mask));
        
        % 估计初始参数
        p_tail_2 = polyfit(c_tail.^2, log_f_tail, 1);
        alpha_2 = -p_tail_2(1);
        
        try
            tail_func = @(params, c) params(1) * exp(-params(2) * c.^params(3));
            params0 = [max(f_data(tail_mask))*10, alpha_2, 2];
            params_fit = lsqcurvefit(tail_func, params0, c_tail, f_data(tail_mask), ...
                [0,0,0.5], [100,100,5], opts);
            beta_fit = params_fit(3);
            
            c_tail_plot = linspace(2, 5, 100);
            tail_fit_curve = tail_func(params_fit, c_tail_plot);
            
            % 画胖尾拟合虚线
            semilogy(c_tail_plot, tail_fit_curve, '--', 'Color', color_map(i,:), 'LineWidth', 2.5, ...
                'DisplayName', sprintf('%s Tail: \\beta=%.2f', phi_labels{i}, beta_fit));
            fprintf('      尾巴拟合参数: β = %.3f\n', beta_fit);
        catch
            fprintf('      尾巴拟合失败。\n');
        end
    end
end

% ========== 4. 图表格式美化 ==========

xlabel('v / \langle|v|\rangle', 'FontSize', 18);
ylabel('Probability Density', 'FontSize', 18);
% title('Speed Distribution Variation with Volume Fraction (\phi)', 'FontSize', 20);

% 整理图例
legend('Location', 'northeast', 'FontSize', 18);
legend boxoff;

set(gca,'FontSize',20,'FontName','Times New Roman','LineWidth',2);
set(gca, 'YScale', 'log'); 

% grid on;
ylim([1e-4, 1e0]);  % 根据需要调整下限，1e-4 能更清晰地展现尾巴的区别
xlim([0, 5]);

fprintf('\n=== 所有数据处理完毕！对比图已生成 ===\n');