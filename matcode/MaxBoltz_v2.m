%% =====================================================================
%  多体积分数批量处理与速率分布对比脚本 (终极自动化版)
%  功能: 通配符自动扫盘 -> 正则提取 phi -> 自动排序 -> 拟合 VDF -> 绘制 Beta 演化图
% =====================================================================

clear all; close all; clc;

%% ========== 1. 用户配置区 (自动化版) ==========

% --- 1.1 数据文件夹路径与通配符模式 ---
% 请确保路径正确，如果不在此文件夹，请修改
data_dir = '/media/gezhuan/M78/Space_Active/Data/Active_Gas/Data/Pos_Vel';
% 你的文件名通配符模式
file_pattern = 'Data*.mat'; 

% --- 1.2 基础系统参数 ---
data_in_cm = true;           % 原始数据是否为 cm 和 cm/s
tag_name = 'Tag_1';          % 颗粒字段名
delta_step = 1;              % 帧步长
n_bins_vdf = 60;             % 直方图 bin 数量

if data_in_cm
    scale_factor = 10.0;     % 转换为 mm
else
    scale_factor = 1.0;
end

speed_edges = linspace(0, 5, n_bins_vdf + 1);
speed_bin_centers = (speed_edges(1:end-1) + speed_edges(2:end)) / 2;

%% ========== 2. 自动检索与正则表达式提取 ==========

fprintf('=== 开始自动检索数据文件 ===\n');
fprintf('  路径: %s\n', data_dir);

file_list = dir(fullfile(data_dir, file_pattern));
if isempty(file_list)
    error('未找到匹配 "%s" 的文件！请检查路径或通配符。', file_pattern);
end

% 使用正则表达式从文件名中提取 phi 值
extract_phi = @(fname) str2double(regexp(fname, '(?<=phi)\d+\.?\d*', 'match', 'once'));

% 使用 Map 来存储 (phi -> filename)，自动去重(如果有同一phi的多个文件，保留最新的)
file_map = containers.Map('KeyType', 'double', 'ValueType', 'char');
for k = 1:length(file_list)
    phi_val = extract_phi(file_list(k).name);
    if ~isnan(phi_val)
        file_map(phi_val) = file_list(k).name;
    end
end

% 提取所有找到的 phi 值，并从大到小排序 (画 VDF 图时从大到小比较好看)
phi_vals_numeric = sort(cell2mat(file_map.keys), 'descend');
n_files = length(phi_vals_numeric);

fprintf('  成功匹配到 %d 组有效数据:\n', n_files);
for k = 1:n_files
    fprintf('    - phi = %-8.4g : %s\n', phi_vals_numeric(k), file_map(phi_vals_numeric(k)));
end

% 准备配色和数据存储容器
color_map = lines(n_files); 
beta_vals_numeric = zeros(n_files, 1);
phi_labels = cell(n_files, 1);

%% ========== 3. 准备画图容器 ==========

figure('Name', 'Speed Distribution Variation with Volume Fraction', ...
    'Position', [150, 150, 950, 750]);
hold on;

fprintf('\n=== 开始批量处理与拟合 ===\n');

%% ========== 4. 循环批量处理与画图 ==========

parfor i = 1:n_files
    current_phi = phi_vals_numeric(i);
    current_file = file_map(current_phi);
    phi_labels{i} = sprintf('\\phi = %g', current_phi); % 自动生成图例标签
    
    fprintf('\n[%d/%d] 正在处理: φ = %g ...\n', i, n_files, current_phi);
    
    try
        loaded = load(fullfile(data_dir, current_file));
        PD = loaded.ExtractedData.(tag_name);
    catch
        warning('文件加载失败或格式不对，跳过。');
        beta_vals_numeric(i) = NaN; 
        continue;
    end
    
    total_frames = length(PD.Pos);
    valid_frames = 1 : delta_step : total_frames;
    num_valid = length(valid_frames);
    
    vel_cell = cell(num_valid, 1);
    for j = 1:num_valid
        fi = valid_frames(j);
        vel_i = PD.Vel{fi};
        if ~isempty(vel_i)
            vel_cell{j} = vel_i * scale_factor;
        end
    end
    all_vel = cell2mat(vel_cell);
    
    nan_mask = any(isnan(all_vel), 2);
    if sum(nan_mask) > 0
        all_vel = all_vel(~nan_mask, :);
    end
    
    speed = sqrt(sum(all_vel.^2, 2));
    avg_speed_global = mean(speed);
    speed_normalized = speed / avg_speed_global;
    
    fprintf('      平均速率 <|v|> = %.4f mm/s\n', avg_speed_global);
    
    counts_speed = histcounts(speed_normalized, speed_edges, 'Normalization', 'pdf');
    valid_bins = counts_speed > 0;
    c_data = speed_bin_centers(valid_bins);
    f_data = counts_speed(valid_bins);
    
    semilogy(speed_bin_centers, counts_speed, 'o', ...
        'MarkerSize', 8, 'MarkerFaceColor', 'none', ...            
        'MarkerEdgeColor', color_map(i,:), 'LineWidth', 1.5, ...   
        'DisplayName', sprintf('%s Data', phi_labels{i}));
        
    mb_func = @(a, c) 4*pi * (a/pi)^(1.5) .* c.^2 .* exp(-a * c.^2);
    opts = optimset('Display', 'off');
    
    try
        a_fit = lsqcurvefit(mb_func, 0.8, c_data, f_data, 0.1, 5, opts);
    catch
        a_fit = 3 / (2 * mean(speed_normalized.^2));
    end
    
    c_fit = linspace(0, 5, 200);
    mb_fit_curve = mb_func(a_fit, c_fit);
    semilogy(c_fit, mb_fit_curve, '-', 'Color', color_map(i,:), 'LineWidth', 1.5, ...
        'HandleVisibility', 'off'); 
    
    tail_mask = c_data > 2;
    if sum(tail_mask) > 5
        c_tail = c_data(tail_mask);
        log_f_tail = log(f_data(tail_mask));
        
        p_tail_2 = polyfit(c_tail.^2, log_f_tail, 1);
        alpha_2 = -p_tail_2(1);
        
        try
            tail_func = @(params, c) params(1) * exp(-params(2) * c.^params(3));
            params0 = [max(f_data(tail_mask))*10, alpha_2, 2];
            params_fit = lsqcurvefit(tail_func, params0, c_tail, f_data(tail_mask), ...
                [0,0,0.5], [100,100,5], opts);
            
            beta_fit = params_fit(3);
            beta_vals_numeric(i) = beta_fit; % 存入数组
            
            c_tail_plot = linspace(2, 5, 100);
            tail_fit_curve = tail_func(params_fit, c_tail_plot);
            
            semilogy(c_tail_plot, tail_fit_curve, '--', 'Color', color_map(i,:), 'LineWidth', 2.5, ...
                'DisplayName', sprintf('%s Tail: \\beta=%.2f', phi_labels{i}, beta_fit));
            fprintf('      尾巴拟合参数: β = %.3f\n', beta_fit);
        catch
            fprintf('      尾巴拟合失败。\n');
            beta_vals_numeric(i) = NaN;
        end
    else
        beta_vals_numeric(i) = NaN;
    end
end

%% ========== 5. VDF 图表格式美化 ==========

xlabel('v / \langle|v|\rangle', 'FontSize', 18);
ylabel('Probability Density', 'FontSize', 18);

legend('Location', 'northeast', 'FontSize', 16);
legend boxoff;

set(gca,'FontSize',20,'FontName','Times New Roman','LineWidth',2);
set(gca, 'YScale', 'log'); 
ylim([1e-4, 1e0]);  
xlim([0, 5]);

%% ========== 6. 绘制 Beta 随体积分数的变化图 ==========

fprintf('\n=== 正在生成 Beta vs \\phi 演化图 ===\n');

% 画趋势图时，横坐标最好是从小到大，所以这里做个升序排序
[phi_sorted, sort_idx] = sort(phi_vals_numeric, 'ascend');
beta_sorted = beta_vals_numeric(sort_idx);

figure('Name', 'Tail Index Beta vs Volume Fraction', 'Position', [250, 250, 650, 500]);

% 画出实际测量的 Beta 数据线
plot(phi_sorted, beta_sorted, 'ro-', 'LineWidth', 2.5, 'MarkerSize', 10, ...
    'MarkerFaceColor', 'w', 'DisplayName', 'Simulated \beta');
hold on;

% 【物理高光】：画出理论参考线
yline(2.0, 'k--', 'Classical MB (\beta=2.0)', 'LineWidth', 2, ...
    'LabelHorizontalAlignment', 'right', 'FontSize', 16, 'HandleVisibility', 'off');
yline(1.5, 'b:', 'Homogeneous Theory (\beta=1.5)', 'LineWidth', 2, ...
    'LabelHorizontalAlignment', 'right', 'FontSize', 16, 'HandleVisibility', 'off');

xlabel('\phi', 'FontSize', 18);
ylabel('\beta', 'FontSize', 18);
% title('Evolution of High-Energy Tail Index', 'FontSize', 18);

set(gca, 'FontSize', 20, 'FontName', 'Times New Roman', 'LineWidth', 2);
% grid on;

% 自适应 Y 轴边界
valid_betas = beta_sorted(~isnan(beta_sorted));
if ~isempty(valid_betas)
    ylim([min(min(valid_betas)*0.8, 1.2), max(max(valid_betas)*1.1, 3.5)]);
end

fprintf('\n=== 🚀 所有自动化处理完毕！ ===\n');