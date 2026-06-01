%% =====================================================================
%  三维颗粒气体密度均匀性 KS 检验 (Kolmogorov-Smirnov Test)
%  只检验 Bulk (体区) 的颗粒空间分布是否服从理想的均匀分布
% =====================================================================
clear all; close all; clc;

%% ========== 1. 用户配置区 ==========
% 数据文件路径 (换成你其中一个代表性的数据，比如 phi=0.01 的文件)
data_file = 'G:\Space_Active\Data\Maxwell_Boltz\data\Pos_Vel_data\Data_FixSpace_R0_04_mode0td_phi0.05_Tag2_Step1_to_3000_20260402_141317.mat'; 
tag_name  = 'Tag_2';

L_box = 4.0;          % 容器边长 [cm]
bulk_ratio = 0.6;     % 取中心 60% 区域作为体区 (可改为 0.5 或 0.7)
half_L_bulk = (L_box * bulk_ratio) / 2;

delta_step = 10;      % 抽样步长 (不需要每帧都算，抽样即可)

fprintf('=== 开始进行空间密度均匀性 KS 检验 ===\n');
fprintf('  检测区域: Bulk (边长 = %.2f cm)\n', L_box * bulk_ratio);

%% ========== 2. 数据加载与 Bulk 过滤 ==========
loaded = load(data_file);
PD = loaded.ExtractedData.(tag_name);
valid_frames = 1 : delta_step : length(PD.Pos);

pos_bulk_all = [];

for j = valid_frames
    pos_j = PD.Pos{j};
    if isempty(pos_j), continue; end
    
    % 清洗 NaN
    nan_mask = any(isnan(pos_j), 2);
    pos_j(nan_mask, :) = [];
    
    % 过滤出位于 Bulk 内部的颗粒 (假设中心在 0,0,0)
    dist_from_center = max(abs(pos_j), [], 2);
    bulk_mask = dist_from_center < half_L_bulk;
    pos_bulk_all = [pos_bulk_all; pos_j(bulk_mask, :)];
end

N_samples = size(pos_bulk_all, 1);
fprintf('  成功提取 Bulk 内有效颗粒坐标点数: %d\n', N_samples);

% ========== 3. KS 检验计算 ==========
% 定义 Bulk 的理论均匀分布边界
lower_bound = -half_L_bulk;
upper_bound = half_L_bulk;

% 构建理论均匀分布对象 (Uniform Distribution)
pd_unif = makedist('Uniform', 'lower', lower_bound, 'upper', upper_bound);

% 分别对 X, Y, Z 三个方向做 One-sample KS test
% 注意：当 N 极大时，p-value 必然趋近于0。
% 因此在物理上，我们主要看最大统计距离 KS Statistic (D 值)，D 越小说明越均匀
[h_x, p_x, ksstat_x] = kstest(pos_bulk_all(:,1), 'CDF', pd_unif);
[h_y, p_y, ksstat_y] = kstest(pos_bulk_all(:,2), 'CDF', pd_unif);
[h_z, p_z, ksstat_z] = kstest(pos_bulk_all(:,3), 'CDF', pd_unif);

fprintf('\n=== KS 检验结果 (统计距离 D) ===\n');
fprintf('  (D 值越接近 0 代表越完美符合均匀分布)\n');
fprintf('  X 方向最大偏差 D_x = %.4f\n', ksstat_x);
fprintf('  Y 方向最大偏差 D_y = %.4f\n', ksstat_y);
fprintf('  Z 方向最大偏差 D_z = %.4f\n', ksstat_z);

% ========== 4. 绘制 CDF 对比图 (极具说服力的文章图) ==========
figure('Name', 'Density Uniformity CDF', 'Position', [150, 150, 600, 500]);
hold on;

% 1. 画理论均匀分布的 CDF (一条对角直线)
x_theory = linspace(lower_bound, upper_bound, 100);
cdf_theory = cdf(pd_unif, x_theory);
plot(x_theory, cdf_theory, 'k--', 'LineWidth', 2.5, 'DisplayName', 'Ideal Uniform');

% 2. 画 X, Y, Z 三个方向的经验 CDF (Empirical CDF)
colors = {'r', 'b', 'g'};
labels = {'X-axis', 'Y-axis', 'Z-axis'};

for dim = 1:3
    [f_emp, x_emp] = ecdf(pos_bulk_all(:, dim));
    
    % 为了防止画图卡顿，如果点太多，降采样画图
    if length(x_emp) > 5000
        idx = round(linspace(1, length(x_emp), 5000));
        x_emp = x_emp(idx);
        f_emp = f_emp(idx);
    end
    
    plot(x_emp, f_emp, '-', 'Color', colors{dim}, 'LineWidth', 1.5, ...
        'DisplayName', sprintf('%s (D=%.3f)', labels{dim}, eval(sprintf('ksstat_%c', lower(labels{dim}(1))))));
end

xlabel('Position in Bulk [cm]', 'FontSize', 16);
ylabel('Cumulative Probability (CDF)', 'FontSize', 16);
% title('Spatial Uniformity Check (KS Test)', 'FontSize', 18);
legend('Location', 'northwest', 'FontSize', 18); legend boxoff;
set(gca, 'FontSize', 24, 'FontName', 'Times New Roman', 'LineWidth', 1.5);


xlim([lower_bound, upper_bound]);
ylim([0, 1]);

%%
% =====================================================================
%  三维颗粒气体密度均匀性 KS 检验 (批量处理演化版)
%  自动遍历多个 \phi 的数据，计算体区 KS 距离 D，并绘制 D_avg vs \phi
% =====================================================================
clear all; close all; clc;

%% ========== 1. 用户配置区 ==========
% 数据文件夹路径与匹配模式
pos_vel_dir = 'G:\Space_Active\Data\Maxwell_Boltz\data\Pos_Vel_data';
pv_pattern  = 'Data_FixSpace_R0_04_mode0td_phi*.mat';
tag_name    = 'Tag_2';

L_box = 4.0;          % 容器边长 [cm]
bulk_ratio = 0.6;     % 取中心 60% 区域作为体区 (可调)
half_L_bulk = (L_box * bulk_ratio) / 2;

delta_step = 10;      % 抽样步长

fprintf('=== 启动空间密度均匀性批量 KS 检验 ===\n');
fprintf('  检测区域: Bulk (边长 = %.2f cm)\n', L_box * bulk_ratio);

%% ========== 2. 自动检索与提取 \phi 值 ==========
pv_files = dir(fullfile(pos_vel_dir, pv_pattern));
if isempty(pv_files)
    error('未找到匹配的文件，请检查路径和通配符！');
end

% 提取文件名中的 phi
extract_phi = @(fname) str2double(regexp(fname, '(?<=phi)\d+\.?\d*', 'match', 'once'));
pv_map = containers.Map('KeyType', 'double', 'ValueType', 'char');
for k = 1:length(pv_files)
    phi_k = extract_phi(pv_files(k).name);
    if ~isnan(phi_k)
        pv_map(phi_k) = pv_files(k).name;
    end
end

phi_vals_all = sort(cell2mat(pv_map.keys));
n_sets = length(phi_vals_all);

%% ========== 3. 批量计算 KS 统计距离 ==========
phi_vals = zeros(n_sets, 1);
D_x_vals = zeros(n_sets, 1);
D_y_vals = zeros(n_sets, 1);
D_z_vals = zeros(n_sets, 1);
D_avg_vals = zeros(n_sets, 1);

% 构建体区 (Bulk) 理论均匀分布对象
pd_unif = makedist('Uniform', 'lower', -half_L_bulk, 'upper', half_L_bulk);

for i = 1:n_sets
    phi_val = phi_vals_all(i);
    pv_file = fullfile(pos_vel_dir, pv_map(phi_val));
    phi_vals(i) = phi_val;

    fprintf('\n[%d/%d] 正在处理 φ = %.4f ...\n', i, n_sets, phi_val);

    % 加载数据
    loaded = load(pv_file);
    PD = loaded.ExtractedData.(tag_name);
    valid_frames = 1 : delta_step : length(PD.Pos);

    pos_bulk_all = [];
    for j = valid_frames
        pos_j = PD.Pos{j};
        if isempty(pos_j), continue; end
        
        % 清洗 NaN
        nan_mask = any(isnan(pos_j), 2);
        pos_j(nan_mask, :) = [];

        % 空间过滤：只保留位于 Bulk 内部的颗粒
        dist_from_center = max(abs(pos_j), [], 2);
        bulk_mask = dist_from_center < half_L_bulk;
        pos_bulk_all = [pos_bulk_all; pos_j(bulk_mask, :)];
    end

    if isempty(pos_bulk_all)
        warning('提取到的 Bulk 颗粒数为空！跳过。');
        continue;
    end

    % --- 核心计算: One-sample KS test ---
    % 返回的 ksstat 就是最大经验分布偏差 D 值
    [~, ~, ksstat_x] = kstest(pos_bulk_all(:,1), 'CDF', pd_unif);
    [~, ~, ksstat_y] = kstest(pos_bulk_all(:,2), 'CDF', pd_unif);
    [~, ~, ksstat_z] = kstest(pos_bulk_all(:,3), 'CDF', pd_unif);

    D_x_vals(i) = ksstat_x;
    D_y_vals(i) = ksstat_y;
    D_z_vals(i) = ksstat_z;
    
    % 计算平均 D 值
    D_avg_vals(i) = (ksstat_x + ksstat_y + ksstat_z) / 3;

    fprintf('  有效点数: %d | D_x=%.4f, D_y=%.4f, D_z=%.4f -> D_avg=%.4f\n', ...
        size(pos_bulk_all,1), ksstat_x, ksstat_y, ksstat_z, D_avg_vals(i));
end

%% ========== 4. 绘制 D_avg vs phi 演化图 ==========
figure('Name', 'KS Statistic vs Volume Fraction', 'Position', [200, 200, 750, 550]);

% 画出各个方向的参考线 (细线，空心点，增加透明感)
semilogx(phi_vals, D_x_vals, 'ro-', 'LineWidth', 1, 'MarkerSize', 5, 'MarkerFaceColor', 'w', 'DisplayName', 'D_x'); hold on;
semilogx(phi_vals, D_y_vals, 'b^-', 'LineWidth', 1, 'MarkerSize', 5, 'MarkerFaceColor', 'w', 'DisplayName', 'D_y');
semilogx(phi_vals, D_z_vals, 'gs-', 'LineWidth', 1, 'MarkerSize', 5, 'MarkerFaceColor', 'w', 'DisplayName', 'D_z');

% 【核心主线】：画出平均 D 值的演化 (黑色粗实线，实心大点)
semilogx(phi_vals, D_avg_vals, 'k-o', 'LineWidth', 3, 'MarkerSize', 10, ...
    'MarkerFaceColor', [0.3 0.3 0.3], 'DisplayName', 'Average D');

% 画一条经验参考线 D=0.05 (在统计物理中，大样本下 D < 0.05 意味着极其优异的均匀性)
yline(0.05, 'k:', 'D = 0.05 (Excellent Uniformity Benchmark)', 'LineWidth', 2, ...
    'LabelHorizontalAlignment', 'left', 'FontSize', 14, 'HandleVisibility', 'off');

xlabel('\phi', 'FontSize', 18);
ylabel('KS Statistic Distance D', 'FontSize', 18);
% title('Spatial Uniformity Evolution in Bulk Area', 'FontSize', 18);
legend('Location', 'northwest', 'FontSize', 18); legend boxoff;
set(gca, 'FontSize', 20, 'FontName', 'Times New Roman', 'LineWidth', 2);
% grid on;

% 优化 Y 轴显示范围 (从 0 开始，上限根据数据自适应)
ylim_curr = ylim;
ylim([0, max(ylim_curr(2)*1.2, 0.1)]);

%% ========== 5. 打印学术汇总表 ==========
fprintf('\n========== 均匀性检测 (KS Test) 汇总 ==========\n');
fprintf('%10s %10s %10s %10s %15s\n', 'phi', 'D_x', 'D_y', 'D_z', 'D_avg');
fprintf('%s\n', repmat('-', 1, 60));
for i = 1:n_sets
    fprintf('%10.4f %10.4f %10.4f %10.4f %15.4f\n', ...
        phi_vals(i), D_x_vals(i), D_y_vals(i), D_z_vals(i), D_avg_vals(i));
end
fprintf('============================================================\n');