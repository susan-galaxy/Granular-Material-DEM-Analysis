%% =====================================================================
%  Knudsen 数 & 平动/转动温度比 (全局 vs 体区 对比版)
%  输入: Pos_Vel 数据 (.mat) + Rotation 数据 (.mat)
%  输出: 双对比趋势图 & 学术汇总表
% =====================================================================
clear all; close all; clc;

%% ========== 1. 用户配置区 ==========

% --- 1.1 数据路径 & 通配符模式 ---
pos_vel_dir = 'G:\Space_Active\Data\Maxwell_Boltz\data\Pos_Vel_data';
rot_dir     = 'G:\Space_Active\Data\Maxwell_Boltz\data\Rotation_data\Rotation_data';

pv_pattern  = 'Data_FixSpace_R0_04_mode0td_phi*_Tag2_*.mat';
rot_pattern = 'RotData_FixSpace_R0_04_mode0td_phi*_Tag2_*.mat';

% --- 1.2 系统参数 ---
tag_name     = 'Tag_2';        
d_particle   = 0.08;           % 颗粒直径 [cm]
R_particle   = d_particle / 2; % 颗粒半径 [cm]
L_box        = 4.0;            % 容器边长 [cm]
V_box        = L_box^3;        % 容器体积 [cm³]

% --- 1.3 采样与 Bulk 设置 ---
delta_step = 5;  
bulk_ratio = 0.5;  % 【关键设置】：体区占全边长的比例 (0.7 代表剔除外围 15% 的边界层)
half_L_bulk = (L_box * bulk_ratio) / 2; % 中心区的半宽

fprintf('=== 启动物理特性计算 (Bulk Ratio = %.2f) ===\n', bulk_ratio);

%% ========== 2. 自动检索与配对数据文件 (逻辑保持不变) ==========

pv_files  = dir(fullfile(pos_vel_dir, pv_pattern));
rot_files = dir(fullfile(rot_dir, rot_pattern));

if isempty(pv_files) || isempty(rot_files)
    error('未找到匹配的文件，请检查路径和通配符！');
end

extract_phi = @(fname) str2double(regexp(fname, '(?<=phi)\d+\.?\d*', 'match', 'once'));

pv_map = containers.Map('KeyType', 'double', 'ValueType', 'char');
for k = 1:length(pv_files)
    phi_k = extract_phi(pv_files(k).name);
    if ~isnan(phi_k), pv_map(phi_k) = pv_files(k).name; end
end

rot_map = containers.Map('KeyType', 'double', 'ValueType', 'char');
for k = 1:length(rot_files)
    phi_k = extract_phi(rot_files(k).name);
    if ~isnan(phi_k), rot_map(phi_k) = rot_files(k).name; end
end

phi_common = sort(intersect(cell2mat(pv_map.keys), cell2mat(rot_map.keys)));
n_sets = length(phi_common);

datasets = cell(n_sets, 4);
for k = 1:n_sets
    datasets{k, 1} = pv_map(phi_common(k));
    datasets{k, 2} = rot_map(phi_common(k));
    datasets{k, 3} = phi_common(k);
end

%% ========== 3. 批量计算 (全局 vs 体区) ==========

phi_vals   = zeros(n_sets, 1);
Kn_vals    = zeros(n_sets, 1);

% 存储全局 (Global) 结果
Tratio_glb = zeros(n_sets, 1);
% 存储体区 (Bulk) 结果
Tratio_blk = zeros(n_sets, 1);

for i = 1:n_sets
    pv_file  = fullfile(pos_vel_dir, datasets{i, 1});
    rot_file = fullfile(rot_dir,     datasets{i, 2});
    phi_val  = datasets{i, 3};
    phi_vals(i) = phi_val;
    
    fprintf('\n[%d/%d] φ = %.4f 正在处理...\n', i, n_sets, phi_val);
    
    pv_data  = load(pv_file);
    rot_data = load(rot_file);
    PD_pv  = pv_data.ExtractedData.(tag_name);
    PD_rot = rot_data.ExtractedData.(tag_name);
    
    valid_frames = 1 : delta_step : length(PD_pv.Pos);
    
    % --- 数据容器 ---
    vel_glb = []; omg_glb = [];
    vel_blk = []; omg_blk = [];
    N_sum = 0; N_count = 0;
    
    for j = valid_frames
        pos_j = PD_pv.Pos{j};
        vel_j = PD_pv.Vel{j};
        omg_j = PD_rot.Omg{j};
        
        if isempty(vel_j) || isempty(omg_j) || isempty(pos_j), continue; end
        
        % 清洗 NaN
        nan_mask = any(isnan(vel_j), 2) | any(isnan(omg_j), 2) | any(isnan(pos_j), 2);
        pos_j(nan_mask, :) = [];
        vel_j(nan_mask, :) = [];
        omg_j(nan_mask, :) = [];
        
        % 1. 全局数据收集
        vel_glb = [vel_glb; vel_j];
        omg_glb = [omg_glb; omg_j];
        N_sum = N_sum + size(vel_j, 1);
        N_count = N_count + 1;
        
        % 2. 空间过滤：计算体区 (Bulk) 掩码
        % 假设坐标以 0 为中心。计算到中心的最大距离 (切比雪夫距离)
        dist_from_center = max(abs(pos_j), [], 2);
        bulk_mask = dist_from_center < half_L_bulk;
        
        % 收集 Bulk 数据
        vel_blk = [vel_blk; vel_j(bulk_mask, :)];
        omg_blk = [omg_blk; omg_j(bulk_mask, :)];
    end
    
    % --- 物理量计算 ---
    N_avg = N_sum / N_count;
    n_density = N_avg / V_box;
    Kn_vals(i) = (1 / (sqrt(2) * n_density * pi * d_particle^2)) / L_box;
    
    % 计算温度比值的内部函数 (复用逻辑)
    calc_T_ratio = @(v_mat, w_mat) (mean(sum((v_mat - mean(v_mat, 1)).^2, 2)) / 3) / ...
                                   ((2.0/5.0) * R_particle^2 * mean(sum(w_mat.^2, 2)) / 3);
    
    % 分别计算全局和体区的温度比
    if ~isempty(vel_glb)
        Tratio_glb(i) = calc_T_ratio(vel_glb, omg_glb);
    end
    if ~isempty(vel_blk)
        Tratio_blk(i) = calc_T_ratio(vel_blk, omg_blk);
    end
    
    fprintf('  Kn = %.3f | 全局 T_tr/T_rot = %.2f | 体区 T_tr/T_rot = %.2f\n', ...
        Kn_vals(i), Tratio_glb(i), Tratio_blk(i));
end

%% ========== 4. 学术级绘图 ==========

% --- 4.1 Kn vs φ (图1保持一致) ---
figure('Name', 'Knudsen Number', 'Position', [100, 200, 550, 500]);
loglog(phi_vals, Kn_vals, 'ko-', 'MarkerSize', 10, 'MarkerFaceColor', [0.2 0.5 0.8], 'LineWidth', 2);
hold on;
yline(1, '--r', 'Kn = 1', 'LineWidth', 1.5, 'FontSize', 14, 'LabelHorizontalAlignment', 'left');
phi_theory = logspace(log10(0.002), log10(0.1), 100);
Kn_theory = (1 ./ (sqrt(2) .* (phi_theory * V_box / (pi/6 * d_particle^3) / V_box) * pi * d_particle^2)) / L_box;
loglog(phi_theory, Kn_theory, 'b--', 'LineWidth', 1.5, 'DisplayName', 'Theory (ideal gas)');
xlabel('\phi', 'FontSize', 18); ylabel('Kn', 'FontSize', 18);
legend('Simulation', 'Theory', 'Location', 'northeast', 'FontSize', 14); legend boxoff;
set(gca, 'FontSize', 22, 'FontName', 'Times New Roman', 'LineWidth', 1.5);

% --- 4.2 T_trans / T_rot vs φ (核心绝杀图：全局 vs 体区) ---
figure('Name', 'Temperature Ratio', 'Position', [700, 200, 600, 500]);

% 画全局 (Global)
semilogx(phi_vals, Tratio_glb, 's--', 'Color', [0.5 0.5 0.5], 'MarkerSize', 8, ...
    'MarkerFaceColor', 'none', 'LineWidth', 1.5, 'DisplayName', 'Global (Whole Box)');
% hold on;
% 画体区 (Bulk)
semilogx(phi_vals, Tratio_blk, 'ro-', 'MarkerSize', 10, 'MarkerFaceColor', [0.8 0.2 0.2], ...
    'LineWidth', 2.5, 'DisplayName', sprintf('Bulk (Inner %.0f%%)', bulk_ratio*100));

yline(1, '--k', 'Equipartition', 'LineWidth', 1.5, 'FontSize', 14, 'LabelHorizontalAlignment', 'left');
xlabel('\phi', 'FontSize', 18); ylabel('T_{trans} / T_{rot}', 'FontSize', 18);
legend('Location', 'northeast', 'FontSize', 16); legend boxoff;
set(gca, 'FontSize', 22, 'FontName', 'Times New Roman', 'LineWidth', 1.5);

%% ========== 5. 打印汇总表格并保存 ==========
fprintf('\n========== 最终物理特性汇总 ==========\n');
fprintf('%10s %10s %15s %15s\n', 'phi', 'Kn', 'T_ratio(Global)', 'T_ratio(Bulk)');
fprintf('%s\n', repmat('-', 1, 55));
for i = 1:n_sets
    fprintf('%10.4f %10.4f %15.4f %15.4f\n', ...
        phi_vals(i), Kn_vals(i), Tratio_glb(i), Tratio_blk(i));
end

%% ========== 6. 解析无量纲穿透深度预测相变 (真实 DEM 数据版) ==========
fprintf('\n=== 正在生成基于真实 Kn 数据的相变预测图 ===\n');

figure('Name', 'Analytical Prediction of Phase Transition (Real Data)', 'Position', [200, 200, 750, 550]);

% 我们依然采用之前推断出的有效耗散常数 (包含极强的切向摩擦耗散)
epsilon_eff = 0.15; 
prefactor = 1.0; 

% --- 1. 计算理论理想曲线 (作为背景参考线) ---
phi_theory_dense = logspace(log10(0.002), log10(0.06), 200);
n_theory_dense = phi_theory_dense / (pi/6 * d_particle^3);
Kn_theory_dense = (1 ./ (sqrt(2) .* n_theory_dense * pi * d_particle^2)) / L_box;
xi_theory = prefactor .* (2 * Kn_theory_dense) ./ sqrt(epsilon_eff);

% --- 2. 【核心修改】直接使用前面老老实实算出来的真实 Kn_vals ! ---
% 这里的 Kn_vals 包含着系统的真实密度涨落和局域聚集效应
xi_sim_real = prefactor .* (2 * Kn_vals) ./ sqrt(epsilon_eff);

% --- 3. 作图 ---
% 背景理论线 (黑线)
loglog(phi_theory_dense, xi_theory, 'k-', 'LineWidth', 2.5, 'DisplayName', 'Analytical \xi (\epsilon_{eff} \approx 0.15)');
hold on;

% 【绝不作假的真实数据点】(蓝点)
loglog(phi_vals, xi_sim_real, 'bo', 'MarkerSize', 10, 'MarkerFaceColor', [0.3 0.6 0.9], ...
    'LineWidth', 2, 'DisplayName', 'Simulated \xi (Real Tracking)');

% 临界阈值
yline(1.0, 'r--', 'Critical Penetration (\xi = 1)', 'LineWidth', 2.5, 'FontSize', 16, 'LabelHorizontalAlignment', 'left', 'Color', 'r');

% 标注聚团区
patch_x = [0.02, 0.06, 0.06, 0.02];
patch_y = [0.1, 0.1, 1, 1];
patch(patch_x, patch_y, 'k', 'FaceAlpha', 0.1, 'EdgeColor', 'none', 'HandleVisibility', 'off');
text(0.025, 0.5, 'Clustering Zone (\xi < 1)', 'FontSize', 16, 'Color', [0.3 0.3 0.3]);

xlabel('\phi', 'Interpreter', 'latex', 'FontSize', 18);
ylabel('$\xi = 2Kn / \sqrt{\epsilon_{eff}}$', 'Interpreter', 'latex', 'FontSize', 18);

% title('Prediction of Clustering Instability (Real DEM Data)', 'FontSize', 18);

legend({'Analytical \xi', 'Simulated \xi (Real Tracking)'}, 'Location', 'northeast', 'FontSize', 14); 
legend boxoff;
set(gca, 'FontSize', 18, 'FontName', 'Times New Roman', 'LineWidth', 2);
ylim([0.1, 10]);
xlim([min(phi_theory_dense), max(phi_theory_dense)]);
% grid on;

fprintf('=== 真实数据相变图生成完毕！请检验蓝点是否依然能命中交点！ ===\n');

%% 保存数据
save_path='/media/gezhuan/M78/Space_Active/Data/Maxwell_Boltz/data/Kn_Tr_Tt';
if ~exist(save_path, 'dir'), mkdir(save_path); end
cd(save_path);
save_name = sprintf('kn_Tr_Tt_bulk%.1f_delta%d.mat', bulk_ratio, delta_step);
save(save_name);
fprintf('\n✅ 数据已保存至: %s\n', fullfile(save_path, save_name));