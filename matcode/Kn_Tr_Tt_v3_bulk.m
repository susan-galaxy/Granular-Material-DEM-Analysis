%% =====================================================================
%  Knudsen 数 & 平动/转动温度比 (真实局域密度进化版)
%  核心修正: 使用 Bulk 内部真实的颗粒数密度来计算局域 Kn 与穿透深度
% =====================================================================
clear all; close all; clc;

%% ========== 1. 用户配置区 ==========
pos_vel_dir = '/media/gezhuan/M78/Space_Active/Data/Active_Gas/Data/Pos_Vel';
rot_dir     = '/media/gezhuan/M78/Space_Active/Data/Active_Gas/Data/Rotations';

pv_pattern  = 'Data*.mat';
rot_pattern = 'RotData*.mat';

tag_name     = 'Tag_2';        
d_particle   = 0.08;           % 颗粒直径 [cm]
R_particle   = d_particle / 2; 
L_box        = 4.0;            % 容器边长 [cm]
V_box        = L_box^3;        % 全局体积 [cm³]

delta_step = 5;  
bulk_ratio = 0.6;  % 取中心 60% 区域
half_L_bulk = (L_box * bulk_ratio) / 2; 
V_bulk = (L_box * bulk_ratio)^3; % 【新增】体区真实体积

fprintf('=== 启动物理特性计算 (Bulk Ratio = %.2f) ===\n', bulk_ratio);

%% ========== 2. 自动检索与配对 ==========
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

%% ========== 3. 批量计算 (全局 vs 真实体区) ==========
phi_vals    = zeros(n_sets, 1);
Kn_glb_vals = zeros(n_sets, 1); % 假想的均匀全局 Kn
Kn_blk_vals = zeros(n_sets, 1); % 【核心】真实的局域 Bulk Kn

Tratio_glb = zeros(n_sets, 1);
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
    
    vel_glb = []; omg_glb = [];
    vel_blk = []; omg_blk = [];
    N_glb_sum = 0; N_blk_sum = 0; N_count = 0;
    
    for j = valid_frames
        pos_j = PD_pv.Pos{j};
        vel_j = PD_pv.Vel{j};
        omg_j = PD_rot.Omg{j};
        
        if isempty(vel_j) || isempty(omg_j) || isempty(pos_j), continue; end
        
        nan_mask = any(isnan(vel_j), 2) | any(isnan(omg_j), 2) | any(isnan(pos_j), 2);
        pos_j(nan_mask, :) = [];
        vel_j(nan_mask, :) = [];
        omg_j(nan_mask, :) = [];
        
        % 1. 全局数据收集
        vel_glb = [vel_glb; vel_j];
        omg_glb = [omg_glb; omg_j];
        N_glb_sum = N_glb_sum + size(vel_j, 1);
        
        % 2. Bulk 空间过滤
        dist_from_center = max(abs(pos_j), [], 2);
        bulk_mask = dist_from_center < half_L_bulk;
        
        vel_blk_current = vel_j(bulk_mask, :);
        vel_blk = [vel_blk; vel_blk_current];
        omg_blk = [omg_blk; omg_j(bulk_mask, :)];
        
        % 【核心】：统计真实落在 Bulk 内部的颗粒数
        N_blk_sum = N_blk_sum + size(vel_blk_current, 1);
        N_count = N_count + 1;
    end
    
    % --- 计算假想的全局 Kn (理论直线) ---
    N_glb_avg = N_glb_sum / N_count;
    n_density_glb = N_glb_avg / V_box;
    Kn_glb_vals(i) = (1 / (sqrt(2) * n_density_glb * pi * d_particle^2)) / L_box;
    
    % --- 【核心】：计算真实的局域 Bulk Kn ---
    N_blk_avg = N_blk_sum / N_count;
    n_density_blk = N_blk_avg / V_bulk;
    lambda_blk = 1 / (sqrt(2) * n_density_blk * pi * d_particle^2);
    Kn_blk_vals(i) = lambda_blk / L_box;
    
    calc_T_ratio = @(v_mat, w_mat) (mean(sum((v_mat - mean(v_mat, 1)).^2, 2)) / 3) / ...
                                   ((2.0/5.0) * R_particle^2 * mean(sum(w_mat.^2, 2)) / 3);
    
    if ~isempty(vel_glb), Tratio_glb(i) = calc_T_ratio(vel_glb, omg_glb); end
    if ~isempty(vel_blk), Tratio_blk(i) = calc_T_ratio(vel_blk, omg_blk); end
    
    fprintf('  理论Kn = %.3f | 真实局域Kn = %.3f\n', Kn_glb_vals(i), Kn_blk_vals(i));
end

%% ========== 4. 学术级绘图：揭示真实的 Kn 偏离 ==========
figure('Name', 'Real Local Knudsen Number', 'Position', [100, 200, 550, 500]);
% 画理论理想气体直线 (基于全局密度)
phi_theory = logspace(log10(0.002), log10(0.1), 100);
Kn_theory = (1 ./ (sqrt(2) .* (phi_theory / (pi/6 * d_particle^3)) * pi * d_particle^2)) / L_box;
loglog(phi_theory, Kn_theory, 'k--', 'LineWidth', 2, 'DisplayName', 'Theory (Uniform Gas)');
hold on;
% 画真实的体区局域 Kn
loglog(phi_vals, Kn_blk_vals, 'bo-', 'MarkerSize', 10, 'MarkerFaceColor', [0.2 0.6 0.9], 'LineWidth', 2, 'DisplayName', 'Simulated Local Kn (Bulk)');

yline(1, '--r', 'Kn = 1', 'LineWidth', 1.5, 'FontSize', 14, 'LabelHorizontalAlignment', 'left');
xlabel('Volume Fraction \phi', 'FontSize', 18); ylabel('Knudsen Number Kn', 'FontSize', 18);
title('Breakdown of Uniformity', 'FontSize', 16);
legend('Location', 'northeast', 'FontSize', 14); legend boxoff;
set(gca, 'FontSize', 18, 'FontName', 'Times New Roman', 'LineWidth', 1.5);

%% ========== 5. T_trans / T_rot vs φ (核心绝杀图) ==========
figure('Name', 'Temperature Ratio', 'Position', [700, 200, 600, 500]);
semilogx(phi_vals, Tratio_glb, 's--', 'Color', [0.5 0.5 0.5], 'MarkerSize', 8, 'MarkerFaceColor', 'none', 'LineWidth', 1.5, 'DisplayName', 'Global (Whole Box)');
hold on;
semilogx(phi_vals, Tratio_blk, 'ro-', 'MarkerSize', 10, 'MarkerFaceColor', [0.8 0.2 0.2], 'LineWidth', 2.5, 'DisplayName', sprintf('Bulk (Inner %.0f%%)', bulk_ratio*100));

yline(1, '--k', 'Equipartition', 'LineWidth', 1.5, 'FontSize', 14, 'LabelHorizontalAlignment', 'left');
xlabel('Volume Fraction \phi', 'FontSize', 18); ylabel('T_{trans} / T_{rot}', 'FontSize', 18);
legend('Location', 'northeast', 'FontSize', 16); legend boxoff;
set(gca, 'FontSize', 18, 'FontName', 'Times New Roman', 'LineWidth', 1.5);

%% ========== 6. 终极解析：基于真实局域 Kn 的无量纲衰减长度 ==========
figure('Name', 'Analytical Prediction (Real Local Data)', 'Position', [300, 150, 750, 550]);

epsilon_eff = 0.15; 
prefactor = 1.0; 

% 1. 理论参考线 (假设绝对均匀)
xi_theory = prefactor .* (2 * Kn_theory) ./ sqrt(epsilon_eff);
loglog(phi_theory, xi_theory, 'k-', 'LineWidth', 2.5, 'DisplayName', 'Analytical \xi (Uniform Assumption)');
hold on;

% 2. 真实追踪的局域穿透深度 (使用刚刚算出来的真实 Kn_blk_vals)
xi_sim_real = prefactor .* (2 * Kn_blk_vals) ./ sqrt(epsilon_eff);
loglog(phi_vals, xi_sim_real, 'bo-', 'MarkerSize', 10, 'MarkerFaceColor', [0.3 0.6 0.9], ...
    'LineWidth', 2.5, 'DisplayName', 'Simulated \xi (Real Local Density)');

% 临界阈值
yline(1.0, 'r--', 'Critical Penetration (\xi = 1)', 'LineWidth', 2.5, 'FontSize', 16, 'LabelHorizontalAlignment', 'left');

% 标注聚团区
patch_x = [0.015, 0.1, 0.1, 0.015]; % 预估临界点在0.015左右
patch_y = [0.01, 0.01, 1, 1];
patch(patch_x, patch_y, 'k', 'FaceAlpha', 0.1, 'EdgeColor', 'none', 'HandleVisibility', 'off');
text(0.02, 0.3, 'Clustering Zone (\xi < 1)', 'FontSize', 16, 'Color', [0.3 0.3 0.3]);

xlabel('Volume Fraction \phi', 'Interpreter', 'latex', 'FontSize', 18);
ylabel('Dimensionless Decay Length \xi = 2Kn_{local} / \sqrt{\epsilon_{eff}}', 'Interpreter', 'latex', 'FontSize', 18);
title('Prediction of Clustering Instability (Real Local Density)', 'FontSize', 18);

legend('Location', 'northeast', 'FontSize', 14); legend boxoff;
set(gca, 'FontSize', 18, 'FontName', 'Times New Roman', 'LineWidth', 2);
ylim([0.1, 10]); xlim([0.002, 0.1]); grid on;