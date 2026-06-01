%% =====================================================================
%  Knudsen 数 & 平动/转动温度比 随体积分数变化
%  输入: Pos_Vel 数据 (.mat) + Rotation 数据 (.mat)
%  输出: Kn(φ) 和 T_trans/T_rot(φ) 的图
% =====================================================================

clear all; close all; clc;

%% ========== 1. 用户配置区 ==========

% --- 1.1 数据路径 & 通配符模式 ---
pos_vel_dir = 'G:\Space_Active\Data\Maxwell_Boltz\data\Pos_Vel_data\Nofric';
rot_dir     = 'G:\Space_Active\Data\Maxwell_Boltz\data\Rotation_data\Nofric';

% 通配符匹配模式 (只需改这里，脚本会自动配对)
pv_pattern  = 'Data_*.mat';
rot_pattern = 'RotData_*.mat';

% --- 1.2 系统参数 ---
tag_name     = 'Tag_2';        % 颗粒字段名
d_particle   = 0.08;           % 颗粒直径 [cm]
R_particle   = d_particle / 2; % 颗粒半径 [cm]
L_box        = 4.0;            % 容器边长 [cm] (特征长度)
V_box        = L_box^3;        % 容器体积 [cm³]

% 3D 硬球平均自由程: λ = 1 / (√2 · n · π · d²)
% Knudsen 数: Kn = λ / L

% 温度比: T_trans/T_rot = m<δv²> / (I<ω²>)
%   对均匀实心球 I = (2/5)mR², 质量约掉:
%   T_trans/T_rot = <δv²> / ((2/5) R² <ω²>)

% --- 1.3 采样设置 ---
delta_step = 5;  % 每隔多少帧取一帧 (设大可加速)

%% ========== 2. 自动检索与配对数据文件 ==========

% 检索位置速度文件
pv_files  = dir(fullfile(pos_vel_dir, pv_pattern));
rot_files = dir(fullfile(rot_dir, rot_pattern));

if isempty(pv_files)
    error('在 %s 中未找到匹配 "%s" 的文件！', pos_vel_dir, pv_pattern);
end
if isempty(rot_files)
    error('在 %s 中未找到匹配 "%s" 的文件！', rot_dir, rot_pattern);
end

% --- 辅助函数: 从文件名中提取 phi 值 ---
% 匹配 "phi0.01" 或 "phi0.004" 这样的模式
extract_phi = @(fname) str2double(regexp(fname, '(?<=phi)\d+\.?\d*', 'match', 'once'));

% 构建 pos_vel 文件的 phi -> filename 映射
pv_map = containers.Map('KeyType', 'double', 'ValueType', 'char');
for k = 1:length(pv_files)
    phi_k = extract_phi(pv_files(k).name);
    if ~isnan(phi_k)
        % 如果同一个 phi 有多个文件，取最新的 (dir 按名称排序，时间戳大的在后)
        pv_map(phi_k) = pv_files(k).name;
    end
end

% 构建 rot 文件的 phi -> filename 映射
rot_map = containers.Map('KeyType', 'double', 'ValueType', 'char');
for k = 1:length(rot_files)
    phi_k = extract_phi(rot_files(k).name);
    if ~isnan(phi_k)
        rot_map(phi_k) = rot_files(k).name;
    end
end

% 找到两边都有的 phi 值，取交集并排序
phi_pv  = cell2mat(pv_map.keys);
phi_rot = cell2mat(rot_map.keys);
phi_common = intersect(phi_pv, phi_rot);
phi_common = sort(phi_common);

if isempty(phi_common)
    error('位置速度数据和转动数据没有匹配的 phi 值！请检查文件命名。');
end

% 组装最终的 datasets 列表
n_sets = length(phi_common);
fprintf('自动匹配到 %d 组数据:\n', n_sets);
datasets = cell(n_sets, 4);  % {pv_file, rot_file, phi, label}
for k = 1:n_sets
    phi_k = phi_common(k);
    datasets{k, 1} = pv_map(phi_k);
    datasets{k, 2} = rot_map(phi_k);
    datasets{k, 3} = phi_k;
    datasets{k, 4} = sprintf('\\phi=%.4g', phi_k);
    fprintf('  φ = %-8.4g | PV: %s\n', phi_k, datasets{k,1});
    fprintf('  %14s | Rot: %s\n', '', datasets{k,2});
end

%% ========== 3. 批量计算 ==========

phi_vals   = zeros(n_sets, 1);
Kn_vals    = zeros(n_sets, 1);
Tratio_vals = zeros(n_sets, 1);

% 额外诊断量
N_particles_avg = zeros(n_sets, 1);
T_trans_vals    = zeros(n_sets, 1);  % 平动温度 (任意单位, ∝ m)
T_rot_vals      = zeros(n_sets, 1);  % 转动温度 (任意单位, ∝ I)

fprintf('=== 开始计算 Kn 和 T_trans/T_rot ===\n');

parfor i = 1:n_sets
    pv_file  = fullfile(pos_vel_dir, datasets{i, 1});
    rot_file = fullfile(rot_dir,     datasets{i, 2});
    phi_val  = datasets{i, 3};
    phi_vals(i) = phi_val;
    
    fprintf('\n[%d/%d] φ = %.4f\n', i, n_sets, phi_val);
    
    % --- 3.1 加载数据 ---
    pv_data  = load(pv_file);
    rot_data = load(rot_file);
    
    PD_pv  = pv_data.ExtractedData.(tag_name);
    PD_rot = rot_data.ExtractedData.(tag_name);
    
    total_frames = length(PD_pv.Pos);
    valid_frames = 1 : delta_step : total_frames;
    
    % --- 3.2 汇总速度和角速度 ---
    vel_all = [];
    omg_all = [];
    N_sum = 0;
    N_count = 0;
    
    for j = valid_frames
        vel_j = PD_pv.Vel{j};
        omg_j = PD_rot.Omg{j};
        
        if isempty(vel_j) || isempty(omg_j)
            continue;
        end
        
        vel_all = [vel_all; vel_j];   % [N_total x 3], cm/s
        omg_all = [omg_all; omg_j];   % [N_total x 3], rad/s
        
        N_sum   = N_sum + size(vel_j, 1);
        N_count = N_count + 1;
    end
    
    N_avg = N_sum / N_count;  % 平均每帧颗粒数
    N_particles_avg(i) = N_avg;
    
    % 数据清洗
    nan_mask = any(isnan(vel_all), 2) | any(isnan(omg_all), 2);
    vel_all(nan_mask, :) = [];
    omg_all(nan_mask, :) = [];
    
    % --- 3.3 Knudsen 数 ---
    n_density = N_avg / V_box;                        % 数密度 [1/cm³]
    lambda_mfp = 1 / (sqrt(2) * n_density * pi * d_particle^2);  % 平均自由程 [cm]
    Kn_vals(i) = lambda_mfp / L_box;
    
    fprintf('  N_avg = %.1f,  n = %.2e /cm³,  λ = %.3f cm,  Kn = %.3f\n', ...
        N_avg, n_density, lambda_mfp, Kn_vals(i));
    
    % --- 3.4 平动温度 (去除质心漂移) ---
    v_mean = mean(vel_all, 1);  % 质心速度 [1x3]
    dv = vel_all - v_mean;      % 热运动速度
    v2_mean = mean(sum(dv.^2, 2));  % <δv²> [cm²/s²]
    
    % 平动粒度温度 (每自由度): T_trans_dof = (1/2) m <δv_α²>
    % 这里只算无量纲比值，不需要 m 的绝对值
    T_trans = v2_mean / 3;  % 单自由度平动温度 ∝ m (此处省略 m)
    T_trans_vals(i) = T_trans;
    
    % --- 3.5 转动温度 ---
    omg2_mean = mean(sum(omg_all.^2, 2));  % <ω²> [rad²/s²]
    
    % 转动粒度温度 (每自由度): T_rot_dof = (1/2) I <ω_α²>
    % I = (2/5) m R²  =>  T_rot_dof ∝ m·R²
    % 无量纲比: T_trans/T_rot = (m <δv²>/3) / ((2/5)mR² <ω²>/3) = <δv²>/((2/5)R²<ω²>)
    T_rot = (2.0/5.0) * R_particle^2 * omg2_mean / 3;  % 单自由度转动温度 ∝ m
    T_rot_vals(i) = T_rot;
    
    Tratio_vals(i) = T_trans / T_rot;
    
    fprintf('  <δv²> = %.4e cm²/s²,  <ω²> = %.4e rad²/s²\n', v2_mean, omg2_mean);
    fprintf('  T_trans/T_rot = %.4f\n', Tratio_vals(i));
end

%% ========== 4. 绘图 ==========

% --- 4.1 Kn vs φ ---
figure('Name', 'Knudsen Number vs Volume Fraction', ...
    'Position', [100, 400, 600, 500]);

loglog(phi_vals, Kn_vals, 'ko-', 'MarkerSize', 10, 'MarkerFaceColor', [0.2 0.5 0.8], ...
    'LineWidth', 2);
hold on;

% 画 Kn = 1 参考线
xline(min(phi_vals)*0.5, 'HandleVisibility', 'off'); % dummy for axis
yline(1, '--r', 'Kn = 1', 'LineWidth', 1.5, 'FontSize', 14, 'LabelHorizontalAlignment', 'left');

% 理论曲线: Kn(φ) = V_box / (√2·N·π·d²·L) ，其中 N = φ·V_box / (π/6·d³)
phi_theory = logspace(log10(0.002), log10(0.1), 100);
N_theory = phi_theory * V_box / (pi/6 * d_particle^3);
n_theory = N_theory / V_box;
lambda_theory = 1 ./ (sqrt(2) .* n_theory * pi * d_particle^2);
Kn_theory = lambda_theory / L_box;
loglog(phi_theory, Kn_theory, 'b--', 'LineWidth', 1.5, 'DisplayName', 'Theory (ideal gas)');

xlabel('\phi', 'FontSize', 18);
ylabel('Kn', 'FontSize', 18);
% title('Knudsen Number vs Volume Fraction', 'FontSize', 16);
legend('Simulation', 'Theory (ideal gas)', 'Location', 'northeast', 'FontSize', 14);
legend boxoff;
set(gca, 'FontSize', 24, 'FontName', 'Times New Roman', 'LineWidth', 1.5);


% --- 4.2 T_trans / T_rot vs φ ---
figure('Name', 'Temperature Ratio vs Volume Fraction', ...
    'Position', [750, 400, 600, 500]);

semilogx(phi_vals, Tratio_vals, 'ks-', 'MarkerSize', 10, 'MarkerFaceColor', [0.8 0.3 0.3], ...
    'LineWidth', 2);
hold on;

% 能量均分参考线 T_trans/T_rot = 1
yline(1, '--k', 'Equipartition', 'LineWidth', 1.5, 'FontSize', 14, ...
    'LabelHorizontalAlignment', 'left');

xlabel('\phi', 'FontSize', 18);
ylabel('T_{trans} / T_{rot}', 'FontSize', 18);
% title('Translational / Rotational Temperature Ratio', 'FontSize', 16);
set(gca, 'FontSize', 24, 'FontName', 'Times New Roman', 'LineWidth', 1.5);


% --- 4.3 打印汇总表格 ---
fprintf('\n========== 汇总结果 ==========\n');
fprintf('%10s %8s %10s %12s %12s %12s\n', 'phi', 'N_avg', 'Kn', 'T_trans', 'T_rot', 'T_tr/T_rot');
fprintf('%s\n', repmat('-', 1, 70));
for i = 1:n_sets
    fprintf('%10.4f %8.1f %10.4f %12.4e %12.4e %12.4f\n', ...
        phi_vals(i), N_particles_avg(i), Kn_vals(i), ...
        T_trans_vals(i), T_rot_vals(i), Tratio_vals(i));
end

fprintf('\n✅ 计算完成！\n');

%%
save_path='/media/gezhuan/M78/Space_Active/Data/Maxwell_Boltz/data/Kn_Tr_Tt';
cd (save_path)
save('kn_Tr_Tt_delta5.mat');