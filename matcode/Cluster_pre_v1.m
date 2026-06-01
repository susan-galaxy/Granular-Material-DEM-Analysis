%% =====================================================================
%  三维颗粒气体温度穿透深度 λ 与相变预测
%  理论基础: Granular Media 教材 §5.4.4 公式 (5.106)-(5.108)
%  
%  核心公式:
%    能量方程: d/dz [ K dT/dz ] = Γ
%    热导率:   K = ρ_p d F4(φ) √T  ∝  ρ_p d √T  (稀薄极限)
%    耗散率:   Γ = (ρ_p/d)(1-e²) F5(φ) T^{3/2}  ∝  (ρ_p/d)(1-e²) φ² T^{3/2}
%    穿透深度: λ ~ √(KT/Γ) ≃ 3SL / (Nπd²√(1-e²))     [1D, 公式5.107]
%    临界粒子数: Nc ≃ 3S / (πd²√(1-e²))                [公式5.108]
%
%  3D 扩展:
%    λ_3D ≈ d / (6φ χ(φ) √(1-e²))     [用平均自由程表达]
%    相变条件: λ_3D < L/2
%
%  同时从模拟数据中直接测量壳层温度剖面并拟合 λ_measured
% =====================================================================

clear all; close all; clc;

%% ========== 1. 系统参数 ==========

L     = 4.0;          % 容器边长 [cm]
d     = 0.08;         % 颗粒直径 [cm] (0.8 mm)
e     = 0.9;          % 法向恢复系数 (亚克力典型值, 可修改)
S     = L^2;          % 壁面面积 [cm²] (一面)
V_box = L^3;          % 容器体积 [cm³]
V_p   = (pi/6)*d^3;   % 单颗粒体积 [cm³]

fprintf('=== 温度穿透深度 λ 与相变预测 ===\n');
fprintf('  容器: %.1f × %.1f × %.1f cm³\n', L, L, L);
fprintf('  颗粒直径: %.2f cm (%.1f mm)\n', d, d*10);
fprintf('  恢复系数: e = %.2f, 1-e² = %.4f\n', e, 1-e^2);

%% ========== 2. 理论计算: λ(φ) 和 Nc ==========

% 扫描体积分数范围
phi_theory = logspace(-4, -0.3, 500)';  % φ 从 0.0001 到 ~0.5

% --- 方法 A: 直接套用教材公式 (5.107) 的 3D 版本 ---
% 教材原式: λ = 3SL / (Nπd²√(1-e²))
% 其中 N = φ V_box / V_p = 6φL³/(πd³)
% 代入得: λ = 3·S·L·πd³ / (6φL³·πd²·√(1-e²))
%          = 3·L²·d / (6φL³·√(1-e²))
%          = d / (2φL·√(1-e²))
% 
% 但这是 1D 的结果 (沿 z 方向). 在 3D 立方体中, 能量从 6 面壁注入,
% 每面壁"负责"的穿透距离是 L/2. 有效公式相同, 只是判据变成 λ vs L/2.

lambda_textbook = d ./ (2 * phi_theory * L * sqrt(1 - e^2));

% --- 方法 B: 用平均自由程的物理图像 ---
% λ_mfp = 1/(n σ_coll χ) = d/(6φχ)  (Enskog 修正)
% χ(φ) = (1 - φ/2)/(1-φ)³           (Carnahan-Starling)
% 每次碰撞损失 (1-e²) 的动能
% 穿透深度 = λ_mfp / √(1-e²)  (因为需要 ~(1-e²)^{-1} 次碰撞才耗尽能量)

chi_CS = (1 - phi_theory/2) ./ (1 - phi_theory).^3;
lambda_mfp = d ./ (6 * phi_theory .* chi_CS);
lambda_kinetic = lambda_mfp ./ sqrt(1 - e^2);

% --- 临界体积分数 (λ = L/2 的解) ---
% 方法 A: d/(2φ_c L √(1-e²)) = L/2  →  φ_c = d/(L²(1-e²))
phi_c_textbook = d / (L^2 * (1 - e^2));

% 方法 B: 数值求解 λ_kinetic(φ_c) = L/2
[~, idx_c] = min(abs(lambda_kinetic - L/2));
phi_c_kinetic = phi_theory(idx_c);

% 教材临界粒子数 (5.108): Nc = 3S/(πd²√(1-e²))
% 在 3D 立方体中, 我们从 6 面壁考虑, 每面贡献 S = L²
Nc_textbook = 3 * S / (pi * d^2 * sqrt(1 - e^2));
% 对应的等效 φ_c:
phi_c_from_Nc = Nc_textbook * V_p / V_box;

fprintf('\n--- 理论预测 ---\n');
fprintf('  教材公式 φ_c (方法A) = %.6f\n', phi_c_textbook);
fprintf('  动力学公式 φ_c (方法B) = %.6f\n', phi_c_kinetic);
fprintf('  教材 Nc (5.108, 单面壁) = %.0f\n', Nc_textbook);
fprintf('  Nc 对应等效 φ = %.6f\n', phi_c_from_Nc);

%% ========== 3. 从模拟数据测量壳层温度剖面 ==========

% --- 数据配置 ---
pos_vel_dir = 'G:\Space_Active\Data\Maxwell_Boltz\data\Pos_Vel_data';
rot_dir     = 'G:\Space_Active\Data\Maxwell_Boltz\data\Rotation_data\Rotation_data';
pv_pattern  = 'Data_FixSpace_R0_04_mode0td_phi*.mat';
rot_pattern = 'RotData_FixSpace_R0_04_mode0td_phi*.mat';
tag_name    = 'Tag_2';
delta_step  = 5;

% --- 壳层设置 ---
n_shells = 20;
half_L = L / 2;
shell_edges = linspace(0, half_L, n_shells + 1);
shell_centers = (shell_edges(1:end-1) + shell_edges(2:end)) / 2;
% 转换为"距壁面距离" (教材的 z 坐标从壁面算起)
dist_from_wall = half_L - shell_centers;  % 壁面=0, 中心=L/2

% --- 自动检索文件 ---
pv_files = dir(fullfile(pos_vel_dir, pv_pattern));
rot_files = dir(fullfile(rot_dir, rot_pattern));

extract_phi = @(fname) str2double(regexp(fname, '(?<=phi)\d+\.?\d*', 'match', 'once'));

pv_map = containers.Map('KeyType', 'double', 'ValueType', 'char');
rot_map = containers.Map('KeyType', 'double', 'ValueType', 'char');

for k = 1:length(pv_files)
    phi_k = extract_phi(pv_files(k).name);
    if ~isnan(phi_k), pv_map(phi_k) = pv_files(k).name; end
end
for k = 1:length(rot_files)
    phi_k = extract_phi(rot_files(k).name);
    if ~isnan(phi_k), rot_map(phi_k) = rot_files(k).name; end
end

phi_data = sort(cell2mat(pv_map.keys));
n_data = length(phi_data);

fprintf('\n  找到 %d 个模拟数据集: ', n_data);
fprintf('%.4f  ', phi_data); fprintf('\n');

% --- 存储 ---
Measured = struct();
Measured.phi = phi_data;
Measured.shell_centers = shell_centers;
Measured.dist_from_wall = dist_from_wall;
Measured.T_tr_profile = zeros(n_data, n_shells);
Measured.T_rot_profile = zeros(n_data, n_shells);
Measured.Tx_profile = zeros(n_data, n_shells);  % 各分量
Measured.Ty_profile = zeros(n_data, n_shells);
Measured.Tz_profile = zeros(n_data, n_shells);
Measured.n_profile = zeros(n_data, n_shells);
Measured.lambda_fit_tr = zeros(n_data, 1);
Measured.lambda_fit_rot = zeros(n_data, 1);
Measured.T_wall_fit = zeros(n_data, 1);
Measured.T_center_fit = zeros(n_data, 1);

%% ========== 4. 逐 φ 提取壳层温度剖面 ==========

for i = 1:n_data
    phi_val = phi_data(i);
    fprintf('\n[%d/%d] φ = %.4f ...\n', i, n_data, phi_val);
    
    % 加载平动数据
    pv_file = fullfile(pos_vel_dir, pv_map(phi_val));
    loaded = load(pv_file);
    PD = loaded.ExtractedData.(tag_name);
    
    total_frames = length(PD.Pos);
    valid_frames = 1:delta_step:total_frames;
    
    % 加载转动数据 (如果有)
    has_rot = rot_map.isKey(phi_val);
    if has_rot
        rot_file = fullfile(rot_dir, rot_map(phi_val));
        loaded_rot = load(rot_file);
        PD_rot = loaded_rot.ExtractedData.(tag_name);
        rot_steps = PD_rot.TimeStep;
    end
    
    % 壳层累加器
    shell_vx2_sum = zeros(n_shells, 1);
    shell_vy2_sum = zeros(n_shells, 1);
    shell_vz2_sum = zeros(n_shells, 1);
    shell_omg2_sum = zeros(n_shells, 1);
    shell_count = zeros(n_shells, 1);
    shell_rot_count = zeros(n_shells, 1);
    
    for j_idx = 1:length(valid_frames)
        j = valid_frames(j_idx);
        
        pos_j = PD.Pos{j};
        vel_j = PD.Vel{j};
        if isempty(pos_j), continue; end
        
        % 清洗
        bad = any(isnan(pos_j), 2) | any(isnan(vel_j), 2);
        pos_j(bad, :) = [];
        vel_j(bad, :) = [];
        if isempty(pos_j), continue; end
        
        % 去宏观流
        mean_v = mean(vel_j, 1);
        vf = vel_j - mean_v;
        
        % Chebyshev 距离 (到中心)
        cheby = max(abs(pos_j), [], 2);
        s_idx = discretize(cheby, shell_edges);
        
        for s = 1:n_shells
            mask = (s_idx == s);
            if any(mask)
                shell_vx2_sum(s) = shell_vx2_sum(s) + sum(vf(mask, 1).^2);
                shell_vy2_sum(s) = shell_vy2_sum(s) + sum(vf(mask, 2).^2);
                shell_vz2_sum(s) = shell_vz2_sum(s) + sum(vf(mask, 3).^2);
                shell_count(s) = shell_count(s) + sum(mask);
            end
        end
        
        % 转动数据
        if has_rot
            rot_idx = find(rot_steps == j, 1);
            if ~isempty(rot_idx) && ~isempty(PD_rot.Omg{rot_idx})
                omg_j = PD_rot.Omg{rot_idx};
                % 同步清洗
                if size(omg_j, 1) == size(pos_j, 1) + sum(bad)
                    omg_j(bad, :) = [];
                end
                if size(omg_j, 1) == size(pos_j, 1)
                    omg2 = sum(omg_j.^2, 2);
                    for s = 1:n_shells
                        mask = (s_idx == s);
                        if any(mask)
                            shell_omg2_sum(s) = shell_omg2_sum(s) + sum(omg2(mask));
                            shell_rot_count(s) = shell_rot_count(s) + sum(mask);
                        end
                    end
                end
            end
        end
    end
    
    % 计算壳层温度 (T = <v²>, 不含 1/2 和 m, 与教材一致)
    for s = 1:n_shells
        if shell_count(s) > 0
            Measured.Tx_profile(i, s) = shell_vx2_sum(s) / shell_count(s);
            Measured.Ty_profile(i, s) = shell_vy2_sum(s) / shell_count(s);
            Measured.Tz_profile(i, s) = shell_vz2_sum(s) / shell_count(s);
            Measured.T_tr_profile(i, s) = (shell_vx2_sum(s) + shell_vy2_sum(s) + shell_vz2_sum(s)) / (3 * shell_count(s));
            Measured.n_profile(i, s) = shell_count(s) / length(valid_frames);
        end
        if shell_rot_count(s) > 0
            Measured.T_rot_profile(i, s) = shell_omg2_sum(s) / (3 * shell_rot_count(s));
        end
    end
    
    % --- 拟合穿透深度 λ ---
    % 教材模型: T(z) ∝ cosh(z/λ) 或近似为指数衰减
    % 更精确: 对称边界条件下 T(r) 的解形如 cosh((r - L/2)/λ) / cosh(L/(2λ))
    % 但简单起见, 先做 log(T) vs dist_from_wall 的线性拟合
    
    T_prof = Measured.T_tr_profile(i, :)';
    valid_s = T_prof > 0 & ~isnan(T_prof);
    
    if sum(valid_s) >= 4
        % 用距壁面的距离 (壁面=0, 越远越深入)
        % 模型: T(z) = T_wall * exp(-z/λ)  其中 z = dist_from_wall
        % 但注意: dist_from_wall 是从壁面向内的距离
        % shell_centers 是从中心向外, 所以:
        % dist_from_wall = half_L - shell_centers (壁面=0)
        % 我们期待: T 随 dist_from_wall 增大(深入内部) 而减小
        
        z_wall = dist_from_wall(valid_s)';  % 从壁面算的距离
        T_valid = T_prof(valid_s);
        
        % 方法1: 简单指数拟合 T = T0 * exp(-z/λ)
        try
            % 非线性拟合更稳健
            ft = fittype('a * exp(-x/b)', 'independent', 'x');
            opts = fitoptions(ft);
            opts.StartPoint = [max(T_valid), L/4];
            opts.Lower = [0, 0.001];
            opts.Upper = [10*max(T_valid), 10*L];
            [fit_result, gof] = fit(z_wall, T_valid, ft, opts);
            
            Measured.lambda_fit_tr(i) = fit_result.b;
            Measured.T_wall_fit(i) = fit_result.a;
            Measured.T_center_fit(i) = fit_result.a * exp(-half_L / fit_result.b);
            
            fprintf('  拟合: T_wall=%.4f, λ_tr=%.4f cm, R²=%.4f\n', ...
                fit_result.a, fit_result.b, gof.rsquare);
        catch
            % Fallback: 线性拟合 log(T)
            p = polyfit(z_wall, log(T_valid), 1);
            if p(1) < 0
                Measured.lambda_fit_tr(i) = -1/p(1);
            else
                Measured.lambda_fit_tr(i) = Inf;
            end
            fprintf('  线性拟合: λ_tr=%.4f cm (fallback)\n', Measured.lambda_fit_tr(i));
        end
    end
    
    % 转动温度的 λ
    T_rot_prof = Measured.T_rot_profile(i, :)';
    valid_r = T_rot_prof > 0 & ~isnan(T_rot_prof);
    if sum(valid_r) >= 4
        try
            ft = fittype('a * exp(-x/b)', 'independent', 'x');
            opts = fitoptions(ft);
            opts.StartPoint = [max(T_rot_prof(valid_r)), L/4];
            opts.Lower = [0, 0.001];
            opts.Upper = [10*max(T_rot_prof(valid_r)), 10*L];
            [fit_rot, ~] = fit(dist_from_wall(valid_r)', T_rot_prof(valid_r), ft, opts);
            Measured.lambda_fit_rot(i) = fit_rot.b;
            fprintf('  转动穿透深度: λ_rot=%.4f cm\n', fit_rot.b);
        catch
            Measured.lambda_fit_rot(i) = Inf;
        end
    end
end

%% ========== 5. 绘图 ==========

cmap = lines(n_data);

% ------ Fig 1: 理论 λ(φ) 曲线 + 模拟测量值 ------
figure('Name', 'Fig.1: Penetration Depth λ vs φ', ...
    'Position', [100, 100, 800, 550]);

loglog(phi_theory, lambda_textbook, 'b-', 'LineWidth', 2, ...
    'DisplayName', '\lambda_{textbook} = d/(2\phi L\surd(1-e^2))');
hold on;
loglog(phi_theory, lambda_kinetic, 'r--', 'LineWidth', 2, ...
    'DisplayName', '\lambda_{kinetic} (Enskog)');

% 容器半长参考线
yline(half_L, 'k:', 'LineWidth', 2, 'HandleVisibility', 'off');
text(phi_theory(end)*0.3, half_L*1.3, 'L/2 = 2 cm', 'FontSize', 14, 'Color', 'k');

% 叠加模拟测量的 λ
valid_meas = Measured.lambda_fit_tr > 0 & isfinite(Measured.lambda_fit_tr);
if any(valid_meas)
    loglog(phi_data(valid_meas), Measured.lambda_fit_tr(valid_meas), ...
        'ko', 'MarkerSize', 12, 'LineWidth', 2, 'MarkerFaceColor', [0.3 0.8 0.3], ...
        'DisplayName', '\lambda_{measured} (T_{tr})');
end
valid_rot = Measured.lambda_fit_rot > 0 & isfinite(Measured.lambda_fit_rot);
if any(valid_rot)
    loglog(phi_data(valid_rot), Measured.lambda_fit_rot(valid_rot), ...
        'ms', 'MarkerSize', 12, 'LineWidth', 2, 'MarkerFaceColor', [0.8 0.3 0.8], ...
        'DisplayName', '\lambda_{measured} (T_{rot})');
end

% 标注临界 φ
xline(phi_c_kinetic, 'r-.', 'LineWidth', 1.5, 'HandleVisibility', 'off');
text(phi_c_kinetic*1.2, 0.5, sprintf('\\phi_c \\approx %.4f', phi_c_kinetic), ...
    'FontSize', 13, 'Color', 'r');

xlabel('\phi (Volume Fraction)', 'FontSize', 16);
ylabel('Penetration Depth \lambda [cm]', 'FontSize', 16);
legend('Location', 'best', 'FontSize', 13); legend boxoff;
set(gca, 'FontSize', 20, 'FontName', 'Times New Roman', 'LineWidth', 2);
ylim([0.001, 100]);

% ------ Fig 2: 壳层温度剖面 T(z) (z = dist_from_wall) ------
figure('Name', 'Fig.2: Temperature Profile T(z) from Wall', ...
    'Position', [150, 100, 850, 550]);

for i = 1:n_data
    T_prof = Measured.T_tr_profile(i, :);
    valid_s = T_prof > 0;
    
    % 数据点
    plot(dist_from_wall(valid_s), T_prof(valid_s), 'o', 'Color', cmap(i,:), ...
        'MarkerSize', 7, 'MarkerFaceColor', cmap(i,:), 'HandleVisibility', 'off');
    hold on;
    
    % 拟合曲线
    if isfinite(Measured.lambda_fit_tr(i)) && Measured.lambda_fit_tr(i) > 0
        z_dense = linspace(0, half_L, 200);
        T_fit = Measured.T_wall_fit(i) * exp(-z_dense / Measured.lambda_fit_tr(i));
        plot(z_dense, T_fit, '-', 'Color', cmap(i,:), 'LineWidth', 1.8, ...
            'DisplayName', sprintf('\\phi=%.3f, \\lambda=%.2f cm', phi_data(i), Measured.lambda_fit_tr(i)));
    else
        plot(dist_from_wall(valid_s), T_prof(valid_s), '-', 'Color', cmap(i,:), ...
            'LineWidth', 1.5, 'DisplayName', sprintf('\\phi=%.3f (uniform)', phi_data(i)));
    end
end

xlabel('Distance from Wall z [cm]', 'FontSize', 16);
ylabel('T_{tr} (Translational Temperature) [cm²/s²]', 'FontSize', 16);
legend('Location', 'best', 'FontSize', 13); legend boxoff;
set(gca, 'FontSize', 20, 'FontName', 'Times New Roman', 'LineWidth', 2);

% ------ Fig 3: 温度各向异性 Tx, Ty, Tz 壳层分布 ------
% 这个图是 3D 体系独有的: 壁面振动是各向同性的 (6面壁), 
% 所以理论上 Tx=Ty=Tz. 如果不等, 说明壁面驱动模式有各向异性.
figure('Name', 'Fig.3: Temperature Anisotropy by Shell (representative φ)', ...
    'Position', [200, 100, 800, 550]);

i_rep = ceil(n_data / 2);  % 取中间 phi 做代表
Tx = Measured.Tx_profile(i_rep, :);
Ty = Measured.Ty_profile(i_rep, :);
Tz = Measured.Tz_profile(i_rep, :);

plot(dist_from_wall, Tx, 'ro-', 'LineWidth', 1.5, 'MarkerSize', 6, ...
    'MarkerFaceColor', 'r', 'DisplayName', 'T_x');
hold on;
plot(dist_from_wall, Ty, 'b^-', 'LineWidth', 1.5, 'MarkerSize', 6, ...
    'MarkerFaceColor', 'b', 'DisplayName', 'T_y');
plot(dist_from_wall, Tz, 'gs-', 'LineWidth', 1.5, 'MarkerSize', 6, ...
    'MarkerFaceColor', 'g', 'DisplayName', 'T_z');

xlabel('Distance from Wall z [cm]', 'FontSize', 16);
ylabel('Temperature Component [cm²/s²]', 'FontSize', 16);
title(sprintf('Temperature Anisotropy (\\phi = %.3f)', phi_data(i_rep)), 'FontSize', 16);
legend('Location', 'best', 'FontSize', 14); legend boxoff;
set(gca, 'FontSize', 20, 'FontName', 'Times New Roman', 'LineWidth', 2);

% ------ Fig 4: T_tr 和 T_rot 同时画在壳层剖面上 ------
figure('Name', 'Fig.4: Translational vs Rotational T Profile', ...
    'Position', [250, 100, 800, 550]);

for i = 1:n_data
    T_tr = Measured.T_tr_profile(i, :);
    T_rot = Measured.T_rot_profile(i, :);
    
    valid_tr = T_tr > 0;
    valid_rot = T_rot > 0;
    
    h_tr = plot(dist_from_wall(valid_tr), T_tr(valid_tr), '-o', 'Color', cmap(i,:), ...
        'LineWidth', 1.5, 'MarkerSize', 6, 'MarkerFaceColor', cmap(i,:), ...
        'DisplayName', sprintf('T_{tr}, \\phi=%.3f', phi_data(i)));
    hold on;
    
    if any(valid_rot)
        plot(dist_from_wall(valid_rot), T_rot(valid_rot), '--s', 'Color', cmap(i,:), ...
            'LineWidth', 1.5, 'MarkerSize', 6, ...
            'DisplayName', sprintf('T_{rot}, \\phi=%.3f', phi_data(i)));
    end
end

xlabel('Distance from Wall z [cm]', 'FontSize', 16);
ylabel('Temperature [cm²/s²]', 'FontSize', 16);
legend('Location', 'best', 'FontSize', 11); legend boxoff;
set(gca, 'FontSize', 20, 'FontName', 'Times New Roman', 'LineWidth', 2);

% ------ Fig 5: 归一化温度剖面 T(z)/T(0) → 直接可见衰减 ------
figure('Name', 'Fig.5: Normalized Temperature Profile T(z)/T_wall', ...
    'Position', [300, 100, 800, 550]);

for i = 1:n_data
    T_prof = Measured.T_tr_profile(i, :);
    valid_s = T_prof > 0;
    
    if any(valid_s)
        % T_wall 取最外层 (dist_from_wall 最小处)
        T_wall_meas = T_prof(find(valid_s, 1, 'last'));  % 最外壳层
        T_norm = T_prof / T_wall_meas;
        
        plot(dist_from_wall(valid_s), T_norm(valid_s), '-o', 'Color', cmap(i,:), ...
            'LineWidth', 1.5, 'MarkerSize', 6, 'MarkerFaceColor', cmap(i,:), ...
            'DisplayName', sprintf('\\phi=%.3f', phi_data(i)));
        hold on;
    end
end

% 理论参考线
z_ref = linspace(0, half_L, 100);
for lam_ref = [0.5, 1.0, 2.0, 5.0]
    plot(z_ref, exp(-z_ref / lam_ref), 'k:', 'LineWidth', 1, 'HandleVisibility', 'off');
    text(half_L*0.75, exp(-half_L*0.75/lam_ref)*1.05, ...
        sprintf('\\lambda=%.1f', lam_ref), 'FontSize', 10, 'Color', [0.5 0.5 0.5]);
end

xlabel('Distance from Wall z [cm]', 'FontSize', 16);
ylabel('T(z) / T_{wall}', 'FontSize', 16);
legend('Location', 'best', 'FontSize', 13); legend boxoff;
set(gca, 'FontSize', 20, 'FontName', 'Times New Roman', 'LineWidth', 2);
ylim([0, 1.2]);

% ------ Fig 6: 相图: λ/L vs φ + 理论曲线 ------
figure('Name', 'Fig.6: Phase Diagram λ/L vs φ', ...
    'Position', [350, 100, 750, 550]);

% 理论线
semilogx(phi_theory, lambda_textbook / L, 'b-', 'LineWidth', 2, ...
    'DisplayName', 'Theory (textbook)');
hold on;
semilogx(phi_theory, lambda_kinetic / L, 'r--', 'LineWidth', 2, ...
    'DisplayName', 'Theory (Enskog)');

% 测量值
if any(valid_meas)
    semilogx(phi_data(valid_meas), Measured.lambda_fit_tr(valid_meas) / L, ...
        'ko', 'MarkerSize', 12, 'LineWidth', 2, 'MarkerFaceColor', [0.3 0.8 0.3], ...
        'DisplayName', 'Simulation');
end

% λ/L = 0.5 → 相变线
yline(0.5, 'k:', 'LineWidth', 2, 'HandleVisibility', 'off');
text(phi_theory(1)*3, 0.55, '\lambda = L/2 (Phase Transition)', 'FontSize', 13);

% 着色区域
fill([phi_theory(1), phi_theory(end), phi_theory(end), phi_theory(1)], ...
    [0, 0, 0.5, 0.5], 'r', 'FaceAlpha', 0.08, 'EdgeColor', 'none', ...
    'HandleVisibility', 'off');
text(phi_theory(end)*0.3, 0.2, 'Clustering Region', 'FontSize', 14, ...
    'Color', [0.8 0 0], 'FontWeight', 'bold');

fill([phi_theory(1), phi_theory(end), phi_theory(end), phi_theory(1)], ...
    [0.5, 0.5, 5, 5], 'b', 'FaceAlpha', 0.05, 'EdgeColor', 'none', ...
    'HandleVisibility', 'off');
text(phi_theory(1)*3, 2.5, 'Homogeneous Gas', 'FontSize', 14, ...
    'Color', [0 0 0.7], 'FontWeight', 'bold');

xlabel('\phi', 'FontSize', 16);
ylabel('\lambda / L', 'FontSize', 16);
legend('Location', 'northeast', 'FontSize', 13); legend boxoff;
set(gca, 'FontSize', 20, 'FontName', 'Times New Roman', 'LineWidth', 2);
ylim([0, 5]);

%% ========== 6. 教材公式 (5.106) 的数值解 (3D 版本) ==========
%
% 精确求解热方程: d/dr [K(T) dT/dr] = Γ(T, φ)
% 边界条件: T(r=L/2) = T0 (壁面温度), dT/dr|_{r=0} = 0 (对称)
% 这里 r 是 Chebyshev 距离

fprintf('\n=== 数值求解能量方程 (代表性 φ) ===\n');

% 对几个代表性 φ 求解
phi_solve = [0.001, 0.005, 0.01, 0.02, 0.05, 0.1];
T0_wall = 1.0;  % 壁面温度 (归一化)

figure('Name', 'Fig.7: Numerical Solution of Energy Equation', ...
    'Position', [400, 100, 800, 550]);
cmap_solve = jet(length(phi_solve));

for ip = 1:length(phi_solve)
    phi_s = phi_solve(ip);
    
    % 构建 ODE: d/dz [K dT/dz] = Γ
    % K = C_K * √T,  Γ = C_Γ * φ² T^{3/2}
    % → C_K √T T'' + C_K/(2√T) (T')² = C_Γ φ² T^{3/2}
    % 令 y1 = T, y2 = dT/dz
    % dy1/dz = y2
    % dy2/dz = [C_Γ φ² y1^{3/2} - C_K/(2√y1) y2²] / (C_K √y1)
    %        = (C_Γ/C_K) φ² y1 - y2²/(2 y1)
    
    % 比值 C_Γ/C_K 决定了衰减速率
    % C_Γ/C_K ∝ (1-e²)/d² (从教材中 K 和 Γ 的表达式)
    ratio = (1 - e^2) / d^2;
    
    odefun = @(z, y) [y(2); ratio * phi_s^2 * y(1) - y(2)^2 / (2*y(1))];
    
    % 边界条件: T(L/2) = T0, dT/dz(0) = 0
    % 用 shooting method: 猜 T(0) = T_center, 从中心积分到壁面
    % 调整 T_center 使得 T(L/2) = T0
    
    z_span = [0, half_L];
    
    % 二分法搜索 T_center
    T_lo = 0.001 * T0_wall;
    T_hi = T0_wall;
    
    for iter = 1:50
        T_guess = (T_lo + T_hi) / 2;
        try
            [z_sol, y_sol] = ode45(odefun, z_span, [T_guess, 0], ...
                odeset('RelTol', 1e-8, 'AbsTol', 1e-10));
            T_at_wall = y_sol(end, 1);
            
            if T_at_wall < T0_wall
                T_lo = T_guess;
            else
                T_hi = T_guess;
            end
            
            if abs(T_at_wall - T0_wall) / T0_wall < 1e-6
                break;
            end
        catch
            T_hi = T_guess;
        end
    end
    
    % 绘制 (翻转为 dist_from_wall)
    z_from_wall = half_L - z_sol;
    plot(z_from_wall, y_sol(:,1) / T0_wall, '-', 'Color', cmap_solve(ip,:), ...
        'LineWidth', 2, 'DisplayName', sprintf('\\phi = %.3f', phi_s));
    hold on;
    
    fprintf('  φ=%.4f: T_center/T_wall = %.4f\n', phi_s, T_guess/T0_wall);
end

xlabel('Distance from Wall z [cm]', 'FontSize', 16);
ylabel('T(z) / T_0', 'FontSize', 16);
legend('Location', 'best', 'FontSize', 13); legend boxoff;
set(gca, 'FontSize', 20, 'FontName', 'Times New Roman', 'LineWidth', 2);
ylim([0, 1.1]);

%% ========== 7. 汇总 ==========

fprintf('\n============================================================\n');
fprintf('  温度穿透深度 λ 汇总\n');
fprintf('============================================================\n');
fprintf('%8s | %12s | %12s | %12s | %12s | %8s\n', ...
    'phi', 'λ_theory', 'λ_meas(tr)', 'λ_meas(rot)', 'T_center/T_wall', 'λ<L/2?');
fprintf('%s\n', repmat('-', 1, 78));

for i = 1:n_data
    phi_val = phi_data(i);
    
    % 理论值
    chi_val = (1 - phi_val/2) / (1 - phi_val)^3;
    lam_th = d / (6 * phi_val * chi_val * sqrt(1-e^2));
    
    % 测量值
    lam_tr = Measured.lambda_fit_tr(i);
    lam_rot = Measured.lambda_fit_rot(i);
    
    % T_center/T_wall 比
    T_prof = Measured.T_tr_profile(i, :);
    valid_s = T_prof > 0;
    if any(valid_s)
        T_ratio = T_prof(find(valid_s, 1)) / T_prof(find(valid_s, 1, 'last'));
    else
        T_ratio = NaN;
    end
    
    cluster_flag = '';
    if lam_tr < half_L, cluster_flag = ' <<<'; end
    
    fprintf('%8.4f | %12.4f | %12.4f | %12.4f | %12.4f | %8s\n', ...
        phi_val, lam_th, lam_tr, lam_rot, T_ratio, ...
        [ternary(lam_tr < half_L, 'YES', 'no'), cluster_flag]);
end

fprintf('============================================================\n');
fprintf('  理论临界 φ_c ≈ %.6f (λ = L/2 = %.1f cm)\n', phi_c_kinetic, half_L);
fprintf('  当 λ < L/2 时, 中心区域冷却 → 密度增大 → 团簇成核\n');
fprintf('============================================================\n');

%% 辅助函数
function result = ternary(cond, a, b)
    if cond, result = a; else, result = b; end
end