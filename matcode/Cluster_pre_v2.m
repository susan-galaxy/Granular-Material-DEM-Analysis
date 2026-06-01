%% =====================================================================
%  转动-平动耦合修正的穿透深度与相变预测
%  Rotational Dissipation Corrected λ and Phase Transition
%
%  基于模拟实际参数:
%    d_bath = 2*Rb = 0.08 cm,  ρ_s = 1.08 g/cm³
%    e_n (粒-粒) ≈ 0.94 (Gn = -0.94)
%    e_n (粒-壁) ≈ 0.99 (Gnw = -0.99)
%    β = 0.1  (切向恢复系数, 非常粗糙!)
%    μ = 0.68 (粒-粒摩擦), μ_w = 0.2 (粒-壁摩擦)
%    振动: Ax=0.1, Ay=0.08, Az=0.07 cm; fx=5, fy=6, fz=7 Hz
%
%  理论框架:
%  1. Lun (1991), Herbst et al. (2000): 粗糙球碰撞的能量耗散
%     每次碰撞平动能量损失: ΔE_tr = (1-e²)/2 + κ(1+β)²/(1+κ)² · (1-...)
%     其中 κ = 2I/(md²) = 2/5 (均质球)
%  2. 稳态能量均分: T_rot/T_tr = f(e, β, κ)
%  3. 修正穿透深度: λ_corr = d / [6φχ √(ξ_eff)]
%     ξ_eff = (1-e²) + α_rot(1-β²)  [等效耗散系数]
% =====================================================================

clear all; close all; clc;

%% ========== 1. 模拟参数 (从 input-space.inp 提取) ==========

% --- 颗粒参数 ---
Rb    = 0.04;         % 浴颗粒半径 [cm]
d     = 2 * Rb;       % 浴颗粒直径 [cm] = 0.08 cm = 0.8 mm
rho_s = 1.08;         % 颗粒密度 [g/cm³]
m     = rho_s * (pi/6) * d^3;  % 颗粒质量 [g]
I_mom = (2/5) * m * (d/2)^2;   % 转动惯量 [g·cm²]
kappa = 2 * I_mom / (m * (d/2)^2);  % = 2/5 for uniform sphere

% --- 碰撞参数 ---
e_pp  = 0.94;         % 粒-粒法向恢复系数 (|Gn| = 0.94)
e_pw  = 0.99;         % 粒-壁法向恢复系数 (|Gnw| = 0.99)
beta  = 0.1;          % 切向恢复系数 (非常粗糙! β→-1为完全粗糙, β→1为光滑)
mu_pp = 0.5;         % 粒-粒摩擦系数
mu_pw = 0.2;          % 粒-壁摩擦系数

% --- 容器参数 ---
L     = 4.0;          % 容器边长 [cm]
half_L = L / 2;

% --- 驱动参数 (各向异性!) ---
Amp = [0.1, 0.08, 0.07];    % 振幅 [cm]
freq = [5, 6, 7];           % 频率 [Hz]
v_wall = Amp .* (2*pi*freq); % 壁面峰值速度 [cm/s]
E_wall = 0.5 * m * v_wall.^2; % 各方向壁面注入能量尺度

fprintf('=== 模拟参数汇总 ===\n');
fprintf('  颗粒: d=%.2f cm (%.1f mm), ρ=%.2f g/cm³, m=%.4e g\n', d, d*10, rho_s, m);
fprintf('  碰撞: e_pp=%.2f, e_pw=%.2f, β=%.2f, μ_pp=%.2f, μ_pw=%.2f\n', ...
    e_pp, e_pw, beta, mu_pp, mu_pw);
fprintf('  κ = 2I/(md²/4) = %.2f (均质球=0.4)\n', kappa);
fprintf('  驱动: v_wall = [%.2f, %.2f, %.2f] cm/s\n', v_wall);
fprintf('  驱动各向异性: Ax·fx/Az·fz = %.2f\n', (Amp(1)*freq(1))/(Amp(3)*freq(3)));

%% ========== 2. 粗糙球碰撞理论: T_rot/T_tr 预测 ==========
%
% 参考: Herbst, Cafiero, Zippelius et al., Phys. Fluids (2000)
%        Goldhirsch, Noskowicz, Bar-Lev, PRL (2005)
%        Brilliantov et al., "Kinetic Theory of Granular Gases" Ch.8
%
% 对于均质粗糙球 (κ = 2/5):
% 在均匀稳态下, 碰撞的角动量传递导致:
%
% T_rot/T_tr 的理论预测 (Goldhirsch 2005, Eq. for rough spheres):
%
% 定义 η = (1+β)/(1+κ)  (切向碰撞的有效传递系数, κ=2/5)
%   → η = (1+0.1)/(1+0.4) = 1.1/1.4 = 0.786
%
% 稳态 T_rot/T_tr 由详细平衡给出:
%   T_rot/T_tr = κ·η·(2-η) / (1 - η(2κ-κη)/(1+κ))
% 
% 但更常用的是 Herbst et al. 的近似:
%   在均匀加热系统中:
%   T_rot/T_tr ≈ κ(1+β)² / [(1+κ)² - κ(1+β)²]
%              = κ η² (1+κ)² / [(1+κ)² - κ(1+β)²]   ... 简化:

fprintf('\n=== 粗糙球碰撞理论 ===\n');

eta = (1 + beta) / (1 + kappa);
fprintf('  有效传递系数 η = (1+β)/(1+κ) = %.4f\n', eta);

% --- 方法1: Goldhirsch-Noskowicz-Bar-Lev 理论 ---
% 对均质球 (κ=2/5), 稳态条件:
% R_ratio = T_rot/T_tr
% 来自碰撞积分的平衡:
%   平动冷却率 ζ_tr ∝ (1-e²) + κ(1+β)²/(1+κ)² · [1 - T_rot/(κ T_tr)]
%   转动冷却率 ζ_rot ∝ (1+β)²/(1+κ)² · [T_rot/(κ T_tr) - 1]
% 在稳态 ζ_rot = 0 (无转动能量注入) 时:
%   T_rot/(κ T_tr) = 1  →  T_rot/T_tr = κ = 2/5

% 但在壁面驱动系统中, 壁面同时注入平动和转动能量!
% 如果壁面粗糙 (μ_w > 0), 碰壁时摩擦力产生力矩 → 注入转动能量
% 这就解释了为什么你的 T_rot >> T_tr

% --- 壁面注入的转动能量估算 ---
% 碰壁时, 如果颗粒-壁面有摩擦 μ_w = 0.2:
% 切向冲量 J_t = min(μ_w · J_n, -(1+β_w)·m_eff·v_t/(1+1/κ))
% 其中 J_n 是法向冲量, v_t 是切向相对速度
% 转动能量增量: ΔE_rot = J_t² · R² / (2I) = J_t² / (2·(2/5)·m)

% 简化估算: 壁面注入的转动能量 ∝ μ_w² · v_wall² · m
% 而平动能量 ∝ (1+e_pw) · v_wall² · m

% 壁面转动能量注入率 / 平动能量注入率:
%   ~ [μ_w · (1+beta_w)]² / [(1+e_pw)² · (1+κ)] 
% 这里的 beta_w 需要从壁面碰撞参数推断

% 更实际地, 直接从数据中提取 T_rot/T_tr 的实测值

%% ========== 3. 从数据提取 T_rot/T_tr(φ) ==========

pos_vel_dir = 'G:\Space_Active\Data\Maxwell_Boltz\data\Pos_Vel_data';
rot_dir     = 'G:\Space_Active\Data\Maxwell_Boltz\data\Rotation_data\Rotation_data';
pv_pattern  = 'Data_FixSpace_R0_04_mode0td_phi*.mat';
rot_pattern = 'RotData_FixSpace_R0_04_mode0td_phi*.mat';
tag_name    = 'Tag_2';
delta_step  = 5;

pv_files = dir(fullfile(pos_vel_dir, pv_pattern));
rot_files = dir(fullfile(rot_dir, rot_pattern));
extract_phi = @(fname) str2double(regexp(fname, '(?<=phi)\d+\.?\d*', 'match', 'once'));

pv_map = containers.Map('KeyType','double','ValueType','char');
rot_map = containers.Map('KeyType','double','ValueType','char');
for k = 1:length(pv_files)
    phi_k = extract_phi(pv_files(k).name);
    if ~isnan(phi_k), pv_map(phi_k) = pv_files(k).name; end
end
for k = 1:length(rot_files)
    phi_k = extract_phi(rot_files(k).name);
    if ~isnan(phi_k), rot_map(phi_k) = rot_files(k).name; end
end

phi_pv = cell2mat(pv_map.keys);
phi_rot = cell2mat(rot_map.keys);
phi_both = sort(intersect(phi_pv, phi_rot));
n_data = length(phi_both);

% 存储
Data = struct();
Data.phi = phi_both;
Data.T_tr_global = zeros(n_data, 1);
Data.T_rot_global = zeros(n_data, 1);
Data.T_tr_xyz = zeros(n_data, 3);      % [Tx, Ty, Tz]
Data.T_rot_xyz = zeros(n_data, 3);     % [ωx², ωy², ωz²]
Data.ratio_rot_tr = zeros(n_data, 1);  % T_rot / T_tr
Data.lambda_meas = zeros(n_data, 1);   % 从之前拟合得到 (或重新算)

% 壳层参数
n_shells = 20;
shell_edges = linspace(0, half_L, n_shells + 1);
shell_centers = (shell_edges(1:end-1) + shell_edges(2:end)) / 2;
dist_from_wall = half_L - shell_centers;

Data.T_tr_shell = zeros(n_data, n_shells);
Data.T_rot_shell = zeros(n_data, n_shells);

for i = 1:n_data
    phi_val = phi_both(i);
    fprintf('[%d/%d] φ=%.4f ... ', i, n_data, phi_val);
    
    % 平动
    loaded = load(fullfile(pos_vel_dir, pv_map(phi_val)));
    PD = loaded.ExtractedData.(tag_name);
    
    % 转动
    loaded_rot = load(fullfile(rot_dir, rot_map(phi_val)));
    PD_rot = loaded_rot.ExtractedData.(tag_name);
    rot_steps = PD_rot.TimeStep;
    
    total_frames = length(PD.Pos);
    valid_frames = 1:delta_step:total_frames;
    
    % 全局累加器
    sum_vx2 = 0; sum_vy2 = 0; sum_vz2 = 0; count_tr = 0;
    sum_ox2 = 0; sum_oy2 = 0; sum_oz2 = 0; count_rot = 0;
    
    % 壳层累加器
    sh_tr_sum = zeros(n_shells, 1); sh_tr_cnt = zeros(n_shells, 1);
    sh_rot_sum = zeros(n_shells, 1); sh_rot_cnt = zeros(n_shells, 1);
    
    for j_idx = 1:length(valid_frames)
        j = valid_frames(j_idx);
        
        pos_j = PD.Pos{j}; vel_j = PD.Vel{j};
        if isempty(pos_j), continue; end
        
        bad = any(isnan(pos_j),2) | any(isnan(vel_j),2);
        pos_j(bad,:) = []; vel_j(bad,:) = [];
        if isempty(pos_j), continue; end
        
        mean_v = mean(vel_j, 1);
        vf = vel_j - mean_v;
        N_j = size(pos_j, 1);
        
        sum_vx2 = sum_vx2 + sum(vf(:,1).^2);
        sum_vy2 = sum_vy2 + sum(vf(:,2).^2);
        sum_vz2 = sum_vz2 + sum(vf(:,3).^2);
        count_tr = count_tr + N_j;
        
        cheby = max(abs(pos_j), [], 2);
        s_idx = discretize(cheby, shell_edges);
        
        % 转动
        rot_idx = find(rot_steps == j, 1);
        omg_j = [];
        if ~isempty(rot_idx) && ~isempty(PD_rot.Omg{rot_idx})
            omg_j = PD_rot.Omg{rot_idx};
            if size(omg_j,1) == N_j + sum(bad)
                omg_j(bad,:) = [];
            end
        end
        
        has_omg = ~isempty(omg_j) && size(omg_j,1) == N_j;
        
        if has_omg
            sum_ox2 = sum_ox2 + sum(omg_j(:,1).^2);
            sum_oy2 = sum_oy2 + sum(omg_j(:,2).^2);
            sum_oz2 = sum_oz2 + sum(omg_j(:,3).^2);
            count_rot = count_rot + N_j;
        end
        
        % 壳层统计
        for s = 1:n_shells
            mask = (s_idx == s);
            if any(mask)
                sh_tr_sum(s) = sh_tr_sum(s) + sum(sum(vf(mask,:).^2, 2));
                sh_tr_cnt(s) = sh_tr_cnt(s) + sum(mask);
                
                if has_omg
                    sh_rot_sum(s) = sh_rot_sum(s) + sum(sum(omg_j(mask,:).^2, 2));
                    sh_rot_cnt(s) = sh_rot_cnt(s) + sum(mask);
                end
            end
        end
    end
    
    % 全局温度 (T = <v²> per DOF)
    Data.T_tr_xyz(i,:) = [sum_vx2, sum_vy2, sum_vz2] / count_tr;
    Data.T_tr_global(i) = (sum_vx2 + sum_vy2 + sum_vz2) / (3 * count_tr);
    
    if count_rot > 0
        Data.T_rot_xyz(i,:) = [sum_ox2, sum_oy2, sum_oz2] / count_rot;
        % 转动温度: T_rot = I <ω²> / 3  (每个转动自由度)
        % 但为了与 T_tr = m<v²>/3 可比, 需要 I/m 的换算
        % T_rot_energy = I·<ω²>/3, T_tr_energy = m·<v²>/3
        % 比值: T_rot/T_tr = (I<ω²>) / (m<v²>)
        T_rot_energy = I_mom * (sum_ox2 + sum_oy2 + sum_oz2) / (3 * count_rot);
        T_tr_energy = m * (sum_vx2 + sum_vy2 + sum_vz2) / (3 * count_tr);
        Data.T_rot_global(i) = T_rot_energy;
        Data.ratio_rot_tr(i) = T_rot_energy / T_tr_energy;
    end
    
    % 壳层温度
    for s = 1:n_shells
        if sh_tr_cnt(s) > 0
            Data.T_tr_shell(i,s) = m * sh_tr_sum(s) / (3 * sh_tr_cnt(s));
        end
        if sh_rot_cnt(s) > 0
            Data.T_rot_shell(i,s) = I_mom * sh_rot_sum(s) / (3 * sh_rot_cnt(s));
        end
    end
    
    fprintf('T_rot/T_tr = %.4f\n', Data.ratio_rot_tr(i));
end

%% ========== 4. 从 T_rot/T_tr 反推有效 β_eff ==========
%
% 理论 (Herbst et al. 2000, Goldhirsch 2005):
% 在自由冷却的均匀稳态中, 转动-平动温度比取决于 β:
%
% 对均质球 (κ = 2/5):
% 转动冷却率 ∝ (1+β)²/(1+κ)² · [T_rot/(κ·T_tr) - 1]
% 平动冷却率 ∝ (1-e²) + κ(1+β)²/(1+κ)² · [1 - T_rot/(κ·T_tr)]
%
% 稳态 (冷却率相等): T_rot/(κ·T_tr) = 1  →  T_rot/T_tr = κ = 2/5
%
% 但这是**自由冷却**的结果. 在壁面驱动系统中, 壁面额外注入转动能量,
% 破坏了这个均分. 如果壁面主要注入平动能量, T_rot/T_tr < κ;
% 如果壁面通过摩擦高效注入转动能量, T_rot/T_tr > κ.
%
% 定义有效切向恢复系数 β_eff:
% 在一个等效的自由加热系统中, 要产生相同的 T_rot/T_tr 比值,
% 需要的 β 值是多少?
%
% 更实用的方法: 直接从能量平衡反推
% 
% 壁面注入: 转动功率 P_rot_in, 平动功率 P_tr_in
% 碰撞耗散: 转动耗散率 Γ_rot, 平动耗散率 Γ_tr
% 稳态: P_rot_in = Γ_rot, P_tr_in = Γ_tr
%
% 碰撞耗散的解析形式 (Brilliantov Ch.8):
%   Γ_tr = n²σd² √(πT_tr/m) · [(1-e²)T_tr + κ(1+β)²/(1+κ)² · (T_tr - T_rot/κ)]
%   Γ_rot = n²σd² √(πT_tr/m) · [(1+β)²/(1+κ)² · (T_rot/κ - T_tr)]
%
% 从 Γ_rot = P_rot_in 和 Γ_tr = P_tr_in 的比值:
%   P_rot_in/P_tr_in = [(1+β)²/(1+κ)² · (T_rot/κ - T_tr)] / 
%                       [(1-e²)T_tr + κ(1+β)²/(1+κ)² · (T_tr - T_rot/κ)]

fprintf('\n=== 反推有效 β ===\n');
fprintf('  输入参数: e_pp=%.2f, κ=%.2f, β_input=%.2f\n', e_pp, kappa, beta);
fprintf('  理论均分预测: T_rot/T_tr = κ = %.2f (自由冷却极限)\n\n', kappa);

% 对每个 φ, 用测量的 T_rot/T_tr 反推 β_eff
% 
% 在粒-粒碰撞中:
% 稳态条件 (忽略壁面): Γ_rot = 0 要求 T_rot/(κ T_tr) = 1
% 如果 T_rot/(κ T_tr) > 1, 说明有额外的转动能量注入源 (壁面摩擦)
%
% 我们定义 "等效β_eff" 使得:
% 在均匀驱动系统中, T_rot/T_tr 的测量值与理论自洽:
%   R ≡ T_rot/T_tr
%   R = κ · [1 + P_rot/(P_tr · A)] 
% 其中 A = (1-e²)·(1+κ)² / [κ(1+β)²] + κ - 1
%
% 实际上, 更直接的做法: 把测量的 R 代入能量平衡方程, 解出 β_eff

Data.beta_eff = zeros(n_data, 1);
Data.xi_eff = zeros(n_data, 1);       % 等效总耗散系数
Data.xi_tr_only = zeros(n_data, 1);   % 纯平动耗散 (1-e²)
Data.xi_rot_only = zeros(n_data, 1);  % 转动耗散贡献

for i = 1:n_data
    R = Data.ratio_rot_tr(i);  % T_rot / T_tr (能量单位)
    
    if R <= 0 || isnan(R)
        Data.beta_eff(i) = beta;
        continue;
    end
    
    % 方法: 从碰撞微观力学出发
    % 每次粒-粒碰撞的总能量损失:
    %   ΔE_total = ΔE_normal + ΔE_tangential
    %   ΔE_normal = -(1-e²)/4 · m · (v_n)²  (法向)
    %   ΔE_tangential = -κ(1-β²)/[4(1+κ)] · m · (v_t)²  (切向)
    %
    % 平均后:
    %   <ΔE_normal> ∝ (1-e²) · T_tr
    %   <ΔE_tangential> ∝ κ(1-β²)/(1+κ) · (T_tr + T_rot/κ)
    %
    % 总耗散系数:
    %   ξ_eff = (1-e²) + κ(1-β²)/(1+κ) · (1 + R/κ)
    %         = (1-e²) + (1-β²)/(1+κ) · (κ + R)
    
    xi_normal = 1 - e_pp^2;
    xi_tang = (1 - beta^2) / (1 + kappa) * (kappa + R);
    
    Data.xi_tr_only(i) = xi_normal;
    Data.xi_rot_only(i) = xi_tang;
    Data.xi_eff(i) = xi_normal + xi_tang;
    
    % 反推 β_eff: 如果我们只用 (1-β_eff²) 来参数化总耗散
    % ξ_eff = (1-β_eff²)  → β_eff = √(1 - ξ_eff)
    if Data.xi_eff(i) < 1
        Data.beta_eff(i) = sqrt(1 - Data.xi_eff(i));
    else
        Data.beta_eff(i) = 0;  % 完全非弹性
    end
    
    fprintf('  φ=%.4f: R=T_rot/T_tr=%.4f, ξ_norm=%.4f, ξ_tang=%.4f, ξ_eff=%.4f\n', ...
        phi_both(i), R, xi_normal, xi_tang, Data.xi_eff(i));
end

%% ========== 5. 修正穿透深度 λ_corr ==========
%
% 原始教材公式: λ = d / [6φχ √(1-e²)]
% 修正公式:     λ_corr = d / [6φχ √(ξ_eff)]
%   ξ_eff = (1-e²) + (1-β²)/(1+κ) · (κ + T_rot/T_tr)

phi_theory = logspace(-4, -0.5, 500)';

% Carnahan-Starling
chi_theory = (1 - phi_theory/2) ./ (1 - phi_theory).^3;

% 原始 λ (只有法向耗散)
lambda_original = d ./ (6 * phi_theory .* chi_theory * sqrt(1 - e_pp^2));

% 修正 λ (加入转动耗散)
% 用所有 φ 的平均 ξ_eff (或者内插)
xi_eff_avg = mean(Data.xi_eff(Data.xi_eff > 0));
lambda_corrected = d ./ (6 * phi_theory .* chi_theory * sqrt(xi_eff_avg));

% φ-dependent 修正 (ξ_eff 可能随 φ 变化)
lambda_corr_phi = zeros(n_data, 1);
lambda_orig_phi = zeros(n_data, 1);

for i = 1:n_data
    phi_val = phi_both(i);
    chi_val = (1 - phi_val/2) / (1 - phi_val)^3;
    
    lambda_orig_phi(i) = d / (6 * phi_val * chi_val * sqrt(1 - e_pp^2));
    
    if Data.xi_eff(i) > 0
        lambda_corr_phi(i) = d / (6 * phi_val * chi_val * sqrt(Data.xi_eff(i)));
    end
end

% 临界 φ
[~, idx_c_orig] = min(abs(lambda_original - half_L));
phi_c_orig = phi_theory(idx_c_orig);
[~, idx_c_corr] = min(abs(lambda_corrected - half_L));
phi_c_corr = phi_theory(idx_c_corr);

fprintf('\n=== 穿透深度修正 ===\n');
fprintf('  平均等效耗散: ξ_eff = %.4f (vs 纯法向 1-e² = %.4f)\n', xi_eff_avg, 1-e_pp^2);
fprintf('  耗散增强因子: ξ_eff/(1-e²) = %.2f 倍\n', xi_eff_avg / (1-e_pp^2));
fprintf('  原始 φ_c (法向) = %.6f\n', phi_c_orig);
fprintf('  修正 φ_c (含转动) = %.6f\n', phi_c_corr);
fprintf('  φ_c 降低了 %.1f%%\n', (1 - phi_c_corr/phi_c_orig)*100);

%% ========== 6. 拟合壳层温度剖面获取 λ_measured ==========

Data.lambda_meas_tr = zeros(n_data, 1);
Data.lambda_meas_rot = zeros(n_data, 1);

for i = 1:n_data
    % 平动 λ
    T_prof = Data.T_tr_shell(i, :)';
    valid = T_prof > 0 & ~isnan(T_prof);
    if sum(valid) >= 4
        try
            ft = fittype('a*exp(-x/b)', 'independent', 'x');
            opts = fitoptions(ft);
            opts.StartPoint = [max(T_prof(valid)), L/4];
            opts.Lower = [0, 0.001]; opts.Upper = [100*max(T_prof(valid)), 50];
            [fr, ~] = fit(dist_from_wall(valid)', T_prof(valid), ft, opts);
            Data.lambda_meas_tr(i) = fr.b;
        catch
            p = polyfit(dist_from_wall(valid)', log(T_prof(valid)), 1);
            Data.lambda_meas_tr(i) = -1/p(1);
        end
    end
    
    % 转动 λ
    T_rot_prof = Data.T_rot_shell(i, :)';
    valid_r = T_rot_prof > 0 & ~isnan(T_rot_prof);
    if sum(valid_r) >= 4
        try
            ft = fittype('a*exp(-x/b)', 'independent', 'x');
            opts = fitoptions(ft);
            opts.StartPoint = [max(T_rot_prof(valid_r)), L/4];
            opts.Lower = [0, 0.001]; opts.Upper = [100*max(T_rot_prof(valid_r)), 50];
            [fr_rot, ~] = fit(dist_from_wall(valid_r)', T_rot_prof(valid_r), ft, opts);
            Data.lambda_meas_rot(i) = fr_rot.b;
        catch
            Data.lambda_meas_rot(i) = Inf;
        end
    end
end

%% ========== 7. 绘图 ==========

cmap = lines(n_data);

% ------ Fig 1: 修正 λ(φ) 相图 ------
figure('Name', 'Fig.1: Corrected λ vs φ', 'Position', [100, 100, 850, 600]);

loglog(phi_theory, lambda_original, 'b-', 'LineWidth', 2, ...
    'DisplayName', sprintf('\\lambda_{orig} (1-e^2=%.3f)', 1-e_pp^2));
hold on;
loglog(phi_theory, lambda_corrected, 'r-', 'LineWidth', 2.5, ...
    'DisplayName', sprintf('\\lambda_{corr} (\\xi_{eff}=%.3f)', xi_eff_avg));

% 测量值
valid_m = Data.lambda_meas_tr > 0 & isfinite(Data.lambda_meas_tr);
loglog(phi_both(valid_m), Data.lambda_meas_tr(valid_m), 'ko', ...
    'MarkerSize', 12, 'LineWidth', 2, 'MarkerFaceColor', [0.3 0.8 0.3], ...
    'DisplayName', '\lambda_{meas} (T_{tr})');

valid_r = Data.lambda_meas_rot > 0 & isfinite(Data.lambda_meas_rot);
if any(valid_r)
    loglog(phi_both(valid_r), Data.lambda_meas_rot(valid_r), 'ms', ...
        'MarkerSize', 12, 'LineWidth', 2, 'MarkerFaceColor', [0.8 0.3 0.8], ...
        'DisplayName', '\lambda_{meas} (T_{rot})');
end

yline(half_L, 'k:', 'LineWidth', 2, 'HandleVisibility', 'off');
text(0.3, half_L*1.3, 'L/2 = 2 cm', 'FontSize', 14);

xline(phi_c_orig, 'b-.', 'LineWidth', 1.5, 'HandleVisibility', 'off');
xline(phi_c_corr, 'r-.', 'LineWidth', 1.5, 'HandleVisibility', 'off');
text(phi_c_orig*1.3, 0.03, sprintf('\\phi_c^{orig}=%.4f', phi_c_orig), ...
    'FontSize', 12, 'Color', 'b');
text(phi_c_corr*0.3, 0.03, sprintf('\\phi_c^{corr}=%.4f', phi_c_corr), ...
    'FontSize', 12, 'Color', 'r');

xlabel('\phi', 'FontSize', 16);
ylabel('Penetration Depth \lambda [cm]', 'FontSize', 16);
legend('Location', 'southwest', 'FontSize', 13); legend boxoff;
set(gca, 'FontSize', 20, 'FontName', 'Times New Roman', 'LineWidth', 2);
ylim([0.01, 100]);

% ------ Fig 2: T_rot/T_tr vs φ ------
figure('Name', 'Fig.2: T_rot/T_tr vs φ', 'Position', [150, 100, 750, 550]);

valid_ratio = Data.ratio_rot_tr > 0;
semilogx(phi_both(valid_ratio), Data.ratio_rot_tr(valid_ratio), 'ko-', ...
    'LineWidth', 2, 'MarkerSize', 10, 'MarkerFaceColor', [0.3 0.3 0.3]);
hold on;

yline(kappa, 'r--', 'LineWidth', 2, 'Label', ...
    sprintf('Equipartition (\\kappa = 2/5 = %.1f)', kappa), ...
    'LabelHorizontalAlignment', 'left', 'FontSize', 14);
yline(1, 'b:', 'LineWidth', 1.5, 'Label', 'T_{rot} = T_{tr}', ...
    'LabelHorizontalAlignment', 'left', 'FontSize', 14);

xlabel('\phi', 'FontSize', 16);
ylabel('T_{rot} / T_{tr}  (energy units)', 'FontSize', 16);
set(gca, 'FontSize', 20, 'FontName', 'Times New Roman', 'LineWidth', 2);

% ------ Fig 3: 耗散分解 ξ_eff = ξ_normal + ξ_tangential ------
figure('Name', 'Fig.3: Dissipation Decomposition', 'Position', [200, 100, 750, 550]);

bar_data = [Data.xi_tr_only, Data.xi_rot_only];
b = bar(1:n_data, bar_data, 'stacked');
b(1).FaceColor = [0.2 0.4 0.8];
b(2).FaceColor = [0.9 0.3 0.2];

set(gca, 'XTick', 1:n_data, 'XTickLabel', ...
    arrayfun(@(x) sprintf('%.3f', x), phi_both, 'UniformOutput', false));
xtickangle(45);

xlabel('\phi', 'FontSize', 16);
ylabel('Effective Dissipation Coefficient \xi_{eff}', 'FontSize', 16);
legend({'Normal (1-e^2)', 'Tangential \propto (1-\beta^2)'}, ...
    'FontSize', 14, 'Location', 'northwest'); legend boxoff;
set(gca, 'FontSize', 16, 'FontName', 'Times New Roman', 'LineWidth', 2);

% 标注各个 bar 上的 ξ_eff 值
for i = 1:n_data
    text(i, Data.xi_eff(i) + 0.01, sprintf('%.3f', Data.xi_eff(i)), ...
        'HorizontalAlignment', 'center', 'FontSize', 9);
end

% ------ Fig 4: 修正后的数值解 vs 实测温度剖面 ------
figure('Name', 'Fig.4: Corrected ODE Solution vs Data', 'Position', [250, 100, 850, 600]);

% 选几个代表性 φ 对比
phi_show = [0.01, 0.02, 0.03, 0.05];
colors_show = {'b', 'r', [0 0.6 0], 'm'};

for ip = 1:length(phi_show)
    [~, i_data] = min(abs(phi_both - phi_show(ip)));
    phi_val = phi_both(i_data);
    xi_val = Data.xi_eff(i_data);
    
    % 数据点
    T_prof = Data.T_tr_shell(i_data, :);
    valid = T_prof > 0;
    T_wall_meas = T_prof(find(valid, 1, 'last'));
    
    plot(dist_from_wall(valid), T_prof(valid)/T_wall_meas, 'o', ...
        'Color', colors_show{ip}, 'MarkerSize', 8, ...
        'MarkerFaceColor', colors_show{ip}, 'HandleVisibility', 'off');
    hold on;
    
    % 数值解 (修正后)
    ratio_ode = xi_val / d^2;  % C_Γ/C_K 比值, 用 ξ_eff 替换 (1-e²)
    odefun = @(z, y) [y(2); ratio_ode * phi_val^2 * y(1) - y(2)^2/(2*y(1))];
    
    T_lo = 0.001; T_hi = 1.0;
    for iter = 1:60
        T_guess = (T_lo + T_hi) / 2;
        try
            [z_sol, y_sol] = ode45(odefun, [0, half_L], [T_guess, 0], ...
                odeset('RelTol',1e-8,'AbsTol',1e-12));
            T_at_wall = y_sol(end, 1);
            if T_at_wall < 1.0, T_lo = T_guess; else, T_hi = T_guess; end
            if abs(T_at_wall - 1.0) < 1e-6, break; end
        catch
            T_hi = T_guess;
        end
    end
    
    z_from_wall = half_L - z_sol;
    plot(z_from_wall, y_sol(:,1), '-', 'Color', colors_show{ip}, 'LineWidth', 2, ...
        'DisplayName', sprintf('\\phi=%.3f, \\xi_{eff}=%.3f', phi_val, xi_val));
    
    % 也画原始的 (只有法向)
    ratio_orig = (1-e_pp^2) / d^2;
    odefun_orig = @(z, y) [y(2); ratio_orig * phi_val^2 * y(1) - y(2)^2/(2*y(1))];
    T_lo2 = 0.001; T_hi2 = 1.0;
    for iter = 1:60
        T_g2 = (T_lo2 + T_hi2)/2;
        try
            [z2, y2] = ode45(odefun_orig, [0, half_L], [T_g2, 0], ...
                odeset('RelTol',1e-8,'AbsTol',1e-12));
            if y2(end,1) < 1, T_lo2 = T_g2; else, T_hi2 = T_g2; end
            if abs(y2(end,1)-1) < 1e-6, break; end
        catch
            T_hi2 = T_g2;
        end
    end
    plot(half_L - z2, y2(:,1), '--', 'Color', colors_show{ip}, 'LineWidth', 1.5, ...
        'HandleVisibility', 'off');
end

xlabel('Distance from Wall z [cm]', 'FontSize', 16);
ylabel('T(z) / T_{wall}', 'FontSize', 16);
legend('Location', 'best', 'FontSize', 13); legend boxoff;
set(gca, 'FontSize', 20, 'FontName', 'Times New Roman', 'LineWidth', 2);
ylim([0, 1.1]);

% 手动添加注释
text(1.5, 0.95, 'Solid: corrected (\xi_{eff})', 'FontSize', 13);
text(1.5, 0.88, 'Dashed: original (1-e^2 only)', 'FontSize', 13, 'Color', [0.5 0.5 0.5]);

% ------ Fig 5: λ/L 相图 (修正版) ------
figure('Name', 'Fig.5: Corrected Phase Diagram', 'Position', [300, 100, 800, 600]);

semilogx(phi_theory, lambda_original/L, 'b-', 'LineWidth', 2, ...
    'DisplayName', 'Original (1-e^2 only)');
hold on;
semilogx(phi_theory, lambda_corrected/L, 'r-', 'LineWidth', 2.5, ...
    'DisplayName', 'Corrected (+rotation)');

if any(valid_m)
    semilogx(phi_both(valid_m), Data.lambda_meas_tr(valid_m)/L, 'ko', ...
        'MarkerSize', 12, 'LineWidth', 2, 'MarkerFaceColor', [0.3 0.8 0.3], ...
        'DisplayName', 'Simulation');
end

yline(0.5, 'k:', 'LineWidth', 2, 'HandleVisibility', 'off');
text(1e-4*3, 0.55, '\lambda = L/2 (Clustering Onset)', 'FontSize', 13);

% 着色
fill([phi_theory(1), phi_theory(end), phi_theory(end), phi_theory(1)], ...
    [0, 0, 0.5, 0.5], 'r', 'FaceAlpha', 0.08, 'EdgeColor', 'none', 'HandleVisibility', 'off');
fill([phi_theory(1), phi_theory(end), phi_theory(end), phi_theory(1)], ...
    [0.5, 0.5, 5, 5], 'b', 'FaceAlpha', 0.05, 'EdgeColor', 'none', 'HandleVisibility', 'off');

text(0.2, 0.15, 'Clustering', 'FontSize', 15, 'Color', [0.8 0 0], 'FontWeight', 'bold');
text(3e-4, 3, 'Homogeneous Gas', 'FontSize', 15, 'Color', [0 0 0.7], 'FontWeight', 'bold');

xlabel('\phi', 'FontSize', 16);
ylabel('\lambda / L', 'FontSize', 16);
legend('Location', 'northeast', 'FontSize', 13); legend boxoff;
set(gca, 'FontSize', 20, 'FontName', 'Times New Roman', 'LineWidth', 2);
ylim([0, 5]);

% ------ Fig 6: 壳层 T_rot/T_tr 剖面 (空间分辨的能量均分) ------
figure('Name', 'Fig.6: Shell-resolved T_rot/T_tr', 'Position', [350, 100, 800, 550]);

for i = 1:n_data
    ratio_shell = zeros(1, n_shells);
    for s = 1:n_shells
        if Data.T_tr_shell(i,s) > 0 && Data.T_rot_shell(i,s) > 0
            ratio_shell(s) = Data.T_rot_shell(i,s) / Data.T_tr_shell(i,s);
        end
    end
    valid_s = ratio_shell > 0;
    if any(valid_s)
        plot(dist_from_wall(valid_s), ratio_shell(valid_s), '-o', ...
            'Color', cmap(i,:), 'LineWidth', 1.5, 'MarkerSize', 5, ...
            'MarkerFaceColor', cmap(i,:), ...
            'DisplayName', sprintf('\\phi=%.3f', phi_both(i)));
        hold on;
    end
end

yline(kappa, 'r--', 'LineWidth', 2, 'HandleVisibility', 'off');
text(0.1, kappa*1.1, '\kappa = 2/5 (equipartition)', 'FontSize', 13, 'Color', 'r');

xlabel('Distance from Wall z [cm]', 'FontSize', 16);
ylabel('T_{rot}(z) / T_{tr}(z)', 'FontSize', 16);
legend('Location', 'best', 'FontSize', 11); legend boxoff;
set(gca, 'FontSize', 20, 'FontName', 'Times New Roman', 'LineWidth', 2);

%% ========== 8. 汇总表 ==========

fprintf('\n======================================================================\n');
fprintf('  完整汇总: 转动修正的穿透深度与相变预测\n');
fprintf('======================================================================\n');
fprintf('  模拟参数: e=%.2f, β=%.2f, μ=%.2f, κ=%.2f\n', e_pp, beta, mu_pp, kappa);
fprintf('----------------------------------------------------------------------\n');
fprintf('%6s | %8s | %8s | %8s | %8s | %8s | %8s\n', ...
    'phi', 'T_rot/T_tr', 'xi_norm', 'xi_tang', 'xi_eff', 'lam_orig', 'lam_corr');
fprintf('%s\n', repmat('-', 1, 70));

for i = 1:n_data
    chi_val = (1 - phi_both(i)/2) / (1 - phi_both(i))^3;
    lam_o = d / (6*phi_both(i)*chi_val*sqrt(1-e_pp^2));
    lam_c = d / (6*phi_both(i)*chi_val*sqrt(max(Data.xi_eff(i), 0.001)));
    
    fprintf('%6.4f | %8.4f | %8.4f | %8.4f | %8.4f | %8.3f | %8.3f\n', ...
        phi_both(i), Data.ratio_rot_tr(i), Data.xi_tr_only(i), ...
        Data.xi_rot_only(i), Data.xi_eff(i), lam_o, lam_c);
end

fprintf('======================================================================\n');
fprintf('  关键结论:\n');
fprintf('  1. 转动耗散使等效 ξ 增大 %.1f 倍 → λ 缩短 %.1f 倍\n', ...
    xi_eff_avg/(1-e_pp^2), sqrt(xi_eff_avg/(1-e_pp^2)));
fprintf('  2. φ_c 从 %.5f 降至 %.5f (降低 %.0f%%)\n', ...
    phi_c_orig, phi_c_corr, (1-phi_c_corr/phi_c_orig)*100);
fprintf('  3. β=%.2f (粗糙) 是转动过热 T_rot>>T_tr 的根本原因\n', beta);
fprintf('======================================================================\n');

%% 保存
save_path = fullfile(pos_vel_dir, 'Rotation_Corrected_Lambda_Results.mat');
save(save_path, 'Data', 'phi_both', 'lambda_original', 'lambda_corrected', ...
    'phi_theory', 'phi_c_orig', 'phi_c_corr', 'xi_eff_avg', '-v7.3');
fprintf('结果已保存: %s\n', save_path);