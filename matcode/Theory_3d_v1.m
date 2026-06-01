%% =========================================================================
%  三维受限颗粒气体相变：DEM 数据提取 + 3D 修正理论 + 拟合对比
%  (纯脚本版；数据提取部分严格沿用 Cluster_pre_v3.m 的流程)
%
%  数据来源：
%    - pos_vel_dir/Data_Fix*.mat     逐帧位置 + 速度（PD.Pos{j}, PD.Vel{j}）
%    - rot_dir/RotData_...phi*.mat    逐帧角速度（PD_rot.Omg{ri}）
%
%  输出：
%    R 结构体（壳层分析结果，同 Cluster_pre_v3 风格）
%    + Tc_over_Tw_DEM, sigma2_n_DEM, kappa_DEM
%    + Tc_theory, sigma2_theory, u_profiles（理论扫描结果）
%    + 论文候选图 6 张
%
%  物理参数、核心 ODE 见《完整理论推导_最终版.md》
% =========================================================================

clear all; close all; clc;

%% ========== 1. 系统参数（理论 + 模拟共用）==================================
d         = 0.08;             % 颗粒直径 [cm]
L         = 4.0;              % 容器边长 [cm]
half_L    = L/2;
rho_s     = 1.08;             % 颗粒密度 [g/cm³]
m         = rho_s*(pi/6)*d^3;
I_mom     = (2/5)*m*(d/2)^2;  % 实心球转动惯量

e_n       = 0.94;             % 法向恢复系数
beta      = 0.1;              % 切向恢复系数
q2        = 2/5;              % 4I/(md²) for uniform solid sphere
gamma0    = 12/sqrt(pi);
kappa0_th = sqrt(pi)/6;

% 旋转温度比（会用 DEM 实测值自动更新）
kappa_rot = 0.16;
xi_smooth = 1 - e_n^2;
xi_rot    = (8/7)*(1-beta^2)*(1 - kappa_rot/q2);
xi_eff    = xi_smooth + xi_rot;

% 源强度（拟合主参数）—— 在本脚本末尾可以调整或做 grid search
s0        = 0.8;

fprintf('====== 理论耗散参数 ======\n');
fprintf('  ξ_smooth = %.4f\n', xi_smooth);
fprintf('  ξ_rot    = %.4f\n', xi_rot);
fprintf('  ξ_eff    = %.4f\n', xi_eff);
fprintf('  δ_rot    = %.2f （旋转耗散放大因子）\n', xi_rot/xi_smooth);


%% ========== 2. 数据配置（与 Cluster_pre_v3 完全一致）=======================
pos_vel_dir = 'G:\Space_Active\Data\Maxwell_Boltz\data\Pos_Vel_data';
rot_dir     = 'G:\Space_Active\Data\Maxwell_Boltz\data\Rotation_data\Rotation_data';
pv_pattern  = 'Data_Fix*.mat';
rot_pattern = 'RotData_FixSpace_R0_04_mode0td_phi*.mat';
tag_name    = 'Tag_2';
delta_step  = 5;
n_shells    = 20;


%% ========== 3. 文件检索 ==================================================
pv_files  = dir(fullfile(pos_vel_dir, pv_pattern));
rot_files = dir(fullfile(rot_dir, rot_pattern));
extract_phi = @(fname) str2double(regexp(fname, '(?<=phi)\d+\.?\d*', 'match', 'once'));

pv_map  = containers.Map('KeyType','double','ValueType','char');
rot_map = containers.Map('KeyType','double','ValueType','char');
for k=1:length(pv_files)
    p_=extract_phi(pv_files(k).name);
    if ~isnan(p_), pv_map(p_)=pv_files(k).name; end
end
for k=1:length(rot_files)
    p_=extract_phi(rot_files(k).name);
    if ~isnan(p_), rot_map(p_)=rot_files(k).name; end
end

phi_data = sort(cell2mat(pv_map.keys));
n_data   = length(phi_data);
fprintf('\n找到 %d 个体积分数\n', n_data);


%% ========== 4. 壳层分析（与 Cluster_pre_v3 相同的流程）====================
shell_edges    = linspace(0, half_L, n_shells+1);
shell_centers  = (shell_edges(1:end-1)+shell_edges(2:end))/2;
dist_from_wall = half_L - shell_centers;          % shell 1 = 中心，shell end = 墙
x_DEM_half     = dist_from_wall / L;              % ∈ [near 0.5, near 0]

R = struct();
R.phi            = phi_data;
R.T_tr_shell     = zeros(n_data, n_shells);
R.T_rot_shell    = zeros(n_data, n_shells);
R.Tx_shell       = zeros(n_data, n_shells);
R.Ty_shell       = zeros(n_data, n_shells);
R.Tz_shell       = zeros(n_data, n_shells);
R.n_shell        = zeros(n_data, n_shells);
R.T_tr_global    = zeros(n_data, 1);
R.T_rot_global   = zeros(n_data, 1);
R.ratio_rot_tr   = zeros(n_data, 1);
R.lambda_meas    = zeros(n_data, 1);
R.T_wall_fit     = zeros(n_data, 1);
R.sigma2_n       = zeros(n_data, 1);

for i = 1:n_data
    phi_val = phi_data(i);
    fprintf('[%d/%d] φ=%.4f\n', i, n_data, phi_val);

    loaded = load(fullfile(pos_vel_dir, pv_map(phi_val)));
    PD = loaded.ExtractedData.(tag_name);

    has_rot = rot_map.isKey(phi_val);
    if has_rot
        lr = load(fullfile(rot_dir, rot_map(phi_val)));
        PD_rot = lr.ExtractedData.(tag_name);
        rot_steps = PD_rot.TimeStep;
    end

    valid_frames = 1:delta_step:length(PD.Pos);

    sh_vx2  = zeros(n_shells,1); sh_vy2  = zeros(n_shells,1); sh_vz2 = zeros(n_shells,1);
    sh_omg2 = zeros(n_shells,1); sh_cnt  = zeros(n_shells,1); sh_rcnt = zeros(n_shells,1);
    g_vx2=0; g_vy2=0; g_vz2=0; g_ox2=0; g_oy2=0; g_oz2=0;
    g_cnt_tr=0; g_cnt_rot=0;

    frame_shell_counts = zeros(length(valid_frames), n_shells);

    for j_idx = 1:length(valid_frames)
        j = valid_frames(j_idx);
        pos_j = PD.Pos{j}; vel_j = PD.Vel{j};
        if isempty(pos_j), continue; end
        bad = any(isnan(pos_j),2)|any(isnan(vel_j),2);
        pos_j(bad,:)=[]; vel_j(bad,:)=[];
        if isempty(pos_j), continue; end

        N_j = size(pos_j,1);
        mean_v = mean(vel_j,1);
        vf = vel_j - mean_v;                       % 去除质心漂移
        cheby = max(abs(pos_j),[],2);              % Chebyshev 距离（到中心）
        s_idx = discretize(cheby, shell_edges);

        g_vx2=g_vx2+sum(vf(:,1).^2);
        g_vy2=g_vy2+sum(vf(:,2).^2);
        g_vz2=g_vz2+sum(vf(:,3).^2);
        g_cnt_tr=g_cnt_tr+N_j;

        omg_j = [];
        if has_rot
            ri = find(rot_steps==j,1);
            if ~isempty(ri) && ~isempty(PD_rot.Omg{ri})
                omg_j = PD_rot.Omg{ri};
                if size(omg_j,1)==N_j+sum(bad), omg_j(bad,:)=[]; end
            end
        end
        has_omg = ~isempty(omg_j) && size(omg_j,1)==N_j;
        if has_omg
            g_ox2=g_ox2+sum(omg_j(:,1).^2);
            g_oy2=g_oy2+sum(omg_j(:,2).^2);
            g_oz2=g_oz2+sum(omg_j(:,3).^2);
            g_cnt_rot=g_cnt_rot+N_j;
        end

        for s=1:n_shells
            mask = (s_idx==s);
            if any(mask)
                sh_vx2(s) = sh_vx2(s) + sum(vf(mask,1).^2);
                sh_vy2(s) = sh_vy2(s) + sum(vf(mask,2).^2);
                sh_vz2(s) = sh_vz2(s) + sum(vf(mask,3).^2);
                sh_cnt(s) = sh_cnt(s) + sum(mask);
                frame_shell_counts(j_idx,s) = sum(mask);
                if has_omg
                    sh_omg2(s) = sh_omg2(s) + sum(sum(omg_j(mask,:).^2,2));
                    sh_rcnt(s) = sh_rcnt(s) + sum(mask);
                end
            end
        end
    end

    % 壳层温度
    for s=1:n_shells
        if sh_cnt(s)>0
            R.Tx_shell(i,s)    = sh_vx2(s)/sh_cnt(s);
            R.Ty_shell(i,s)    = sh_vy2(s)/sh_cnt(s);
            R.Tz_shell(i,s)    = sh_vz2(s)/sh_cnt(s);
            R.T_tr_shell(i,s)  = (sh_vx2(s)+sh_vy2(s)+sh_vz2(s))/(3*sh_cnt(s));
            R.n_shell(i,s)     = sh_cnt(s)/length(valid_frames);
        end
        if sh_rcnt(s)>0
            R.T_rot_shell(i,s) = I_mom*sh_omg2(s)/(3*sh_rcnt(s));
        end
    end

    % 全局温度
    R.T_tr_global(i) = (g_vx2+g_vy2+g_vz2)/(3*g_cnt_tr);
    if g_cnt_rot>0
        R.T_rot_global(i) = I_mom*(g_ox2+g_oy2+g_oz2)/(3*g_cnt_rot);
        R.ratio_rot_tr(i) = R.T_rot_global(i)/(m*R.T_tr_global(i));
    end

    % 密度涨落
    total_per_frame = sum(frame_shell_counts, 2);
    valid_f = total_per_frame > 0;
    if any(valid_f)
        n_frames_valid = total_per_frame(valid_f);
        R.sigma2_n(i) = var(n_frames_valid) / mean(n_frames_valid);
    end

    % 拟合 λ：T(z) = T_wall · exp(-z/λ),  z = dist_from_wall
    T_prof = R.T_tr_shell(i,:)';
    valid_s = T_prof>0 & ~isnan(T_prof);
    if sum(valid_s)>=4
        try
            ft = fittype('a*exp(-x/b)','independent','x');
            opts = fitoptions(ft);
            opts.StartPoint = [max(T_prof(valid_s)), L/4];
            opts.Lower = [0, 0.001]; opts.Upper = [100*max(T_prof(valid_s)), 50];
            [fr,~] = fit(dist_from_wall(valid_s)', T_prof(valid_s), ft, opts);
            R.lambda_meas(i) = fr.b;
            R.T_wall_fit(i)  = fr.a;
        catch
            p_fit = polyfit(dist_from_wall(valid_s)', log(T_prof(valid_s)), 1);
            R.lambda_meas(i) = -1/p_fit(1);
            R.T_wall_fit(i)  = exp(p_fit(2));
        end
    end
end


%% ========== 5. 从 DEM 导出的关键量 =======================================
% T_c/T_w：中心壳层温度 / 指数拟合的 T_wall
Tc_over_Tw_DEM = nan(n_data, 1);
for i = 1:n_data
    if R.T_wall_fit(i) > 0 && R.T_tr_shell(i,1) > 0 && isfinite(R.T_wall_fit(i))
        Tc_over_Tw_DEM(i) = R.T_tr_shell(i,1) / R.T_wall_fit(i);
    end
end

% DEM 测得的 κ = T_rot/T_tr
valid_k    = R.ratio_rot_tr > 0 & isfinite(R.ratio_rot_tr);
kappa_DEM  = mean(R.ratio_rot_tr(valid_k));
kappa_DEM_std = std(R.ratio_rot_tr(valid_k));

fprintf('\n====== DEM 衍生量 ======\n');
fprintf('  <T_rot/T_tr> = %.4f ± %.4f （理论默认 0.16）\n', kappa_DEM, kappa_DEM_std);
fprintf('  T_c/T_w 范围: [%.3f, %.3f]\n', ...
        min(Tc_over_Tw_DEM(~isnan(Tc_over_Tw_DEM))), ...
        max(Tc_over_Tw_DEM(~isnan(Tc_over_Tw_DEM))));

% 将 DEM 实测 κ 回代给理论（如果显著偏离 0.16）
if ~isnan(kappa_DEM)
    kappa_rot = kappa_DEM;
    xi_rot    = (8/7)*(1-beta^2)*(1 - kappa_rot/q2);
    xi_eff    = xi_smooth + xi_rot;
    fprintf('  已更新理论 κ = %.4f  →  ξ_eff = %.4f\n', kappa_rot, xi_eff);
end


%% ========== 6. 理论扫描：自洽打靶法求解 u''=C·[1/u - s0·g(x)] ===============
% 生成扫描网格（加密相变区 + 包含 DEM 的每一个 φ）
phi_th_list = sort(unique([linspace(max(1e-3, 0.7*min(phi_data)), 0.012, 10), ...
                           linspace(0.012, 0.025, 22), ...
                           linspace(0.025, min(0.08, 1.3*max(phi_data)), 12), ...
                           phi_data(:)']));
n_th = length(phi_th_list);

Tc_theory       = nan(n_th, 1);
sigma2_theory   = nan(n_th, 1);
nwall_theory    = nan(n_th, 1);
lambda_b_theory = zeros(n_th, 1);
A_coef_theory   = zeros(n_th, 1);
Lambda_ctrl     = nan(n_th, 1);
converged_th    = false(n_th, 1);
x_profiles      = cell(n_th, 1);
u_profiles      = cell(n_th, 1);

ode_opts = odeset('RelTol',1e-6,'AbsTol',1e-9,'MaxStep',0.01);
uc_prev  = 0.99;

fprintf('\n====== 理论扫描 (s0=%.2f, e_n=%.3f) ======\n', s0, e_n);
for ii = 1:n_th
    phi = phi_th_list(ii);
    chi      = (1-phi/2)/(1-phi)^3;
    lambda_b = d / (6*sqrt(2)*phi*chi);
    Kn       = lambda_b / L;
    nmean    = phi * 6/(pi*d^3);
    n_w      = nmean;
    lambda_b_theory(ii) = lambda_b;

    converged = false;
    uc_root   = NaN;
    u_full    = [];  x_full = [];  u2 = NaN;

    for sc_it = 1:20   % 自洽 n_w 迭代
        A    = L^2 * gamma0 * n_w * d / (2*kappa0_th);
        coef = A * chi * xi_eff;

        % --- 打靶第1步：粗扫 uc 找 u(0)=1 的括号 ---
        uc_grid = [0.005, 0.02, 0.05:0.05:0.99];
        u0_vals = nan(size(uc_grid));
        for kg = 1:length(uc_grid)
            try
                [~, Y] = ode45(@(x,y) [y(2); ...
                    coef*(1/max(y(1),1e-4) - s0*(exp(-x/Kn)+exp(-(1-x)/Kn)+4*exp(-0.5/Kn)))], ...
                    [0.5, 0.0], [uc_grid(kg); 0], ode_opts);
                if all(Y(:,1) > 0) && all(isfinite(Y(:,1)))
                    u0_vals(kg) = Y(end,1);
                end
            catch, end
        end

        % 找所有根的括号
        bracket_pairs = [];
        for kg = 1:length(uc_grid)-1
            if isnan(u0_vals(kg)) || isnan(u0_vals(kg+1)), continue; end
            if (u0_vals(kg)-1)*(u0_vals(kg+1)-1) <= 0
                bracket_pairs(end+1,:) = [uc_grid(kg), uc_grid(kg+1)]; %#ok<SAGROW>
            end
        end
        if isempty(bracket_pairs), break; end  % 无解 → 相变

        % 选离 uc_prev 最近的括号（分支延续）
        mid_pts = mean(bracket_pairs, 2);
        [~, idx_best] = min(abs(mid_pts - uc_prev));
        uc_lo = bracket_pairs(idx_best, 1);
        uc_hi = bracket_pairs(idx_best, 2);

        % --- 打靶第2步：二分法细化 ---
        for bi = 1:25
            uc_mid = 0.5*(uc_lo + uc_hi);
            try
                [~, Y] = ode45(@(x,y) [y(2); ...
                    coef*(1/max(y(1),1e-4) - s0*(exp(-x/Kn)+exp(-(1-x)/Kn)+4*exp(-0.5/Kn)))], ...
                    [0.5, 0.0], [uc_mid; 0], ode_opts);
                if ~all(Y(:,1) > 0) || ~all(isfinite(Y(:,1)))
                    break;
                end
                u0_mid = Y(end, 1);
            catch, break; end
            if u0_mid > 1, uc_hi = uc_mid; else, uc_lo = uc_mid; end
            if abs(u0_mid - 1) < 1e-5, break; end
        end
        uc_root = 0.5*(uc_lo + uc_hi);

        % --- 用 uc_root 积分完整剖面 ---
        try
            [x_out, Y_out] = ode45(@(x,y) [y(2); ...
                coef*(1/max(y(1),1e-4) - s0*(exp(-x/Kn)+exp(-(1-x)/Kn)+4*exp(-0.5/Kn)))], ...
                linspace(0.5, 0.0, 51), [uc_root; 0], ode_opts);
        catch, break; end

        x_half = fliplr(x_out');
        u_half = fliplr(Y_out(:,1)');
        if any(u_half < 1e-3) || any(~isfinite(u_half)), break; end

        % 对称展开
        x_full = [x_half, 1 - fliplr(x_half(1:end-1))];
        u_full = [u_half, fliplr(u_half(1:end-1))];

        % 自洽更新 n_w
        u2        = trapz(x_full, u_full.^(-2));
        n_w_new   = nmean / u2;
        rel_err   = abs(n_w_new - n_w) / n_w;
        n_w       = 0.5*n_w + 0.5*n_w_new;
        uc_prev   = uc_root;

        if rel_err < 1e-4, converged = true; break; end
    end

    converged_th(ii)  = converged;
    nwall_theory(ii)  = n_w;
    A_coef_theory(ii) = L^2 * gamma0 * n_w * d / (2*kappa0_th);
    g_center          = 2*exp(-0.5/Kn) + 4*exp(-0.5/Kn);
    Lambda_ctrl(ii)   = A_coef_theory(ii)*chi*xi_eff / (s0*g_center);

    if converged
        Tc_theory(ii)     = uc_root^2;
        x_profiles{ii}    = x_full;
        u_profiles{ii}    = u_full;
        u4                = trapz(x_full, u_full.^(-4));
        sigma2_theory(ii) = n_w * (u4 - u2^2) / u2;
    else
        Tc_theory(ii) = 0;
        uc_prev = 0.2;   % 相变后重置初值
    end

    if mod(ii,5)==0 || ii==1 || ii==n_th
        fprintf('  φ=%.4f  Kn=%.3f  A=%6.1f  Λ=%6.2f  T_c/T_w=%.4f  conv=%d\n', ...
                phi, Kn, A_coef_theory(ii), Lambda_ctrl(ii), Tc_theory(ii), converged);
    end
end

% 理论相变点
idx_trans = find(~converged_th, 1, 'first');
if ~isempty(idx_trans)
    phi_c_theory = phi_th_list(idx_trans);
    fprintf('>>> 理论预测相变点 φ_c ≈ %.4f\n', phi_c_theory);
else
    phi_c_theory = NaN;
    fprintf('>>> 扫描范围内未出现相变（减小 s0 或扩大 φ 上限）\n');
end


%% ========== 7. 论文级绘图（沿用 Cluster_pre_v3 风格）======================
set(0, 'DefaultAxesFontName', 'Times New Roman');
set(0, 'DefaultAxesFontSize', 18);
set(0, 'DefaultLineLineWidth', 1.5);

% 颜色表（蓝→红）
n_c = n_data;
cmap_custom = zeros(n_c, 3);
for ci = 1:n_c
    t = (ci-1)/(n_c-1);
    cmap_custom(ci,:) = (1-t)*[0.1 0.2 0.7] + t*[0.8 0.1 0.15];
end

% ===== Fig. 1: 主图 T_c/T_w + σ²_n 对 φ （双轴） =====
fig1 = figure('Name','Fig1_MainPhaseTransition','Position',[50,80,900,560]);

yyaxis left;
vt = ~isnan(Tc_theory);
plot(phi_th_list(vt), Tc_theory(vt), '-', 'Color',[0.1 0.2 0.7], 'LineWidth',2.5); hold on;
valid_TcD = ~isnan(Tc_over_Tw_DEM);
plot(phi_data(valid_TcD), Tc_over_Tw_DEM(valid_TcD), 'o', ...
     'MarkerSize',10, 'LineWidth',1.8, ...
     'MarkerEdgeColor',[0.05 0.1 0.4], 'MarkerFaceColor',[0.3 0.5 0.9]);
ylabel('$T_\mathrm{center}/T_\mathrm{wall}$','Interpreter','latex','Color',[0.1 0.2 0.7]);
ylim([0 1.15]); set(gca, 'YColor', [0.1 0.2 0.7]);

yyaxis right;
v2 = ~isnan(sigma2_theory) & sigma2_theory>0;
semilogy(phi_th_list(v2), sigma2_theory(v2), '-', 'Color',[0.8 0.1 0.15], 'LineWidth',2.5); hold on;
vd = R.sigma2_n > 0 & isfinite(R.sigma2_n);
semilogy(phi_data(vd), R.sigma2_n(vd), 's', ...
         'MarkerSize',10, 'LineWidth',1.8, ...
         'MarkerEdgeColor',[0.5 0.05 0.1], 'MarkerFaceColor',[0.9 0.3 0.3]);
ylabel('$\sigma_n^2/\langle n\rangle$','Interpreter','latex','Color',[0.8 0.1 0.15]);
set(gca, 'YColor', [0.8 0.1 0.15]);
set(gca,'YScale','log')

% 相变点标记
if ~isnan(phi_c_theory)
    xline(phi_c_theory, 'k--', 'LineWidth', 1.5, ...
          'Label', sprintf('$\\phi_c=%.4f$', phi_c_theory), ...
          'Interpreter','latex','LabelVerticalAlignment','top','FontSize',14);
end
xlabel('$\phi$', 'Interpreter','latex');
title(sprintf('$s_0=%.2f,\\ e_n=%.3f,\\ \\xi_{\\rm eff}=%.3f$', s0, e_n, xi_eff), ...
      'Interpreter','latex','FontSize',14);
 box on;

% ===== Fig. 2: λ 比较（理论 vs DEM 指数拟合） =====
fig2 = figure('Name','Fig2_Lambda','Position',[80,100,800,560]);

lambda_noir_theory = d./(6*phi_th_list.*(1-phi_th_list/2)./((1-phi_th_list).^3)*sqrt(xi_smooth));
loglog(phi_th_list, lambda_b_theory, '-',  'Color',[0.2 0.2 0.7], 'LineWidth',2.5); hold on;
loglog(phi_th_list, lambda_noir_theory, '--', 'Color',[0.0 0.5 0.7], 'LineWidth',2);
valid_L = R.lambda_meas > 0 & isfinite(R.lambda_meas);
loglog(phi_data(valid_L), R.lambda_meas(valid_L), 'o', ...
       'MarkerSize',11, 'LineWidth',1.8, ...
       'MarkerEdgeColor','k', 'MarkerFaceColor',[0.3 0.7 0.3]);
yline(half_L, ':', 'Color',[0.3 0.3 0.3], 'LineWidth',2);
text(1.2*min(phi_th_list), half_L*1.15, '$L/2$', 'Interpreter','latex','FontSize',18,'Color',[0.3 0.3 0.3]);

xlabel('$\phi$', 'Interpreter','latex');
ylabel('$\lambda$ [cm]', 'Interpreter','latex');
legend({'$\lambda_b = d/(6\sqrt{2}\phi\chi)$ （本工作）', ...
        '$\lambda_{\rm Noir} = d/(6\phi\chi\sqrt{1-e^2})$', ...
        '$\lambda_{\rm meas}$ (DEM 指数拟合)'}, ...
       'Interpreter','latex', 'Location','southwest', 'FontSize',14);
legend boxoff;
set(gca, 'LineWidth',2, 'TickDir','in'); grid on; box on;

% ===== Fig. 3: 温度剖面对比 T(x)/T_wall =====
fig3 = figure('Name','Fig3_ProfileComparison','Position',[120,120,900,600]);

% 选 4–6 个代表性 φ
N_sel = min(6, n_data);
idx_sel = round(linspace(1, n_data, N_sel));
col_prof = jet(N_sel);

for k = 1:N_sel
    iD = idx_sel(k);
    phi = phi_data(iD);

    % DEM 剖面
    T_D = R.T_tr_shell(iD, :);
    Tw  = R.T_wall_fit(iD);
    if Tw <= 0 || ~isfinite(Tw), continue; end
    T_ratio_D = T_D / Tw;
    vd_k      = T_D > 0 & isfinite(T_D);

    % 对应的理论剖面
    [~, iT] = min(abs(phi_th_list - phi));
    if converged_th(iT)
        xT = x_profiles{iT};
        uT = u_profiles{iT};
        mask_half = xT <= 0.5;
        plot(xT(mask_half), uT(mask_half).^2, '-', ...
             'Color', col_prof(k,:), 'LineWidth',2.5, ...
             'DisplayName', sprintf('$\\phi=%.4f$ 理论', phi)); hold on;
    end
    plot(x_DEM_half(vd_k), T_ratio_D(vd_k), 'o', ...
         'Color', col_prof(k,:), 'MarkerSize',9, ...
         'MarkerFaceColor','w', 'LineWidth',1.8, ...
         'DisplayName', sprintf('$\\phi=%.4f$ DEM', phi));
end
xlabel('$x = (\mathrm{dist\,from\,wall})/L$', 'Interpreter','latex');
ylabel('$T(x)/T_\mathrm{wall}$', 'Interpreter','latex');
legend('Interpreter','latex', 'Location','eastoutside', 'FontSize',12); legend boxoff;
set(gca, 'LineWidth',2, 'TickDir','in');
xlim([0 0.5]); ylim([0 1.15]); grid on; box on;

% ===== Fig. 4: κ = T_rot/T_tr 的 φ 独立性（核心假设验证） =====
fig4 = figure('Name','Fig4_KappaValidation','Position',[160,140,750,500]);

semilogx(phi_data(valid_k), R.ratio_rot_tr(valid_k), 'o', ...
         'MarkerSize',11, 'LineWidth',1.8, ...
         'MarkerEdgeColor','k', 'MarkerFaceColor',[0.4 0.4 0.4]); hold on;
yline(kappa_DEM, '-', 'Color',[0.8 0.1 0.15], 'LineWidth',2.5);
yline(kappa_DEM+kappa_DEM_std, ':', 'Color',[0.8 0.1 0.15], 'LineWidth',1.2);
yline(kappa_DEM-kappa_DEM_std, ':', 'Color',[0.8 0.1 0.15], 'LineWidth',1.2);

text(min(phi_data(valid_k))*1.3, kappa_DEM+0.015, ...
     sprintf('$\\kappa = %.3f \\pm %.3f$', kappa_DEM, kappa_DEM_std), ...
     'Interpreter','latex', 'Color',[0.8 0.1 0.15], 'FontSize',16);

% χ²/ndf 独立性检验
chi2_kappa = sum( ((R.ratio_rot_tr(valid_k) - kappa_DEM)/kappa_DEM_std).^2 );
ndf_kappa  = sum(valid_k) - 1;
text(min(phi_data(valid_k))*1.3, kappa_DEM-0.05, ...
     sprintf('$\\chi^2/{\\rm ndf} = %.2f$', chi2_kappa/max(ndf_kappa,1)), ...
     'Interpreter','latex', 'FontSize',14);

xlabel('$\phi$', 'Interpreter','latex');
ylabel('$T_\mathrm{rot}/T_\mathrm{tr}$', 'Interpreter','latex');
title('$\kappa$ 与 $\phi$ 的独立性（模型闭合核心前提）', 'Interpreter','latex','FontSize',14);
set(gca, 'LineWidth',2, 'TickDir','in');
ylim([0, max(0.5, 1.8*max(R.ratio_rot_tr(valid_k)))]); grid on; box on;

% ===== Fig. 5: 密度涨落 σ²/⟨n⟩ 单独放大 =====
fig5 = figure('Name','Fig5_DensityFluctuation','Position',[200,160,700,500]);

semilogy(phi_th_list(v2), sigma2_theory(v2), '-', ...
         'Color',[0.1 0.2 0.7], 'LineWidth',2.5); hold on;
semilogy(phi_data(vd), R.sigma2_n(vd), 'o', ...
         'MarkerSize',11, 'LineWidth',1.8, ...
         'MarkerEdgeColor','k', 'MarkerFaceColor',[0.3 0.5 0.9]);
yline(1, '--', 'Color',[0.5 0.5 0.5], 'LineWidth',1.2);
text(0.9*max(phi_data), 1.4, 'Ideal Gas ($\sigma^2/\langle n\rangle=1$)', ...
     'Interpreter','latex', 'FontSize',14, 'Color',[0.5 0.5 0.5]);
if ~isnan(phi_c_theory)
    xline(phi_c_theory, 'k--', 'LineWidth',1.5);
end
xlabel('$\phi$', 'Interpreter','latex');
ylabel('$\sigma_n^2/\langle n\rangle$', 'Interpreter','latex');
legend({'理论', 'DEM'}, 'Interpreter','latex', 'Location','best'); legend boxoff;
set(gca, 'LineWidth',2, 'TickDir','in'); grid on; box on;

% ===== Fig. 6: 诊断图（λ_b/L, A, 3D vs 1D, Λ） =====
fig6 = figure('Name','Fig6_Diagnostics','Position',[240,180,1200,750]);

subplot(2,2,1);
semilogy(phi_th_list, lambda_b_theory/L, 'k-', 'LineWidth',2); grid on;
xlabel('$\phi$','Interpreter','latex'); ylabel('Kn $=\lambda_b/L$','Interpreter','latex');
title('(a) Knudsen 数','Interpreter','latex','FontSize',14);

subplot(2,2,2);
plot(phi_th_list, A_coef_theory, 'k-', 'LineWidth',2); grid on;
xlabel('$\phi$','Interpreter','latex'); ylabel('$A(\phi)$','Interpreter','latex');
title('(b) 无量纲耗散系数','Interpreter','latex','FontSize',14);

subplot(2,2,3);
g3D_ = zeros(size(phi_th_list));
g1D_ = zeros(size(phi_th_list));
for k = 1:length(phi_th_list)
    Kn_k = lambda_b_theory(k)/L;
    g3D_(k) = 2*exp(-0.5/Kn_k) + 4*exp(-0.5/Kn_k);
    g1D_(k) = 2*exp(-0.5/Kn_k);
end
semilogy(phi_th_list, g3D_, 'r-', phi_th_list, g1D_, 'b--', 'LineWidth',2); grid on;
xlabel('$\phi$','Interpreter','latex'); ylabel('$g(x{=}1/2;\phi)$','Interpreter','latex');
legend({'3D 六面墙','1D 二面墙'},'Interpreter','latex','Location','best','FontSize',12);
legend boxoff;
title('(c) 中心源强度：3D 托底','Interpreter','latex','FontSize',14);

subplot(2,2,4);
vL = ~isnan(Lambda_ctrl);
semilogy(phi_th_list(vL), Lambda_ctrl(vL), 'k-', 'LineWidth',2); grid on;
xlabel('$\phi$','Interpreter','latex'); ylabel('$\Lambda(\phi)$','Interpreter','latex');
title('(d) 相变控制参数','Interpreter','latex','FontSize',14);


%% ========== 8. 汇总表 ====================================================
fprintf('\n============================================================\n');
fprintf('  φ        T_c/T_w_DEM  T_c/T_w_thy  λ_meas    κ=T_rot/T_tr  σ²/⟨n⟩\n');
fprintf('------------------------------------------------------------\n');
for i = 1:n_data
    [~, iT] = min(abs(phi_th_list - phi_data(i)));
    fprintf('%7.4f   %9.4f    %9.4f    %7.3f     %8.4f    %8.3f\n', ...
            phi_data(i), Tc_over_Tw_DEM(i), Tc_theory(iT), ...
            R.lambda_meas(i), R.ratio_rot_tr(i), R.sigma2_n(i));
end
fprintf('============================================================\n');
fprintf('  理论 φ_c     = %.4f\n', phi_c_theory);
fprintf('  <T_rot/T_tr> = %.4f ± %.4f  (χ²/ndf = %.2f)\n', ...
        kappa_DEM, kappa_DEM_std, chi2_kappa/max(ndf_kappa,1));
fprintf('  ξ_eff / ξ_smooth = %.2f  （相变点相对光滑球理论的偏移因子）\n', xi_eff/xi_smooth);
fprintf('============================================================\n');


%% ========== 9. 保存结果 ==================================================
save_name = fullfile(pos_vel_dir, 'GranularGas_3D_Result.mat');
save(save_name, 'R', 'phi_data', 'dist_from_wall', 'shell_centers', ...
     'Tc_over_Tw_DEM', 'kappa_DEM', 'kappa_DEM_std', ...
     'phi_th_list', 'Tc_theory', 'sigma2_theory', ...
     'lambda_b_theory', 'A_coef_theory', 'Lambda_ctrl', ...
     'x_profiles', 'u_profiles', 'converged_th', ...
     'phi_c_theory', 's0', 'e_n', 'xi_eff', 'kappa_rot', '-v7.3');
fprintf('结果已保存至：%s\n', save_name);


%% ========== 10. 手动调参提示 ============================================
% 若理论与 DEM 不重合，按以下方式手动调节：
%
%  - 相变太晚（理论 φ_c 偏大）：减小 s0（例如 0.5, 0.3, 0.2）
%  - 相变太早（理论 φ_c 偏小）：增大 s0（例如 1.2, 1.5, 2.0）
%  - 低 φ 区域 T_c/T_w 偏高：轻微减小 e_n（0.92, 0.90），使 ξ_smooth 增加
%  - 低 φ 区域 T_c/T_w 偏低：增大 e_n（0.96, 0.98）
%
% 修改脚本开头的 s0 和 e_n 后重跑即可。若需自动拟合，可把 s0-loop 写成：
%
%   best_resid = Inf;  best_s0 = s0;
%   for s0_try = linspace(0.3, 2.0, 10)
%       s0 = s0_try;
%       [运行 Section 6 的核心代码...]
%       resid_tmp = sum((Tc_theory(id_DEM)-Tc_over_Tw_DEM).^2);
%       if resid_tmp < best_resid
%           best_resid = resid_tmp;  best_s0 = s0_try;
%       end
%   end