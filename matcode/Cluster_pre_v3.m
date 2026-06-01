%% =====================================================================
%  三维颗粒气体能量穿透深度与气液相变预测 (最终整合版)
%  Energy Penetration Depth λ and Gas-Liquid Transition in 3D Granular Gas
%
%  理论: Granular Media §5.4.4, Noirhomme et al. PRL 126, 128002 (2021)
%  λ ≃ d / [6φ χ(φ) √(1-e²)]
%  相变条件: λ < L/2
%
%  模拟参数: d=0.08cm, e=0.94, β=0.1, μ=0.68, L=4cm
% =====================================================================
clear all; close all; clc;

%% ========== 1. 系统参数 ==========
d     = 0.08;    % 颗粒直径 [cm]
e_n   = 1;    % 法向恢复系数
beta  = 1;     % 切向恢复系数
L     = 4.0;     % 容器边长 [cm]
half_L = L/2;
rho_s = 1.08;    % 密度 [g/cm³]
m     = rho_s*(pi/6)*d^3;
I_mom = (2/5)*m*(d/2)^2;
kappa = 0.4;     % 2/5 for uniform sphere

%% ========== 2. 数据配置 ==========
pos_vel_dir = 'G:\Space_Active\Data\Maxwell_Boltz\data\Pos_Vel_data\Nofric';
rot_dir     = 'G:\Space_Active\Data\Maxwell_Boltz\data\Rotation_data\Nofric';
pv_pattern  = 'Data_*.mat';
rot_pattern = 'RotData_*.mat';
tag_name    = 'Tag_2';
delta_step  = 5;
n_shells    = 20;

%% ========== 3. 文件检索 ==========
pv_files = dir(fullfile(pos_vel_dir, pv_pattern));
rot_files = dir(fullfile(rot_dir, rot_pattern));
extract_phi = @(fname) str2double(regexp(fname, '(?<=phi)\d+\.?\d*', 'match', 'once'));

pv_map = containers.Map('KeyType','double','ValueType','char');
rot_map = containers.Map('KeyType','double','ValueType','char');
for k=1:length(pv_files), p=extract_phi(pv_files(k).name); if ~isnan(p), pv_map(p)=pv_files(k).name; end; end
for k=1:length(rot_files), p=extract_phi(rot_files(k).name); if ~isnan(p), rot_map(p)=rot_files(k).name; end; end

phi_data = sort(cell2mat(pv_map.keys));
n_data = length(phi_data);
fprintf('找到 %d 个体积分数\n', n_data);

%% ========== 4. 壳层分析 ==========
shell_edges = linspace(0, half_L, n_shells+1);
shell_centers = (shell_edges(1:end-1)+shell_edges(2:end))/2;
dist_from_wall = half_L - shell_centers;

R = struct();
R.phi = phi_data;
R.T_tr_shell = zeros(n_data, n_shells);
R.T_rot_shell = zeros(n_data, n_shells);
R.Tx_shell = zeros(n_data, n_shells);
R.Ty_shell = zeros(n_data, n_shells);
R.Tz_shell = zeros(n_data, n_shells);
R.n_shell = zeros(n_data, n_shells);
R.T_tr_global = zeros(n_data,1);
R.T_rot_global = zeros(n_data,1);
R.ratio_rot_tr = zeros(n_data,1);
R.lambda_meas = zeros(n_data,1);
R.T_wall_fit = zeros(n_data,1);
R.sigma2_n = zeros(n_data,1);  % 密度涨落

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
    
    sh_vx2=zeros(n_shells,1); sh_vy2=zeros(n_shells,1); sh_vz2=zeros(n_shells,1);
    sh_omg2=zeros(n_shells,1); sh_cnt=zeros(n_shells,1); sh_rcnt=zeros(n_shells,1);
    g_vx2=0; g_vy2=0; g_vz2=0; g_ox2=0; g_oy2=0; g_oz2=0;
    g_cnt_tr=0; g_cnt_rot=0;
    
    % 密度涨落: 逐帧统计每个壳层的颗粒数
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
        vf = vel_j - mean_v;
        cheby = max(abs(pos_j),[],2);
        s_idx = discretize(cheby, shell_edges);
        
        g_vx2=g_vx2+sum(vf(:,1).^2); g_vy2=g_vy2+sum(vf(:,2).^2);
        g_vz2=g_vz2+sum(vf(:,3).^2); g_cnt_tr=g_cnt_tr+N_j;
        
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
            g_ox2=g_ox2+sum(omg_j(:,1).^2); g_oy2=g_oy2+sum(omg_j(:,2).^2);
            g_oz2=g_oz2+sum(omg_j(:,3).^2); g_cnt_rot=g_cnt_rot+N_j;
        end
        
        for s=1:n_shells
            mask=(s_idx==s);
            if any(mask)
                sh_vx2(s)=sh_vx2(s)+sum(vf(mask,1).^2);
                sh_vy2(s)=sh_vy2(s)+sum(vf(mask,2).^2);
                sh_vz2(s)=sh_vz2(s)+sum(vf(mask,3).^2);
                sh_cnt(s)=sh_cnt(s)+sum(mask);
                frame_shell_counts(j_idx,s)=sum(mask);
                if has_omg
                    sh_omg2(s)=sh_omg2(s)+sum(sum(omg_j(mask,:).^2,2));
                    sh_rcnt(s)=sh_rcnt(s)+sum(mask);
                end
            end
        end
    end
    
    % 壳层温度
    for s=1:n_shells
        if sh_cnt(s)>0
            R.Tx_shell(i,s) = sh_vx2(s)/sh_cnt(s);
            R.Ty_shell(i,s) = sh_vy2(s)/sh_cnt(s);
            R.Tz_shell(i,s) = sh_vz2(s)/sh_cnt(s);
            R.T_tr_shell(i,s) = (sh_vx2(s)+sh_vy2(s)+sh_vz2(s))/(3*sh_cnt(s));
            R.n_shell(i,s) = sh_cnt(s)/length(valid_frames);
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
    
    % 密度涨落 σ²/⟨n⟩
    total_per_frame = sum(frame_shell_counts, 2);
    valid_f = total_per_frame > 0;
    if any(valid_f)
        n_frames_valid = total_per_frame(valid_f);
        R.sigma2_n(i) = var(n_frames_valid) / mean(n_frames_valid);
    end
    
    % 拟合 λ
    T_prof = R.T_tr_shell(i,:)';
    valid_s = T_prof>0 & ~isnan(T_prof);
    if sum(valid_s)>=4
        try
            ft = fittype('a*exp(-x/b)','independent','x');
            opts = fitoptions(ft);
            opts.StartPoint = [max(T_prof(valid_s)), L/4];
            opts.Lower = [0,0.001]; opts.Upper = [100*max(T_prof(valid_s)),50];
            [fr,~] = fit(dist_from_wall(valid_s)', T_prof(valid_s), ft, opts);
            R.lambda_meas(i) = fr.b;
            R.T_wall_fit(i) = fr.a;
        catch
            p = polyfit(dist_from_wall(valid_s)', log(T_prof(valid_s)), 1);
            R.lambda_meas(i) = -1/p(1);
            R.T_wall_fit(i) = exp(p(2));
        end
    end
end

%% ========== 5. 理论曲线 ==========
phi_th = logspace(-4, -0.3, 500)';
chi_th = (1-phi_th/2)./(1-phi_th).^3;
lambda_th = d./(6*phi_th.*chi_th*sqrt(1-e_n^2));
[~,ic] = min(abs(lambda_th - half_L));
phi_c = phi_th(ic);

%% ========== 6. 出版级绘图 ==========
% 统一风格设置
set(0, 'DefaultAxesFontName', 'Times New Roman');
set(0, 'DefaultAxesFontSize', 18);
set(0, 'DefaultLineLineWidth', 1.5);

% 自定义色表 (16色, 从蓝到红)
n_c = n_data;
cmap_custom = zeros(n_c, 3);
for ci = 1:n_c
    t = (ci-1)/(n_c-1);
    cmap_custom(ci,:) = (1-t)*[0.1 0.2 0.7] + t*[0.8 0.1 0.15];
end

% ===== Figure 1: λ(φ) 相图 (核心图) =====
fig1 = figure('Name','Fig1_Lambda_vs_phi','Position',[50,100,700,550]);

loglog(phi_th, lambda_th, '-', 'Color', [0.2 0.2 0.7], 'LineWidth', 2.5); hold on;
valid_m = R.lambda_meas>0 & isfinite(R.lambda_meas);
loglog(phi_data(valid_m), R.lambda_meas(valid_m), 'o', ...
    'MarkerSize', 10, 'LineWidth', 1.8, ...
    'MarkerEdgeColor', [0.1 0.1 0.1], 'MarkerFaceColor', [0.3 0.7 0.3]);
yline(half_L, ':', 'Color', [0.3 0.3 0.3], 'LineWidth', 2);
text(0.25, half_L*1.35, '$L/2$', 'Interpreter', 'latex', 'FontSize', 20, 'Color', [0.3 0.3 0.3]);
xline(phi_c, '-.', 'Color', [0.7 0.2 0.2], 'LineWidth', 1.5);
text(phi_c*1.5, 0.06, sprintf('$\\phi_c \\approx %.3f$', phi_c), ...
    'Interpreter', 'latex', 'FontSize', 18, 'Color', [0.7 0.2 0.2]);

xlabel('$\phi$', 'Interpreter', 'latex', 'FontSize', 22);
ylabel('$\lambda$ [cm]', 'Interpreter', 'latex', 'FontSize', 22);
legend({'Theory: $\lambda = d/[6\phi\chi\sqrt{1-e^2}]$', 'DEM Simulation'}, ...
    'Interpreter', 'latex', 'FontSize', 16, 'Location', 'southwest'); legend boxoff;
set(gca, 'LineWidth', 2, 'TickDir', 'in', 'TickLength', [0.02 0.02]);
ylim([0.01 100]);
box on;

% ===== Figure 2: T_tr(z) 壳层温度剖面 =====
fig2 = figure('Name','Fig2_Shell_T_profile','Position',[100,100,750,550]);

for i = 1:n_data
    T_prof = R.T_tr_shell(i,:);
    vs = T_prof > 0;
    plot(dist_from_wall(vs), T_prof(vs), 'o', 'Color', cmap_custom(i,:), ...
        'MarkerSize', 5, 'MarkerFaceColor', cmap_custom(i,:), 'HandleVisibility','off');
    hold on;
    
    if R.lambda_meas(i)>0 && isfinite(R.lambda_meas(i))
        zd = linspace(0, half_L, 150);
        Tfit = R.T_wall_fit(i)*exp(-zd/R.lambda_meas(i));
        plot(zd, Tfit, '-', 'Color', cmap_custom(i,:), 'LineWidth', 1.3, ...
            'DisplayName', sprintf('$\\phi$=%.3f, $\\lambda$=%.2f', phi_data(i), R.lambda_meas(i)));
    end
end
xlabel('Distance from wall $z$ [cm]', 'Interpreter', 'latex', 'FontSize', 22);
ylabel('$T_\mathrm{tr}$ [cm$^2$/s$^2$]', 'Interpreter', 'latex', 'FontSize', 22);
legend('Interpreter','latex','FontSize',11,'Location','northeast','NumColumns',2); legend boxoff;
set(gca, 'LineWidth', 2, 'TickDir', 'in');
box on;

% ===== Figure 3: 归一化温度 T(z)/T_wall =====
fig3 = figure('Name','Fig3_Normalized_T','Position',[150,100,700,550]);

for i = 1:n_data
    T_prof = R.T_tr_shell(i,:);
    vs = T_prof>0;
    if any(vs)
        Tw = T_prof(find(vs,1,'last'));
        plot(dist_from_wall(vs), T_prof(vs)/Tw, '-o', 'Color', cmap_custom(i,:), ...
            'MarkerSize', 4, 'MarkerFaceColor', cmap_custom(i,:), 'LineWidth', 1.2, ...
            'DisplayName', sprintf('%.3f', phi_data(i)));
        hold on;
    end
end
% 参考线
zr = linspace(0, half_L, 100);
for lv = [0.5, 1.0, 2.0, 5.0]
    plot(zr, exp(-zr/lv), 'k:', 'LineWidth', 0.8, 'HandleVisibility', 'off');
    text(1.7, exp(-1.7/lv)+0.03, sprintf('$\\lambda$=%.1f', lv), ...
        'Interpreter','latex','FontSize',11,'Color',[0.5 0.5 0.5]);
end
xlabel('$z$ [cm]', 'Interpreter', 'latex', 'FontSize', 22);
ylabel('$T(z)/T_\mathrm{wall}$', 'Interpreter', 'latex', 'FontSize', 22);
lgd = legend('Location','southwest','FontSize',10,'NumColumns',4);
title(lgd, '$\phi$', 'Interpreter', 'latex'); legend boxoff;
set(gca, 'LineWidth', 2, 'TickDir', 'in');
ylim([0, 1.15]); box on;

% ===== Figure 4: T_rot/T_tr vs φ =====
fig4 = figure('Name','Fig4_Trot_Ttr','Position',[200,100,650,500]);

valid_r = R.ratio_rot_tr > 0;
semilogx(phi_data(valid_r), R.ratio_rot_tr(valid_r), 'ko-', ...
    'LineWidth', 2, 'MarkerSize', 9, 'MarkerFaceColor', [0.3 0.3 0.3]);
hold on;
yline(kappa, '--', 'Color', [0.8 0.1 0.1], 'LineWidth', 2);
text(phi_data(1)*1.2, kappa+0.03, '$\kappa = 2/5$ (equipartition)', ...
    'Interpreter','latex','FontSize',16,'Color',[0.8 0.1 0.1]);

xlabel('$\phi$', 'Interpreter', 'latex', 'FontSize', 22);
ylabel('$T_\mathrm{rot}/T_\mathrm{tr}$', 'Interpreter', 'latex', 'FontSize', 22);
set(gca, 'LineWidth', 2, 'TickDir', 'in');
ylim([0, 0.6]); box on;

% ===== Figure 5: 温度各向异性 Tx,Ty,Tz =====
fig5 = figure('Name','Fig5_Anisotropy','Position',[250,100,650,500]);
i_rep = ceil(n_data/2);

plot(dist_from_wall, R.Tx_shell(i_rep,:), 'ro-', 'MarkerSize', 5, ...
    'MarkerFaceColor', 'r', 'DisplayName', '$T_x$'); hold on;
plot(dist_from_wall, R.Ty_shell(i_rep,:), 'b^-', 'MarkerSize', 5, ...
    'MarkerFaceColor', 'b', 'DisplayName', '$T_y$');
plot(dist_from_wall, R.Tz_shell(i_rep,:), 'gs-', 'MarkerSize', 5, ...
    'MarkerFaceColor', 'g', 'DisplayName', '$T_z$');

xlabel('$z$ [cm]', 'Interpreter', 'latex', 'FontSize', 22);
ylabel('Temperature component [cm$^2$/s$^2$]', 'Interpreter', 'latex', 'FontSize', 22);
title(sprintf('$\\phi = %.3f$', phi_data(i_rep)), 'Interpreter', 'latex', 'FontSize', 20);
legend('Interpreter','latex','FontSize',16,'Location','northeast'); legend boxoff;
set(gca, 'LineWidth', 2, 'TickDir', 'in'); box on;

% ===== Figure 6: λ/L 相图 =====
fig6 = figure('Name','Fig6_Phase_Diagram','Position',[300,100,750,550]);

semilogx(phi_th, lambda_th/L, '-', 'Color', [0.2 0.2 0.7], 'LineWidth', 2.5); hold on;
if any(valid_m)
    semilogx(phi_data(valid_m), R.lambda_meas(valid_m)/L, 'o', ...
        'MarkerSize', 11, 'LineWidth', 1.8, ...
        'MarkerEdgeColor', [0.1 0.1 0.1], 'MarkerFaceColor', [0.3 0.7 0.3]);
end
yline(0.5, ':', 'Color', [0.3 0.3 0.3], 'LineWidth', 2);

% 着色区域
fill([phi_th(1) phi_th(end) phi_th(end) phi_th(1)], [0 0 0.5 0.5], ...
    [0.9 0.85 0.85], 'EdgeColor', 'none', 'FaceAlpha', 0.5, 'HandleVisibility', 'off');
fill([phi_th(1) phi_th(end) phi_th(end) phi_th(1)], [0.5 0.5 5 5], ...
    [0.85 0.85 0.95], 'EdgeColor', 'none', 'FaceAlpha', 0.5, 'HandleVisibility', 'off');
text(0.15, 0.15, '\textbf{Clustering}', 'Interpreter', 'latex', 'FontSize', 18, 'Color', [0.7 0.1 0.1]);
text(3e-4, 3, '\textbf{Homogeneous Gas}', 'Interpreter', 'latex', 'FontSize', 18, 'Color', [0.1 0.1 0.6]);

xlabel('$\phi$', 'Interpreter', 'latex', 'FontSize', 22);
ylabel('$\lambda / L$', 'Interpreter', 'latex', 'FontSize', 22);
legend({'Theory', 'Simulation'}, 'Interpreter', 'latex', 'FontSize', 16, ...
    'Location', 'northeast'); legend boxoff;
set(gca, 'LineWidth', 2, 'TickDir', 'in');
ylim([0, 5]); box on;

% ===== Figure 7: 密度涨落 σ²/⟨n⟩ =====
fig7 = figure('Name','Fig7_Density_Fluctuation','Position',[350,100,650,500]);

semilogy(phi_data, R.sigma2_n, 'bo-', 'LineWidth', 2, 'MarkerSize', 9, ...
    'MarkerFaceColor', [0.2 0.4 0.8]);
hold on;
yline(1, '--', 'Color', [0.4 0.4 0.4], 'LineWidth', 2);
text(phi_data(end)*0.7, 1.3, 'Ideal Gas', 'FontSize', 16, 'Color', [0.4 0.4 0.4]);

xlabel('$\phi$', 'Interpreter', 'latex', 'FontSize', 22);
ylabel('$\sigma_n^2 / \langle n \rangle$', 'Interpreter', 'latex', 'FontSize', 22);
set(gca, 'LineWidth', 2, 'TickDir', 'in'); box on;

% ===== Figure 8: 壳层 T_rot/T_tr(z) =====
fig8 = figure('Name','Fig8_Shell_Ratio','Position',[400,100,700,500]);

for i = 1:n_data
    ratio_sh = zeros(1, n_shells);
    for s = 1:n_shells
        if R.T_tr_shell(i,s)>0 && R.T_rot_shell(i,s)>0
            ratio_sh(s) = R.T_rot_shell(i,s) / (m*R.T_tr_shell(i,s));
        end
    end
    vs = ratio_sh>0;
    if any(vs)
        plot(dist_from_wall(vs), ratio_sh(vs), '-o', 'Color', cmap_custom(i,:), ...
            'MarkerSize', 4, 'MarkerFaceColor', cmap_custom(i,:), 'LineWidth', 1.2, ...
            'DisplayName', sprintf('%.3f', phi_data(i)));
        hold on;
    end
end
yline(kappa, '--', 'Color', [0.8 0.1 0.1], 'LineWidth', 2, 'HandleVisibility', 'off');
text(0.1, kappa+0.02, '$\kappa=2/5$', 'Interpreter','latex','FontSize',16,'Color',[0.8 0.1 0.1]);

xlabel('$z$ [cm]', 'Interpreter', 'latex', 'FontSize', 22);
ylabel('$T_\mathrm{rot}(z)/T_\mathrm{tr}(z)$', 'Interpreter', 'latex', 'FontSize', 22);
lgd = legend('Location','southeast','FontSize',10,'NumColumns',4);
title(lgd, '$\phi$', 'Interpreter', 'latex'); legend boxoff;
set(gca, 'LineWidth', 2, 'TickDir', 'in'); box on;

%% ========== 7. 汇总 ==========
fprintf('\n============================================================\n');
fprintf('  φ        T_tr      T_rot/T_tr    λ [cm]    λ/L     σ²/⟨n⟩\n');
fprintf('------------------------------------------------------------\n');
for i=1:n_data
    fprintf('%7.4f  %9.4f  %10.4f  %8.3f  %7.3f  %8.2f\n', ...
        phi_data(i), R.T_tr_global(i), R.ratio_rot_tr(i), ...
        R.lambda_meas(i), R.lambda_meas(i)/L, R.sigma2_n(i));
end
fprintf('============================================================\n');
fprintf('  φ_c (theory) ≈ %.4f\n', phi_c);
fprintf('  T_rot/T_tr ≈ %.2f (const) — 能量均分显著破缺\n', mean(R.ratio_rot_tr(valid_r)));
fprintf('============================================================\n');

%% 保存
save(fullfile(pos_vel_dir, 'Final_Lambda_Analysis.mat'), 'R', 'phi_data', ...
    'lambda_th', 'phi_th', 'phi_c', 'dist_from_wall', 'shell_centers', '-v7.3');
fprintf('结果已保存。\n');