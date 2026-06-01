%% =====================================================================
%  微重力颗粒气体：从底层 DEM 数据到相变分岔图的全自动闭环预测
%  1. 自动读取 Pos_Vel_data (.mat)
%  2. 提取序参量 T_center/T_wall 
%  3. 提取 3D 空间局部密度涨落 \sigma_n^2 / <n>
%  4. 求解双流体理论 ODE 并自动寻找相变分岔点
%  5. 绘制出版级双 Y 轴相图
% =====================================================================
clear all; close all; clc;

%% ========== 1. 用户配置区 (匹配 Kn_Tr_Tt_v1 逻辑) ==========

% --- 1.1 数据路径 & 通配符 ---
pos_vel_dir = 'G:\Space_Active\Data\Maxwell_Boltz\data\Pos_Vel_data\NoDissNofric';
pv_pattern  = 'Data*.mat';

% --- 1.2 系统与理论参数 ---
tag_name     = 'Tag_2';
d            = 0.08;           % 颗粒直径 [cm]
e_n          = 0.94;           % 法向恢复系数
L            = 4.0;            % 盒子尺寸 [cm]
half_L       = L / 2;
S0           = 0.25;            % 理论源项注入强度 (根据之前提取设定)
lambda_b_real= 0.6;            % 弹道穿透深度 [cm] (根据之前提取设定)

% --- 1.3 数据提取设置 ---
delta_step   = 1;              % 抽样帧间隔
n_shells     = 20;             % 温度壳层数量
shell_edges  = linspace(0, half_L, n_shells+1);
grid_bins    = 10;             % 用于密度涨落的 3D 网格数 (10x10x10 = 1000个子空间)
grid_edges   = linspace(-half_L, half_L, grid_bins+1); 

%% ========== 2. 自动检索与配对数据文件 ==========
fprintf('=== 第一步: 检索并解析 DEM 数据 ===\n');

pv_files  = dir(fullfile(pos_vel_dir, pv_pattern));
if isempty(pv_files)
    error('在 %s 中未找到匹配 "%s" 的文件！', pos_vel_dir, pv_pattern);
end

% 提取 phi
extract_phi = @(fname) str2double(regexp(fname, '(?<=phi)\d+\.?\d*', 'match', 'once'));
pv_map = containers.Map('KeyType', 'double', 'ValueType', 'char');
for k = 1:length(pv_files)
    phi_k = extract_phi(pv_files(k).name);
    if ~isnan(phi_k), pv_map(phi_k) = pv_files(k).name; end
end

phi_dem_all = sort(cell2mat(pv_map.keys));
n_sets = length(phi_dem_all);
fprintf('自动匹配到 %d 组体积分数数据。\n', n_sets);

% 预分配数组
Tc_dem = zeros(1, n_sets);
fluc_dem = zeros(1, n_sets);

%% ========== 3. 核心循环：提取温度与密度涨落 ==========
parfor i = 1:n_sets
    phi_val = phi_dem_all(i);
    fname = pv_map(phi_val);
    fprintf('[%d/%d] 正在处理 φ = %.4f ... ', i, n_sets, phi_val);
    
    ld = load(fullfile(pos_vel_dir, fname));
    PD = ld.ExtractedData.(tag_name);
    
    total_frames = length(PD.Pos);
    valid_frames = 1:delta_step:total_frames;
    
    sv2 = zeros(n_shells, 1);
    sc  = zeros(n_shells, 1);
    fluc_frames = zeros(length(valid_frames), 1);
    
    frame_idx = 1;
    for j = valid_frames
            % 强制将 single 转换为 double，满足 histcounts3 和高精度计算的要求
            p = double(PD.Pos{j}); 
            v = double(PD.Vel{j});
    
            if isempty(p), continue; end
            b = any(isnan(p),2) | any(isnan(v),2); 
            p(b,:) = []; v(b,:) = [];
            if isempty(p), continue; end
    
            vf = v - mean(v, 1);
        
        % ... 下面的代码完全保持不变 ...
        
        % --- 3.1 计算壳层温度分布 ---
        ch = max(abs(p), [], 2);
        si = discretize(ch, shell_edges);
        for s = 1:n_shells
            m = (si == s);
            if any(m)
                sv2(s) = sv2(s) + sum(sum(vf(m,:).^2, 2));
                sc(s)  = sc(s) + sum(m);
            end
        end
        
% --- 3.2 计算 3D 空间局部密度涨落 ---
        % MATLAB 没有原生 histcounts3，我们使用 discretize + 线性索引实现极其高效的 3D 统计
        idx_x = discretize(p(:,1), grid_edges);
        idx_y = discretize(p(:,2), grid_edges);
        idx_z = discretize(p(:,3), grid_edges);
        
        % 剔除刚好在边界外产生的 NaN 游离颗粒
        valid_pts = ~isnan(idx_x) & ~isnan(idx_y) & ~isnan(idx_z);
        idx_x = idx_x(valid_pts);
        idx_y = idx_y(valid_pts);
        idx_z = idx_z(valid_pts);
        
        if ~isempty(idx_x)
            % 将 3D 网格坐标转换为 1D 线性索引
            total_bins = grid_bins^3;
            lin_idx = sub2ind([grid_bins, grid_bins, grid_bins], idx_x, idx_y, idx_z);
            
            % 统计每个 3D 网格内的颗粒数 (由于是整数索引，边条设置为 0.5 到 total_bins+0.5)
            N_vec = histcounts(lin_idx, 0.5 : 1 : (total_bins + 0.5));
            
            mean_N = mean(N_vec);
            if mean_N > 0
                % \sigma_n^2 / <n> (泊松分布气相极限下为 1，液相聚束时急剧上升)
                fluc_frames(frame_idx) = var(N_vec) / mean_N;
            end
        end
        frame_idx = frame_idx + 1;
    end
    
    % --- 3.3 归纳当前 phi 的序参量 ---
    T_shell = zeros(n_shells, 1);
    for s = 1:n_shells
        if sc(s) > 0, T_shell(s) = sv2(s) / (3 * sc(s)); end
    end
    
    % 寻找最靠中心和最靠壁面的有效温度
    valid_s = find(sc > 0);
    if ~isempty(valid_s)
        T_center = T_shell(valid_s(1));   % 最靠近中心的壳层
        T_wall   = T_shell(valid_s(end)); % 最靠近壁面的壳层
        Tc_dem(i) = T_center / T_wall;
    end
    
    % 取所有时间帧的密度涨落平均值
    fluc_frames = fluc_frames(fluc_frames > 0);
    if ~isempty(fluc_frames)
        fluc_dem(i) = mean(fluc_frames);
    end
    
    fprintf(' Tc/Tw = %.3f, Fluc = %.2f\n', Tc_dem(i), fluc_dem(i));
end

%% ========== 4. 理论模型扫描预测 ==========
fprintf('\n=== 第二步: 运行连续性双流体 ODE 理论扫描 ===\n');
chi_f = @(phi) (1-phi/2)./(1-phi).^3;
lam_f = @(phi) d./(6*phi.*chi_f(phi)*sqrt(1-e_n^2));

% 自动设定扫描范围覆盖实验数据
phi_theory = logspace(log10(min(phi_dem_all)*0.8), log10(max(phi_dem_all)*1.2), 100);
Tc_theory  = zeros(size(phi_theory));

for k = 1:length(phi_theory)
    lam_0 = lam_f(phi_theory(k));
    xi_max = half_L / lam_0;
    lambda_b_nondim = lambda_b_real / lam_0;
    
    ode_an = @(xi, w) [w(2); 1/(2*max(w(1), 1e-6)) - S0 * (exp(-xi/lambda_b_nondim) + exp(-(2*xi_max-xi)/lambda_b_nondim)) / (2*max(w(1), 1e-6))];
    ul = 1e-6; uh = 1.0; success = false;
    
    for it = 1:80
        ug = (ul + uh) / 2; 
        try
            options = odeset('RelTol', 1e-6, 'AbsTol', 1e-8, 'MaxStep', xi_max/50);
            [~, ws_an] = ode45(ode_an, [xi_max, 0], [ug, 0], options);
            u_w = ws_an(end, 1);
            if ~isreal(u_w) || isnan(u_w), ul = ug; continue; end
            if u_w < 1, ul = ug; else, uh = ug; end
            if abs(u_w - 1) < 1e-4, success = true; break; end
        catch, ul = ug; end
    end
    if success, Tc_theory(k) = ug^2; else, Tc_theory(k) = 0; end
end

phi_c_idx = find(Tc_theory < 1e-4, 1, 'first');
phi_c_theory = phi_theory(phi_c_idx);
fprintf('理论预测临界相变点: phi_c ≈ %.4f\n', phi_c_theory);

%% ========== 5. 绘制顶刊级双 Y 轴相图 ==========
fprintf('\n=== 第三步: 绘制物理相图 ===\n');
figure('Position', [150, 150, 850, 550], 'Color', 'w');

% 【左轴：温度序参量】
yyaxis left;
plot(phi_theory, Tc_theory, 'b-', 'LineWidth', 3, 'DisplayName', 'Theory: Steady-state $T_c$'); hold on;
plot(phi_dem_all, Tc_dem, 'bo', 'MarkerFaceColor', 'b', 'MarkerSize', 8, 'DisplayName', 'DEM: $T_c$');
ylabel('$T_\mathrm{center}/T_\mathrm{wall}$', 'Interpreter', 'latex', 'FontSize', 22);
ylim([-0.05, 1.05]); 
yticks(0:0.2:1.0);  
ax = gca; ax.YColor = 'b'; 

% 【右轴：密度涨落】
yyaxis right;
plot(phi_dem_all, fluc_dem, 'rs-', 'LineWidth', 2.5, 'MarkerFaceColor', 'r', 'MarkerSize', 8, 'DisplayName', 'DEM: Fluctuation $\sigma_n^2/\langle n \rangle$');
ylabel('Density Fluctuation', 'Interpreter', 'latex', 'FontSize', 22);
set(gca, 'YScale', 'log');

% 动态调整右轴显示范围以适配真实数据
max_fluc = max(fluc_dem) * 2;
ylim([0.8, max(10, max_fluc)]); 
ax = gca; ax.YColor = 'r'; 

yline(1, 'r:', 'LineWidth', 2, 'HandleVisibility', 'off');
% text(min(phi_theory)*1.1, 1.2, 'Ideal Gas Limit', 'Interpreter', 'latex', 'FontSize', 14, 'Color', [0.8 0 0]);

% 【公共 X 轴设置】
xlabel('$\phi$', 'Interpreter', 'latex', 'FontSize', 22);
set(gca, 'XScale', 'log');
xlim([min(phi_theory), max(phi_theory)]);

% 动态设置X轴标签，让其对数刻度好看
min_pow = floor(log10(min(phi_theory)));
max_pow = ceil(log10(max(phi_theory)));
xt = [];
for p = min_pow:max_pow
    xt = [xt, 10^p, 2*10^p, 5*10^p];
end
xt = xt(xt >= min(phi_theory) & xt <= max(phi_theory));
xticks(xt);
xticklabels(arrayfun(@num2str, xt, 'UniformOutput', false));

% % 【标注理论临界点】
% if ~isnan(phi_c_theory)
%     xline(phi_c_theory, 'k--', 'LineWidth', 2, 'HandleVisibility', 'off');
%     % 将文字动态靠在虚线左侧或右侧
%     % text(phi_c_theory * 1.05, ax.YLim(2)^0.5, sprintf('Theory $\\phi_c \\approx %.3f$', phi_c_theory), ...
%     %     'Interpreter', 'latex', 'FontSize', 18, 'Color', 'k');
% end

% 【美化与图例】
lgd = legend('Location', 'northwest', 'FontSize', 15, 'Interpreter', 'latex');
legend boxoff;
set(gca, 'FontSize', 18, 'FontName', 'Times New Roman', 'TickDir', 'in', 'LineWidth', 1.5);
box on;

fprintf('\n✅ 处理完成，图表已生成！\n');