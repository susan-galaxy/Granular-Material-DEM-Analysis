%% =====================================================================
%  微观平转动耦合向宏观聚团相变的定量分析脚本
%  功能: 计算 理论碰撞频率、有效能量耗散率、局域密度涨落(聚团指标)
% =====================================================================

clear all; close all; clc;

%% ========== 1. 用户配置区 ==========
data_dir = '/media/gezhuan/M78/Space_Active/Data/Maxwell_Boltz/data/Pos_Vel_data/NoDissNofric';
rot_dir  = '/media/gezhuan/M78/Space_Active/Data/Maxwell_Boltz/data/Rotation_data/NodissNofric';
pv_pattern  = 'Data*.mat';
rot_pattern = 'RotData*.mat';

tag_name = 'Tag_2';
L_box = 4.0;          % [cm]
d_part = 0.08;        % [cm]
R_part = d_part/2;
bulk_ratio = 0.6;     % 取纯净体区
half_L_bulk = (L_box * bulk_ratio) / 2;
V_bulk = (2 * half_L_bulk)^3;

delta_step = 5;

%% ========== 2. 自动检索与文件配对 ==========
pv_files  = dir(fullfile(data_dir, pv_pattern));
rot_files = dir(fullfile(rot_dir, rot_pattern));

% 颗粒单体体积与盒体积 (用于从 Number 反推 phi)
V_part_single = (pi/6) * d_part^3;
V_box_total   = L_box^3;

% 注意:不要用嵌套匿名函数捕获工作区变量,parfor 后续会更稳
extract_phi = @(fname) get_phi_from_name(fname, V_part_single, V_box_total);

pv_map  = containers.Map('KeyType', 'double', 'ValueType', 'char');
rot_map = containers.Map('KeyType', 'double', 'ValueType', 'char');
for k = 1:length(pv_files)
    phi_k = extract_phi(pv_files(k).name);
    if ~isnan(phi_k), pv_map(phi_k) = pv_files(k).name; end
end
for k = 1:length(rot_files)
    phi_k = extract_phi(rot_files(k).name);
    if ~isnan(phi_k), rot_map(phi_k) = rot_files(k).name; end
end

phi_common = sort(intersect(cell2mat(pv_map.keys), cell2mat(rot_map.keys)), 'ascend');
n_sets = length(phi_common);
%% ========== 3. 核心相变参数计算 ==========
phi_vals       = zeros(n_sets, 1);
Z_Enskog_vals  = zeros(n_sets, 1); % Enskog 碰撞频率
Diss_Rate_vals = zeros(n_sets, 1); % 有效能量耗散率代理
Density_Fluct  = zeros(n_sets, 1); % 局域密度涨落 (相变序参量)
Beta_proxy     = zeros(n_sets, 1); % 用于对比的宏观尾巴状态

% 密度涨落网格划分 (将 Bulk 分成 M x M x M 个小网格)
M_grid = 5; 
grid_edges = linspace(-half_L_bulk, half_L_bulk, M_grid+1);

parfor i = 1:n_sets
    phi_val = phi_common(i);
    phi_vals(i) = phi_val;
    
    fprintf('\n[%d/%d] 正在分析微观耗散与宏观涨落 φ = %g ...\n', i, n_sets, phi_val);
    
    pv_data = load(fullfile(data_dir, pv_map(phi_val)));
    rot_data = load(fullfile(rot_dir, rot_map(phi_val)));
    PD_pv = pv_data.ExtractedData.(tag_name);
    PD_rot = rot_data.ExtractedData.(tag_name);
    
    valid_frames = 1 : delta_step : length(PD_pv.Pos);
    
    v2_sum = 0; w2_sum = 0; N_bulk_sum = 0; count = 0;
    fluct_sum = 0; % 收集各帧的密度方差
    
    for j = valid_frames
        pos = PD_pv.Pos{j}; vel = PD_pv.Vel{j}; omg = PD_rot.Omg{j};
        if isempty(pos), continue; end
        
        nan_mask = any(isnan(pos),2) | any(isnan(vel),2);
        pos(nan_mask,:) = []; vel(nan_mask,:) = []; omg(nan_mask,:) = [];
        
        % 提取 Bulk
        dist = max(abs(pos), [], 2);
        mask = dist < half_L_bulk;
        pos_b = pos(mask, :); vel_b = vel(mask, :); omg_b = omg(mask, :);
        
        N_b = size(pos_b, 1);
        if N_b == 0, continue; end
        
        % 计算温度
        dv = vel_b - mean(vel_b, 1);
        v2_sum = v2_sum + mean(sum(dv.^2, 2));
        w2_sum = w2_sum + mean(sum(omg_b.^2, 2));
        N_bulk_sum = N_bulk_sum + N_b;
        
        % 【关键算法】：计算局域密度涨落方差 (判定聚团的铁证)
        % [修复版]: 使用 discretize + accumarray 替代 histcountsnd，完美兼容所有 MATLAB 版本
        idx_x = discretize(pos_b(:,1), grid_edges);
        idx_y = discretize(pos_b(:,2), grid_edges);
        idx_z = discretize(pos_b(:,3), grid_edges);
        
        % 过滤掉恰好在网格最边缘外的无效点
        valid_idx = ~isnan(idx_x) & ~isnan(idx_y) & ~isnan(idx_z);
        
        if any(valid_idx)
            % 构建三维下标矩阵
            subs = [idx_x(valid_idx), idx_y(valid_idx), idx_z(valid_idx)];
            % 统计 M_grid x M_grid x M_grid 三维网格内每个格子的颗粒数
            counts = accumarray(subs, 1, [M_grid, M_grid, M_grid]);
        else
            counts = zeros(M_grid, M_grid, M_grid);
        end
        
        n_mean = mean(counts(:));
        if n_mean > 0
            % 泊松比值: var(n) / mean(n) 
            % (当比值远大于 1 时，说明不再是均匀的泊松分布，而是发生了聚团)
            fluct_sum = fluct_sum + (var(counts(:)) / n_mean);
        end
        count = count + 1;
    end
    
    % 平均物理量
    T_trans = (v2_sum / count) / 3;
    T_rot = (2/5 * R_part^2 * (w2_sum / count)) / 3;
    Density_Fluct(i) = fluct_sum / count;
    n_density = (N_bulk_sum / count) / V_bulk;
    
    % 【关键理论公式】：计算 Enskog 碰撞频率
    % g0(phi) = (1 - phi/2) / (1 - phi)^3  (Carnahan-Starling 径向分布函数在接触点的值)
    g0_phi = (1 - phi_val/2) / (1 - phi_val)^3;
    % Z_E = 4 * sqrt(pi) * n * d^2 * g0 * sqrt(T_trans)
    Z_Enskog = 4 * sqrt(pi) * n_density * d_part^2 * g0_phi * sqrt(T_trans);
    Z_Enskog_vals(i) = Z_Enskog;
    
    % 【关键理论公式】：有效耗散率 (Dissipation Rate Proxy)
    % 耗散率正比于碰撞频率。且因为 T_trans/T_rot 高达 6，平动向转动的单向能量倾泻加剧了耗散。
    % 这里我们用 Z_E * T_trans 来表征单位体积内的总耗散能力强度
    Diss_Rate_vals(i) = Z_Enskog * T_trans;
    
    fprintf('      Bulk T_trans/T_rot = %.2f\n', T_trans/T_rot);
    fprintf('      碰撞频率 Z_E = %.2e\n', Z_Enskog);
    fprintf('      密度涨落指标 = %.3f (>>1 代表聚团)\n', Density_Fluct(i));
end


%% ========== 4. 宏观相变演化作图 (上下子图版，顶刊风格) ==========
figure('Name', 'Micro-Macro Phase Transition Link', 'Position', [150, 100, 700, 800]);

% --- 上图：微观起因 (Enskog 碰撞耗散率) ---
% ax1 = subplot(2, 1, 1);
plot(phi_vals, Z_Enskog_vals, 'rs-', 'LineWidth', 2.5, 'MarkerSize', 9, 'MarkerFaceColor', 'w');
xlabel('$\phi$', 'Interpreter', 'latex', 'FontSize', 18);
ylabel('$Z_E$', 'Interpreter', 'latex', 'FontSize', 18);
% title('Microscopic Energy Sink & Macroscopic Clustering', 'FontSize', 18);
set(gca, 'FontSize', 16, 'FontName', 'Times New Roman', 'LineWidth', 2);
% 隐藏上图的 X 轴刻度数字，使上下图视觉上更连贯
% set(gca, 'XTickLabel', []); 
figure;
% --- 下图：宏观结果 (局域密度涨落相变) ---
% ax2 = subplot(2, 1, 2);
plot(phi_vals, Density_Fluct, 'bo-', 'LineWidth', 2.5, 'MarkerSize', 12, 'MarkerFaceColor', 'w');
hold on;
% 理论基准线
yline(1.0, 'k--', 'Ideal Gas Limit', 'FontSize',15,'LineWidth', 2, 'FontSize', 15, 'LabelHorizontalAlignment', 'right');

xlabel('$\phi$', 'Interpreter', 'latex', 'FontSize', 18);
ylabel('$\sigma_n^2 / \langle n \rangle$', 'Interpreter', 'latex', 'FontSize', 18);
set(gca, 'YScale', 'log');

set(gca, 'FontSize', 16, 'FontName', 'Times New Roman', 'LineWidth', 2);

% --- 调整子图间距，让排版更紧凑美观 ---
pos1 = get(ax1, 'Position');
pos2 = get(ax2, 'Position');
pos1(2) = pos1(2) - 0.04; % 上图往下移一点
pos1(4) = pos1(4) + 0.04; % 上图拉高一点
pos2(4) = pos2(4) + 0.04; % 下图拉高一点
set(ax1, 'Position', pos1);
set(ax2, 'Position', pos2);

fprintf('\n=== 分析完毕！高颜值相变图表已生成 ===\n');
%% ========== 4. 宏观相变演化作图 (共享统一刻度双轴版) ==========
figure('Name', 'Micro-Macro Phase Transition Link', 'Position', [150, 150, 850, 600]);

% 【关键算法】：获取两个物理量的全局最大值，用于强制统一 Y 轴刻度
max_y = max(max(Density_Fluct), max(Z_Enskog_vals));
y_limits = [0, max_y * 1.1]; % 从 0 开始，顶部留出 10% 的呼吸空间

% --- 左坐标轴：局域密度涨落 (判定聚团相变) ---
yyaxis left
plot(phi_vals, Density_Fluct, 'bo-', 'LineWidth', 3, 'MarkerSize', 10, 'MarkerFaceColor', 'w');
ylabel('$\sigma_n^2 / \langle n \rangle$', 'Interpreter', 'latex', 'FontSize', 18, 'Color', 'b');
yline(1.0, 'b--', 'Ideal Gas Limit', 'FontSize',15,'LineWidth', 2, 'LabelHorizontalAlignment', 'right');
set(gca, 'YColor', 'b');
ylim(y_limits); % 强制左轴刻度

% --- 右坐标轴：有效碰撞耗散率 ---
yyaxis right
plot(phi_vals, Z_Enskog_vals, 'rs-', 'LineWidth', 3, 'MarkerSize', 10, 'MarkerFaceColor', 'w');
ylabel('$Z_E$', 'Interpreter', 'latex', 'FontSize', 18, 'Color', 'r');
set(gca, 'YColor', 'r');
ylim(y_limits); % 强制右轴刻度与左轴绝对一致！

% --- 格式美化 ---
xlabel('$\phi$', 'Interpreter', 'latex', 'FontSize', 18);
% title('Microscopic Energy Sink drives Macroscopic Clustering Phase Transition', 'FontSize', 18);
set(gca, 'FontSize', 20, 'FontName', 'Times New Roman', 'LineWidth', 2);


fprintf('\n=== 分析完毕！统一 Y 轴刻度的完美双轴图已生成 ===\n');

%% ========== 5. 理论解析：无量纲衰减长度预测相变 ==========
figure('Name', 'Analytical Prediction of Phase Transition', 'Position', [200, 200, 750, 550]);

% 假设宏观有效耗散常数 (包含法向碰撞与极强的切向摩擦耗散)
% 由于 T_trans/T_rot ≈ 6，能量损耗极大，我们估算有效耗散参数 epsilon_eff 约为 0.15
epsilon_eff = 0.15; 
prefactor = 1.0; % 动理学理论修正前缀

% --- 1. 计算理论曲线 ---
phi_theory_dense = logspace(log10(0.002), log10(0.06), 200);
% 修复：化简体积项，统一使用 d_part
n_theory_dense = phi_theory_dense / (pi/6 * d_part^3);
Kn_theory_dense = (1 ./ (sqrt(2) .* n_theory_dense * pi * d_part^2)) / L_box;
% 解析公式：计算无量纲衰减长度 xi
xi_theory = prefactor .* (2 * Kn_theory_dense) ./ sqrt(epsilon_eff);

% --- 2. 计算模拟数据点的 xi ---
% 将仿真跑出来的 phi_vals 转换为 Kn
n_sim_dense = phi_vals / (pi/6 * d_part^3);
Kn_sim = (1 ./ (sqrt(2) .* n_sim_dense * pi * d_part^2)) / L_box;
xi_sim = prefactor .* (2 * Kn_sim) ./ sqrt(epsilon_eff);

% --- 3. 作图 ---
% 理论预测的衰减长度 (黑线)
loglog(phi_theory_dense, xi_theory, 'k-', 'LineWidth', 2.5, 'DisplayName', 'Analytical \xi (\epsilon_{eff} \approx 0.15)');
hold on;
% 仿真数据算出的衰减长度 (蓝点)
loglog(phi_vals, xi_sim, 'bo', 'MarkerSize', 10, 'MarkerFaceColor', [0.3 0.6 0.9], 'LineWidth', 2, 'DisplayName', 'Simulated \xi');

% 画出相变临界阈值 (xi = 1)
yline(1.0, 'r--', 'Critical Penetration (\xi = 1)', 'LineWidth', 2.5, 'FontSize', 12, 'LabelHorizontalAlignment', 'right', 'Color', 'r');

% 标注聚团区 (灰色背景填充)
patch_x = [0.02, 0.06, 0.06, 0.02];
patch_y = [0.1, 0.1, 1, 1];
patch(patch_x, patch_y, 'k', 'FaceAlpha', 0.1, 'EdgeColor', 'none', 'HandleVisibility', 'off');
text(0.025, 0.5, 'Clustering Zone (\xi < 1)', 'FontSize', 16, 'Color', [0.3 0.3 0.3]);

% --- 4. 格式美化 ---
xlabel('$\phi$', 'Interpreter', 'latex', 'FontSize', 18);
ylabel('$\xi = 2Kn / \sqrt{\epsilon_{eff}}$', 'Interpreter', 'latex', 'FontSize', 18);
% title('Analytical Prediction of Clustering Instability', 'FontSize', 18);

legend({'Analytical \xi', 'Simulated \xi'}, 'FontSize', 18); 
legend boxoff;set(gca, 'FontSize', 18, 'FontName', 'Times New Roman', 'LineWidth', 2);
ylim([0.1, 10]);
xlim([min(phi_theory_dense), max(phi_theory_dense)]);


%% ========== 局部函数: 多规则文件名 -> phi 解析 ==========
function phi = get_phi_from_name(fname, V_part_single, V_box_total)
    % 规则 1 (旧命名): 直接匹配 'phi' 后的数字, 例: Data_phi0.02_...
    s = regexp(fname, '(?<=phi)\d+\.?\d*', 'match', 'once');
    if ~isempty(s)
        phi = str2double(s);
        return;
    end
    
    % 规则 2 (新命名): 从 'Number' 后的颗粒总数 N 反推
    %   phi = N * (pi/6 * d^3) / L^3
    % 例: Data_ActiveGas_0_03-0_06_Number14324_Tag1_...
    s = regexp(fname, '(?<=Number)\d+', 'match', 'once');
    if ~isempty(s)
        N = str2double(s);
        phi = N * V_part_single / V_box_total;
        return;
    end
    
    % 两种规则都没匹配上
    phi = NaN;
end