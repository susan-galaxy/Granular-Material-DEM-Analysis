%% =====================================================================
%  无重力体系特供：基于速度绝对突变的真实平均自由程 (λ) 与 Kn 计算
%  核心思想: 颗粒在无重力下做匀速直线运动。dv ≠ 0 即代表发生碰撞。
% =====================================================================
clear all; close all; clc;

%% ========== 1. 用户配置区 ==========
pos_vel_dir = 'G:\Space_Active\Data\Maxwell_Boltz\data\Pos_Vel_data\NoDissNofric';
pv_pattern  = 'Data_*.mat';

tag_name     = 'Tag_2';        
d_particle   = 0.08;
L_box        = 4.0;            
bulk_ratio   = 0.6;  
half_L_bulk  = (L_box * bulk_ratio) / 2; 

% 【极其重要】必须逐帧比对，否则会漏掉子步里的多次碰撞
delta_step = 1;  

% 速度突变容差 (为了屏蔽浮点数计算误差，设定一个极小的阈值)
% dv^2 > 1e-6 认为发生了物理碰撞
vel_tolerance = 1e-6; 

fprintf('=== 启动无重力体系真实轨迹追踪 (Bulk Ratio = %.2f) ===\n', bulk_ratio);

%% ========== 2. 自动检索数据 ==========
pv_files  = dir(fullfile(pos_vel_dir, pv_pattern));
extract_phi = @(fname) str2double(regexp(fname, '(?<=phi)\d+\.?\d*', 'match', 'once'));


pv_map = containers.Map('KeyType', 'double', 'ValueType', 'char');
for k = 1:length(pv_files)
    phi_k = extract_phi(pv_files(k).name);
    if ~isnan(phi_k), pv_map(phi_k) = pv_files(k).name; end
end

phi_vals = sort(cell2mat(pv_map.keys));
n_sets = length(phi_vals);

%% ========== 3. 核心算法：轨迹累加与动量突变检测 ==========
Kn_true_vals = zeros(n_sets, 1);
lambda_true_vals = zeros(n_sets, 1);

parfor i = 1:n_sets
    phi_val = phi_vals(i);
    fprintf('\n[%d/%d] φ = %.4f | 正在追踪颗粒真实自由飞行轨迹...\n', i, n_sets, phi_val);
    
    pv_data = load(fullfile(pos_vel_dir, pv_map(phi_val)));
    PD = pv_data.ExtractedData.(tag_name);
    
    num_frames = length(PD.Pos);
    if num_frames < 2
        warning('帧数过少，跳过');
        continue;
    end
    
    % 获取最大颗粒数，初始化每个颗粒的"自由飞翔里程表"
    max_particles = size(PD.Pos{1}, 1); 
    accumulated_dist = zeros(max_particles, 1); 
    
    bulk_free_paths = []; % 记录在 Bulk 区内【终结】的有效自由程
    
    % 逐帧比对
    for j = 2 : delta_step : num_frames
        pos_curr = PD.Pos{j};   vel_curr = PD.Vel{j};
        pos_prev = PD.Pos{j-1}; vel_prev = PD.Vel{j-1};
        
        if isempty(pos_curr) || isempty(pos_prev), continue; end
        
        min_N = min(size(pos_curr, 1), size(pos_prev, 1));
        
        % 1. 积分累加本帧的实际空间位移 (惯性飞行距离)
        dist_step = sqrt(sum((pos_curr(1:min_N,:) - pos_prev(1:min_N,:)).^2, 2));
        accumulated_dist(1:min_N) = accumulated_dist(1:min_N) + dist_step;
        
        % 2. 绝杀碰撞检测：速度矢量差的平方大于容差，绝对是撞了！
        dv_sq = sum((vel_curr(1:min_N,:) - vel_prev(1:min_N,:)).^2, 2);
        is_collision = dv_sq > vel_tolerance;
        
        if any(is_collision)
            col_idx = find(is_collision);
            
            % 提取发生碰撞的颗粒的当前位置
            col_pos = pos_curr(col_idx, :);
            
            % 判断这些碰撞是否发生在 Bulk 体区内
            dist_from_center = max(abs(col_pos), [], 2);
            in_bulk = dist_from_center < half_L_bulk;
            
            % 对于在 Bulk 内发生碰撞的颗粒，把它之前的飞行里程存入结果库
            valid_bulk_idx = col_idx(in_bulk);
            if ~isempty(valid_bulk_idx)
                bulk_free_paths = [bulk_free_paths; accumulated_dist(valid_bulk_idx)];
            end
            
            % 【核心重置】：无论碰撞发生在墙壁还是体区，只要撞了，里程表必须清零！
            accumulated_dist(col_idx) = 0;
        end
    end
    
    % 清除由于持续接触(Enduring contacts)造成的微小数值抖动，保留物理有效自由程
    bulk_free_paths = bulk_free_paths(bulk_free_paths > d_particle * 0.01); 
    
    if ~isempty(bulk_free_paths)
        lambda_true = mean(bulk_free_paths);
    else
        lambda_true = NaN; 
    end
    
    lambda_true_vals(i) = lambda_true;
    Kn_true_vals(i) = lambda_true / L_box;
    
    fprintf('  -> 捕捉到 %d 次 Bulk 区有效碰撞\n', length(bulk_free_paths));
    fprintf('  -> 真实局域平均自由程 λ_true = %.4f cm\n', lambda_true);
    fprintf('  -> 真实 Knudsen 数 Kn_true = %.5f\n', Kn_true_vals(i));
end

%% ========== 4. 绘图：理论公式 vs 真实的降维打击 ==========
figure('Name', 'True Lagrangian Mean Free Path', 'Position', [150, 150, 700, 550]);

% 1. 画理论公式 (理想气体均匀假设)
phi_theory = logspace(log10(0.002), log10(0.06), 100);
lambda_theory = 1 ./ (sqrt(2) .* (phi_theory / (pi/6 * d_particle^3)) * pi * d_particle^2);
Kn_theory = lambda_theory / L_box;
loglog(phi_theory, Kn_theory, 'k-', 'LineWidth', 2.5, 'DisplayName', 'Kinetic Theory ');
hold on;

% 2. 画通过轨迹突变追踪出来的纯正真实 Kn
loglog(phi_vals, Kn_true_vals, 'ro-', 'MarkerSize', 10, 'MarkerFaceColor', 'w', ...
    'LineWidth', 2.5, 'DisplayName', 'True Computed $Kn_{bulk}$ ');

xlabel('$\phi$', 'Interpreter', 'latex', 'FontSize', 18);
ylabel('$Kn = \lambda_{true} / L_{box}$', 'Interpreter', 'latex', 'FontSize', 18);
% title('Breakdown of Kinetic Theory in Clustering Phase', 'FontSize', 18);

legend('Interpreter', 'latex', 'Location', 'northeast', 'FontSize', 14); legend boxoff;
set(gca, 'FontSize', 18, 'FontName', 'Times New Roman', 'LineWidth', 2);


%% ========== 5. 终极审判：基于真实 Kn 的相变预测 ==========
figure('Name', 'True Dimensionless Penetration', 'Position', [900, 150, 700, 550]);

epsilon_eff = 0.15; 
prefactor = 1.0; 

xi_theory = prefactor .* (2 * Kn_theory) ./ sqrt(epsilon_eff);
xi_true = prefactor .* (2 * Kn_true_vals) ./ sqrt(epsilon_eff);

% 黑线：基于均匀假设的衰减长度
loglog(phi_theory, xi_theory, 'k-', 'LineWidth', 2.5, 'DisplayName', 'Analytical $\xi$ (Uniform Assumption)');
hold on;
% 红线：真实的衰减长度
loglog(phi_vals, xi_true, 'ro-', 'MarkerSize', 10, 'MarkerFaceColor', 'w', ...
    'LineWidth', 2.5, 'DisplayName', 'True Computed $\xi$ (Real Local Physics)');

% 临界阈值
yline(1.0, 'r--', 'Critical Penetration (\xi = 1)', 'LineWidth', 2.5, 'FontSize', 16, 'LabelHorizontalAlignment', 'right');

% 标注聚团区
patch_x = [0.015, 0.06, 0.06, 0.015];
patch_y = [0.01, 0.01, 1, 1];
patch(patch_x, patch_y, 'k', 'FaceAlpha', 0.1, 'EdgeColor', 'none', 'HandleVisibility', 'off');
text(0.02, 0.3, 'Clustering Zone ($\xi < 1$)', 'Interpreter', 'latex', 'FontSize', 16, 'Color', [0.3 0.3 0.3]);

xlabel('$\phi$', 'Interpreter', 'latex', 'FontSize', 18);
ylabel('$\xi = 2Kn_{true} / \sqrt{\epsilon_{eff}}$', 'Interpreter', 'latex', 'FontSize', 18);
% title('Prediction of Phase Transition (First-Principles)', 'FontSize', 18);

legend('Interpreter', 'latex', 'Location', 'northeast', 'FontSize', 14); legend boxoff;
set(gca, 'FontSize', 18, 'FontName', 'Times New Roman', 'LineWidth', 2);
ylim([0.05, 10]); xlim([0.002, 0.06]); 

fprintf('\n=== 🚀 所有的真实追踪计算已完成！ ===\n');

%%
% =========================================================
% 【新增】自动拟合并标注 phi > 0.02 区间的斜率
% =========================================================

% 1. 提取 phi > 0.02 的数据段
fit_mask = phi_vals > 0.02;
phi_fit = phi_vals(fit_mask);
Kn_fit = Kn_true_vals(fit_mask);

if length(phi_fit) >= 2
    % 2. 在双对数坐标下进行线性拟合: log10(Kn) = k * log10(phi) + b
    log_phi = log10(phi_fit);
    log_Kn = log10(Kn_fit);
    p = polyfit(log_phi, log_Kn, 1);
    slope = p(1); % 提取斜率 (幂律指数)
    
    % 3. 生成拟合直线的绘图数据
    phi_plot = logspace(log10(min(phi_fit)*0.9), log10(max(phi_fit)*1.1), 50);
    Kn_plot = 10.^(polyval(p, log10(phi_plot)));
    
    % 4. 绘制拟合辅助线 (蓝色虚线)
    loglog(phi_plot, Kn_plot, 'b--', 'LineWidth', 2, 'DisplayName', 'Power-law Fit ($\phi > 0.02$)');
    
    % 5. 在图上合适的位置标注斜率文本
    % 动态计算文本坐标，使其位于拟合线附近
    text_x = phi_fit(floor(length(phi_fit)/2)) * 1.1; 
    text_y = Kn_fit(floor(length(Kn_fit)/2)) * 1.5; 
    text(text_x, text_y, sprintf('Slope $\\approx %.2f$', slope), ...
        'Interpreter', 'latex', 'FontSize', 16, 'Color', 'b', 'Rotation', slope*10); 
    
    % 终端输出结果
    fprintf('\n---> 发现拐点! phi > 0.02 区间的对数斜率为: %.3f <---\n', slope);
else
    fprintf('\n---> [提示] phi > 0.02 的数据点不足，无法拟合斜率 <---\n');
end

% 刷新图例 (包含新加的拟合线)
legend('Interpreter', 'latex', 'Location', 'northeast', 'FontSize', 14); legend boxoff;