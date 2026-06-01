% =====================================================================
%  自动化 HDF5 颗粒数据提取脚本 (转动信息) - 并行优化版 v3
%
%  v3 在 v2 基础上新增 (假设颗粒数恒定):
%   ★ PTag 缓存: 只在首帧读一次 /PTag,缓存 idx,跳过 2999 次 /PTag 读取
%   ★ 首末帧 PTag 校验
%
%  累计相对原版的优化:
%   1. parfor 并行读取 HDF5 文件
%   2. PTag 缓存 (节省 ~17% h5read 调用,5个数据集变4个)
%   3. Inertia 维度判断提到 parfor 外
%   4. 逻辑索引向量化过滤文件名
%   5. parfor 友好的预分配切片数组
%   6. 并行池跨文件夹复用
%
%  输出 ExtractedData 结构与 v1 完全一致,下游脚本无需任何修改。
% =====================================================================

clearvars;
close all;
clc;

%% ========== 1. 核心配置区 ==========
% cd /mnt/nvme/scratch/extract_input
folder_pattern = 'NoDiss*'; 

output_folder = '/media/gezhuan/M78/Space_Active/Data/Maxwell_Boltz/data/Rotation_data/NodissNofric';

start_step = 1;
end_step   = 3000;
delta_step = 1;

tags_to_extract = -2;

% 并行工作进程数 (NVMe SSD: 8~12; SATA SSD: 6~8; HDD: 2~4)
num_workers = 4;

% 首末帧 PTag 一致性校验
verify_ptag_consistency = true;

%% ========== 2. 启动并行池 ==========

pool = gcp('nocreate');
if isempty(pool)
    fprintf('启动并行池 (%d workers)...\n', num_workers);
    parpool('local', num_workers);
elseif pool.NumWorkers ~= num_workers
    fprintf('重启并行池: %d -> %d workers\n', pool.NumWorkers, num_workers);
    delete(pool);
    parpool('local', num_workers);
else
    fprintf('复用已有并行池 (%d workers)\n', pool.NumWorkers);
end

%% ========== 3. 文件夹检索 ==========

items = dir(folder_pattern);
is_real_folder = [items.isdir] & ~ismember({items.name}, {'.','..'});
items = items(is_real_folder);

if isempty(items)
    error('未找到匹配模式 "%s" 的文件夹,请检查路径或拼写!', folder_pattern);
end

target_folders = arrayfun(@(s) fullfile(s.folder, s.name), items, ...
                          'UniformOutput', false);
fprintf('成功匹配到 %d 个文件夹,准备开始提取【转动数据】...\n', ...
        length(target_folders));

if ~exist(output_folder, 'dir')
    mkdir(output_folder);
    fprintf('已创建数据输出目录: %s\n', output_folder);
end

n_tags = length(tags_to_extract);

%% ========== 4. 主循环 ==========

for f_idx = 1:length(target_folders)
    current_folder = target_folders{f_idx};
    [~, name_part, ext_part] = fileparts(current_folder);
    folder_name = [name_part, ext_part]; 
    
    fprintf('\n======================================================\n');
    fprintf('▶ 正在处理 [%d/%d]: %s\n', f_idx, length(target_folders), folder_name);
    
    % --- 4.1 获取并过滤 H5 文件 ---
    file_pattern = fullfile(current_folder, 'Space_s_*.h5');
    all_files = dir(file_pattern);
    if isempty(all_files)
        warning('文件夹 %s 中未找到 .h5 文件,跳过', folder_name);
        continue;
    end
    mask = ~contains({all_files.name}, '_bf_');
    df_normal = all_files(mask);
    
    total_files = length(df_normal);
    if total_files == 0
        warning('文件夹 %s 中未找到符合条件的 .h5 文件,跳过', folder_name);
        continue;
    end
    
    calc_end = min(end_step, total_files);
    stepp = start_step:delta_step:calc_end;
    num_steps = length(stepp);
    fprintf('  找到 %d 个文件,计划处理 %d 个时间步 (从 %d 到 %d)\n', ...
            total_files, num_steps, start_step, calc_end);
    
    % 组装文件全路径
    file_names_full = cell(num_steps, 1);
    for k = 1:num_steps
        file_names_full{k} = fullfile(current_folder, df_normal(stepp(k)).name);
    end
    
    % --- 4.2 ★ PTag 缓存 + Inertia 维度判断 (只读首帧) ★ ---
    fprintf('  读取首帧 PTag 并缓存索引...\n');
    first_ptag    = h5read(file_names_full{1}, '/PTag');
    first_inertia = h5read(file_names_full{1}, '/Inertia');
    
    idx_cache = cell(n_tags, 1);
    for t = 1:n_tags
        idx_cache{t} = find(first_ptag == tags_to_extract(t));
        fprintf('    Tag %d: %d 个颗粒\n', tags_to_extract(t), length(idx_cache{t}));
    end
    
    is_3d_inertia = (length(first_inertia) == 3 * length(first_ptag));
    if is_3d_inertia
        fprintf('  转动惯量维度: 3D (三方向张量)\n');
    else
        fprintf('  转动惯量维度: 1D (标量)\n');
    end
    
    % 校验末帧 PTag
    if verify_ptag_consistency && num_steps > 1
        last_ptag = h5read(file_names_full{end}, '/PTag');
        if ~isequal(first_ptag, last_ptag)
            error(['首末帧 PTag 不一致! 该体系颗粒数/Tag 发生变化,不能使用 PTag 缓存。' ...
                   '请使用 v2 版本(每步读 PTag)。']);
        end
        fprintf('  ✔ 首末帧 PTag 一致性校验通过\n');
    end
    
    % --- 4.3 预分配切片输出 ---
    omg_cells = cell(num_steps, n_tags);
    am_cells  = cell(num_steps, n_tags);
    aa_cells  = cell(num_steps, n_tags);
    in_cells  = cell(num_steps, n_tags);
    ts_arr    = zeros(num_steps, n_tags);
    
    stepp_local = stepp;
    is_3d_local = is_3d_inertia;
    
    % --- 4.4 并行读取 (每步只读 4 个数据集,不再读 /PTag) ---
    fprintf('  开始并行读取 (%d workers)...\n', num_workers);
    tic;
    
    parfor step_idx = 1:num_steps
        filename = file_names_full{step_idx};
        
        % [读取转动数据] - 跳过 /PTag (已缓存)
        omg     = h5read(filename, '/PAngVelocity');
        ang_mom = h5read(filename, '/PAngMomentum');
        ang_acc = h5read(filename, '/PAngacceleration');
        inertia = h5read(filename, '/Inertia');
        
        % 数据解包
        omgx = omg(1:3:end);     omgy = omg(2:3:end);     omgz = omg(3:3:end);
        amx  = ang_mom(1:3:end); amy  = ang_mom(2:3:end); amz  = ang_mom(3:3:end);
        aax  = ang_acc(1:3:end); aay  = ang_acc(2:3:end); aaz  = ang_acc(3:3:end);
        
        % parfor 透明性要求: 在条件分支外初始化
        inx = []; iny = []; inz = [];
        if is_3d_local
            inx = inertia(1:3:end); iny = inertia(2:3:end); inz = inertia(3:3:end);
        end
        
        % 直接用缓存的 idx
        for t = 1:n_tags
            idx = idx_cache{t};
            
            if isempty(idx)
                continue;
            end
            
            ts_arr(step_idx, t)    = stepp_local(step_idx);
            omg_cells{step_idx, t} = [omgx(idx), omgy(idx), omgz(idx)];
            am_cells{step_idx, t}  = [amx(idx),  amy(idx),  amz(idx)];
            aa_cells{step_idx, t}  = [aax(idx),  aay(idx),  aaz(idx)];
            
            if is_3d_local
                in_cells{step_idx, t} = [inx(idx), iny(idx), inz(idx)];
            else
                in_cells{step_idx, t} = inertia(idx);
            end
        end
    end
    
    t_elapsed = toc;
    fprintf('  ✔ 转动数据提取完成! 耗时: %.2f 秒 (平均 %.3f 秒/文件)\n', ...
            t_elapsed, t_elapsed/num_steps);
    
    % --- 4.5 组装回与 v1 一致的 struct ---
    ExtractedData = struct();
    for t = 1:n_tags
        tag = tags_to_extract(t);
        t_name = sprintf('Tag_%d', abs(tag));
        ExtractedData.(t_name).TimeStep = ts_arr(:, t);
        ExtractedData.(t_name).Omg      = omg_cells(:, t);
        ExtractedData.(t_name).AngMom   = am_cells(:, t);
        ExtractedData.(t_name).AngAcc   = aa_cells(:, t);
        ExtractedData.(t_name).Inertia  = in_cells(:, t);
    end
    
    % --- 4.6 保存 ---
    tag_str = sprintf('%d_', abs(tags_to_extract));
    tag_str = tag_str(1:end-1);
    time_str = char(datetime('now', 'Format', 'yyyyMMdd_HHmmss'));
    save_name = sprintf('RotData_%s_Tag%s_Step%d_to_%d_%s.mat', ...
                        folder_name, tag_str, start_step, calc_end, time_str);
    save_path = fullfile(output_folder, save_name);
    
    fprintf('  正在打包保存数据到: %s ...\n', save_path);
    tic;
    save(save_path, 'ExtractedData', 'stepp', '-v7.3'); 
    fprintf('  ✔ 保存成功! 耗时: %.2f 秒\n', toc);
end

fprintf('\n✅ 所有匹配的文件夹【转动数据】均已处理完毕!\n');
fprintf('   数据已集中保存在 "%s" 目录下。\n', output_folder);