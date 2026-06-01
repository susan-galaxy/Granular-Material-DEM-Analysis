% =====================================================================
%  自动化 HDF5 颗粒数据提取脚本 (仅位置、速度、半径)
%  支持通配符批量匹配文件夹，并统一输出到指定目录
% =====================================================================

clear all;
close all;
clc;

%% ========== 1. 核心配置区 ==========

% 需要处理的文件夹名称模式 (支持通配符)
% 修改为匹配 FixSpace 开头的文件夹
folder_pattern = 'En0.94*'; 

% 提取出的数据 (.mat文件) 统一保存的输出文件夹路径
% 如果该文件夹不存在，脚本会自动创建
output_folder = '/media/gezhuan/M78/Space_Active/Data/Maxwell_Boltz/data/Pos_Vel_data/En0_94';

% 时间步设置
start_step = 1;       % 起始步数
end_step   = 3000;   % 结束步数
delta_step = 1;       % 步长

% 需要提取的颗粒类型 (Tag)
% -1: 活性颗粒, -3: 亚克力颗粒, -4: 硅胶颗粒
tags_to_extract = -2;

%% ========== 2. 文件夹检索与环境准备 ==========

% 检索符合命名模式的项目
items = dir(folder_pattern);
target_folders = {};

% 筛选出真正的文件夹（排除文件以及 '.' 和 '..'）
for k = 1:length(items)
    if items(k).isdir && ~strcmp(items(k).name, '.') && ~strcmp(items(k).name, '..')
        % 记录文件夹的完整绝对路径或相对路径
        target_folders{end+1} = fullfile(items(k).folder, items(k).name);
    end
end

if isempty(target_folders)
    error('未找到匹配模式 "%s" 的文件夹，请检查路径或拼写！', folder_pattern);
end

fprintf('成功匹配到 %d 个文件夹，准备开始提取数据...\n', length(target_folders));

% 确保输出文件夹存在
if ~exist(output_folder, 'dir')
    mkdir(output_folder);
    fprintf('已创建数据输出目录: %s\n', output_folder);
end

%% ========== 3. 自动化批处理主循环 ==========

for f_idx = 1:length(target_folders)
    current_folder = target_folders{f_idx};
    
    % 【关键修复】分别提取名字和“后缀”，然后重新拼接，防止小数点后面的数字丢失
    [~, name_part, ext_part] = fileparts(current_folder);
    folder_name = [name_part, ext_part]; 
    
    fprintf('\n======================================================\n');
    fprintf('▶ 正在处理 [%d/%d]: %s\n', f_idx, length(target_folders), folder_name);
    
    % --- 3.1 获取并筛选 H5 文件 ---
    file_pattern = fullfile(current_folder, 'Space_s_*.h5');
    all_files = dir(file_pattern);
    % 用 logical indexing 一次性过滤掉包含 '_bf_' 的边界力文件，避免动态数组增长
    if isempty(all_files)
        df_normal = all_files;
    else
        keep_mask = ~contains({all_files.name}, '_bf_');
        df_normal = all_files(keep_mask);
    end
    
    total_files = length(df_normal);
    if total_files == 0
        warning('文件夹 %s 中未找到符合条件的 .h5 文件，跳过此文件夹。', folder_name);
        continue;
    end
    
    % 确定实际要处理的步数范围
    calc_end = min(end_step, total_files);
    stepp = start_step:delta_step:calc_end;
    num_steps = length(stepp);
    fprintf('  找到 %d 个文件，计划处理 %d 个时间步 (从 %d 到 %d)。\n', total_files, num_steps, start_step, calc_end);
    
    % --- 3.2 初始化存储结构体 ---
    % 结构: ExtractedData.Tag_1.Pos{step_idx} = [N x 3] 矩阵
    % 预先计算所有 tag 对应的字段名，避免在内循环中反复调用 sprintf
    num_tags = length(tags_to_extract);
    tag_names = cell(num_tags, 1);
    for t = 1:num_tags
        tag_names{t} = sprintf('Tag_%d', abs(tags_to_extract(t)));
    end
    
    ExtractedData = struct();
    for t = 1:num_tags
        t_name = tag_names{t};
        ExtractedData.(t_name).TimeStep = zeros(num_steps, 1);
        ExtractedData.(t_name).Pos      = cell(num_steps, 1);
        ExtractedData.(t_name).Vel      = cell(num_steps, 1);
        ExtractedData.(t_name).Radius   = cell(num_steps, 1);
    end
    
    % --- 3.3 开始逐个文件提取数据 ---
    tic;
    
    % parfor 不支持向 ExtractedData.(t_name).Field{step_idx} 这种嵌套字段做切片写入，
    % 所以先用一个真正可切片的 cell 数组 step_data 收集每步结果，循环外再串行合并。
    step_data = cell(num_steps, 1);
    
    parfor step_idx = 1:num_steps
        i = stepp(step_idx);
        filename = fullfile(current_folder, df_normal(i).name);
        
        % [读取基础数据]
        pos  = h5read(filename, '/Position');
        vel  = h5read(filename, '/PVelocity');
        ptag = h5read(filename, '/PTag');
        rad  = h5read(filename, '/Radius');
        
        % 用 reshape 把 1D 数据视图重排为 3xN（O(1) 操作，避免三次 stride 索引）
        pos3 = reshape(pos, 3, []);
        vel3 = reshape(vel, 3, []);
        
        % 把本步每个 Tag 的结果打包成一个 1 x num_tags 的本地 cell（空表示该步无此 Tag）
        local = cell(1, num_tags);
        for t = 1:num_tags
            tag = tags_to_extract(t);
            
            % 用 logical mask 一次定位（替代 find，节省一次中间向量分配）
            mask = (ptag == tag);
            if ~any(mask)
                continue; % 当前步没有这个Tag的颗粒，跳过
            end
            
            s = struct();
            s.TimeStep = i;
            s.Pos      = pos3(:, mask).';
            s.Vel      = vel3(:, mask).';
            s.Radius   = rad(mask);
            local{t}   = s;
        end
        step_data{step_idx} = local;
        
        % 打印进度（parfor 下顺序可能乱，但仍能看到推进）
        if mod(step_idx, 100) == 0 || step_idx == num_steps
            fprintf('  进度: %d / %d (%.1f%%)\n', step_idx, num_steps, (step_idx/num_steps)*100);
        end
    end
    
    % 串行合并到 ExtractedData（嵌套字段写入只能在串行 for 里，一次性遍历开销很小）
    for step_idx = 1:num_steps
        local = step_data{step_idx};
        if isempty(local), continue; end
        for t = 1:num_tags
            s = local{t};
            if isempty(s), continue; end
            t_name = tag_names{t};
            ExtractedData.(t_name).TimeStep(step_idx) = s.TimeStep;
            ExtractedData.(t_name).Pos{step_idx}      = s.Pos;
            ExtractedData.(t_name).Vel{step_idx}      = s.Vel;
            ExtractedData.(t_name).Radius{step_idx}   = s.Radius;
        end
        step_data{step_idx} = []; % 边合并边释放内存
    end
    clear step_data;
    
    t_elapsed = toc;
    fprintf('  数据提取完成! 耗时: %.2f 秒\n', t_elapsed);
    
    % --- 3.4 保存当前文件夹的数据为 .mat 文件 ---
    
    % 1. 将本次提取的 Tag 信息拼接成字符串
    tag_str = sprintf('%d_', abs(tags_to_extract));
    tag_str = tag_str(1:end-1); % 移除末尾下划线
    
    % 2. 获取当前时间戳防覆盖
    time_str = char(datetime('now', 'Format', 'yyyyMMdd_HHmmss'));
    
    % 3. 拼接全新文件名 (包含文件夹全名、Tag、步数范围、时间戳)
    save_name = sprintf('Data_%s_Tag%s_Step%d_to_%d_%s.mat', folder_name, tag_str, start_step, calc_end, time_str);
    save_path = fullfile(output_folder, save_name);
    
    fprintf('  正在打包保存数据到: %s ...\n', save_path);
    % 采用 v7.3 格式保存以支持大于 2GB 的大型变量
    save(save_path, 'ExtractedData', 'stepp', '-v7.3'); 
    fprintf('  保存成功!\n');
end

fprintf('\n✅ 所有匹配的文件夹均已处理完毕！数据已集中保存在 "%s" 目录下。\n', output_folder);