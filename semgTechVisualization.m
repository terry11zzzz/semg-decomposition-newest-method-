function [finalres, union_set, tjns1set] = semgTechVisualization(semg, start_pt, end_pt, sample_window, K, para,amplitude_transform_method)
% semgTechVisualization - 对输入的 sEMG 信号进行预处理、白化和脉冲检测，并进行可视化
%
% 语法:
%   [pulsefinal, union_set, tjns1set] = semgTechVisualization(semg, start_pt, end_pt)
%   [pulsefinal, union_set, tjns1set] = semgTechVisualization(semg, start_pt, end_pt, sample_window, K, para)
%
% 输入:
%   semg         - 数值矩阵，尺寸为 (通道数 x 样本点数)，如果输入不是数值矩阵则报错
%   start_pt     - 起始样本点（标量）
%   end_pt       - 结束样本点（标量）
%   sample_window- 初始分解样本点数，默认 15000
%   K            - 每个通道扩展倍数，默认 13
%   para         - 参数向量，默认值为 [32, 16, 28, 40]-sample_window=15000
%   amplitude_transform_method - 幅值变换方法，默认不变换，线性变换就是方法1，分段变换就是方法2
%
% 输出:
%   pulsefinal   - 结构体数组，包含检测到的脉冲信息（每个元素对应一个 motor unit）
%   union_set    - (来自 dynamic_decp_whiten_fun 的输出)
%   tjns1set     - (来自 dynamic_decp_whiten_fun 的输出)
%
% 注意: 
%   此函数依赖函数 dynamic_decp_whiten_fun，请确保该函数在 MATLAB 路径中。
%
% 小贴士: 如果输入格式不对，函数会直接报错哦～
%
% 作者: yun zheng
% 日期: 02/16/2025

    % 检查输入 semg 是否为数值矩阵
    if ~isnumeric(semg) || iscell(semg)
        error('Input semg signal must be a numeric matrix (channels x samples).');
    end

    % 参数默认值设定
    if nargin < 4 || isempty(sample_window)
        sample_window = 15000;
    end
    if nargin < 5 || isempty(K)
        K = 13;
    end
    if nargin < 6 || isempty(para)
        para = [32, 16, 28, 40];
    end
    if nargin < 7 || isempty(amplitude_transform_method)
        amplitude_transform_method=0;
    end
    % 从输入的 semg 信号中截取需要的部分
    if end_pt > size(semg,2) || start_pt < 1 || start_pt >= end_pt
        error('start_pt 和 end_pt 设置有误，请检查！');
    end
    signal = semg(:, start_pt:end_pt);
    if strcmp(amplitude_transform_method, 'linear')
        signal = amplitude_transform(signal, 'linear');
    elseif strcmp(amplitude_transform_method, 'sn')
        signal = amplitude_transform(signal, 'sn');
    end
    % 获取基本参数
    [M, nsamples] = size(signal);
    
    % 构建扩展的脉冲矩阵 xn_exten
    xn_exten = zeros(M*K, nsamples);
    for i = 1:M
        for j = 1:K
            % 每个通道的扩展版本
            xn_exten(K*(i-1)+j, :) = [zeros(1, j-1), signal(i, 1:end-j+1)];
        end
    end
    
    % 去均值操作
    for i = 1:size(xn_exten, 1)
        xn_exten(i,:) = xn_exten(i,:) - mean(xn_exten(i,:));
    end
    
    % 为了避免边缘效应，从第101个样本开始白化
    if size(xn_exten,2) < 101
        error('扩展后的信号样本点太少，无法进行白化。');
    end
    zn_exten = Whiten(xn_exten(:, 101:end));
    
    % 设置 gap 为 sample_window 的一半
    gap = sample_window / 4;
    
    % 检查白化后信号的样本点数是否足够
    if size(zn_exten, 2) < sample_window
        error('白化后的信号样本点数不足，无法达到指定的 sample_window。');
    end
    
    % 调用动态分解函数进行初步分解
    % 注意：请确保 dynamic_decp_whiten_fun 已经在 MATLAB 路径中！
    [pulsefinal, ~, union_set, tjns1set] = dynamic_decp_whiten_fun(zn_exten(:, 1:sample_window), gap, para);
        
    % ---------去除空数组--------------
    % 找出所有pulse字段不为空的元素的索引
    notEmptyIndices = arrayfun(@(x) ~isempty(x.pulse), pulsefinal);
    % 保留那些pulse字段不为空的元素
    pulsefinal = pulsefinal(notEmptyIndices);
    
    %% STEP4: precision_pro 处理
    pulsefinal1 = pulsefinal;
    % 这里假设为 static 数据，sp = start_pt+100
    sp = start_pt + 100;
    ep = sample_window;
    save_all_ctj = zn_exten' * zn_exten;
    for loop2 = 1:length(pulsefinal1)
        if isempty(pulsefinal1(loop2).pulse)
            continue;
        end
        temp_dif = diff(pulsefinal1(loop2).pulse);
        pulsefinal1(loop2).cov = std(temp_dif) / mean(temp_dif);
        pulsefinal1(loop2).ave_gap = mean(temp_dif);
        pulsefinal1(loop2).pulse(pulsefinal1(loop2).pulse < 1) = [];
        pulsefinal1(loop2).pulse(pulsefinal1(loop2).pulse > ep) = [];
        pulsefinal1(loop2).ctsj = mean(save_all_ctj(pulsefinal1(loop2).pulse, :), 1);
        pulsefinal1(loop2).init_pulsetime_num = length(pulsefinal1(loop2).pulse);
    end
    % 删除重复出现次数较少的序列（unionpulsetime 字段由 dynamic_decp_whiten_fun 产生）
    v = 0;
    while 1
        v = v + 1;
        if v > length(pulsefinal1)
            break;
        end
        if ~isfield(pulsefinal1(v), 'unionpulsetime')
            continue;
        end
        if length(pulsefinal1(v).unionpulsetime) <= 1.8 * length(pulsefinal1(v).pulse)
           pulsefinal1(v) = [];
           v = v - 1;
        end
        if v >= length(pulsefinal1)
            break;
        end
    end
    
    % 调用 CalMUAPNonShift 计算 MUAPSet（请确保此函数在路径中）
    [MUAPSet] = CalMUAPNonShift(ep, signal(:,101:end), pulsefinal1, [], size(signal,1));
    
    %% STEP5: 人工筛选 —— 图形界面依次展示每个 MU 的信息
    decisionArr = false(1, length(pulsefinal1));
    for idx = 1:length(pulsefinal1)
        % 调用界面函数，展示：subplot(2,1,1)显示当前 MU 的 MUAP 曲线，subplot(2,1,2)显示 pulsefinal1.ctsj，
        % 并提供“保留”与“删除”按钮供你选择。
        decision = getUserDecisionForPulse(idx,MUAPSet, pulsefinal1(idx).ctsj);
        decisionArr(idx) = decision;
    end
    % 删除用户选择删除的 MU
    pulsefinal1 = pulsefinal1(decisionArr);
    
    fprintf('精细处理和人工筛选完成，最终保留 %d 个脉冲序列。\n', length(pulsefinal1));

    %% STEP6: 自动提取时刻
    %加入gap向量

    finalres = Dynamic_for_real_stric_MUAP1_visual(semg, start_pt, end_pt, pulsefinal1, []);
    finalres = Dynamic_for_real_stric_MUAP2_visual(semg, start_pt, end_pt, finalres);
    % 调用 CalMUAPNonShift 计算 MUAPSet（请确保此函数在路径中）
    [MUAPSet] = CalMUAPNonShift(ep, signal(:,101:end), finalres, [], size(signal,1));
    
    %% STEP7: 人工筛选 —— 图形界面依次展示每个 MU 的信息
    decisionArr = false(1, length(finalres));
    for idx = 1:length(finalres)
        % 调用界面函数，展示：subplot(2,1,1)显示当前 MU 的 MUAP 曲线，subplot(2,1,2)显示 pulsefinal1.ctsj，
        % 并提供“保留”与“删除”按钮供你选择。
        decision = getUserDecisionForPulse(idx,MUAPSet, finalres(idx).ctsj);
        decisionArr(idx) = decision;
    end
    % 删除用户选择删除的 MU
    finalres = finalres(decisionArr);
    
    fprintf('精细处理和人工筛选完成，最终保留 %d 个脉冲序列。\n', length(finalres));
    
    %% 嵌套函数：图形界面供用户选择保留/删除
    function decision = getUserDecisionForPulse(idx,MUAPSet, ctsj)
        % 创建界面，展示当前 MU 的两幅图：第一幅为 MUAP 曲线，第二幅为 ctsj 曲线
        hFig = figure('Position', [100, 100, 800, 1200], 'Name', '人工筛选界面');
        
        subplot(2,1,1);
        % 对于每个 pulsefinal1，联合各通道的 MUAP 集合
        current_allmuap = [];
        for jjj = 1:size(signal,1)
            current_allmuap = union(current_allmuap, MUAPSet{jjj, idx});
            plot(MUAPSet{jjj,idx});hold on
        end
        
        title('当前 MU 的 MUAP 图');
        xlabel('样本点'); ylabel('幅值'); grid on;
        
        subplot(2,1,2);
        plot(ctsj, 'b');
        title('当前 MU 的 ctsj 图');
        xlabel('样本点'); ylabel('ctsj 值'); grid on;
        
        % 添加按钮：保留 和 删除
        uicontrol('Style', 'pushbutton', 'String', '保留', ...
            'Position', [150, 20, 100, 40], 'Callback', @(src,event) buttonCallback(true));
        uicontrol('Style', 'pushbutton', 'String', '删除', ...
            'Position', [300, 20, 100, 40], 'Callback', @(src,event) buttonCallback(false));
        
        % 等待用户作出选择
        uiwait(hFig);
        decision = getappdata(hFig, 'decision');
        if ishandle(hFig)
            close(hFig);
        end
        
        function buttonCallback(dec)
            setappdata(hFig, 'decision', dec);
            uiresume(hFig);
        end
    end



    
end

% --- 子函数：白化 ---
function zn = Whiten(x)
% Whiten - 对矩阵 x 进行白化处理
%
% 输入:
%   x - 输入矩阵，每一行代表一个信号
%
% 输出:
%   zn - 白化后的矩阵

    % 再次去均值，防止遗漏
    for i = 1:size(x, 1)
         x(i,:) = x(i,:) - mean(x(i,:));
    end
    % 计算协方差矩阵
    C = x*x' / size(x,2);
    [u,s,~] = svd(C);
    % 计算白化矩阵 W
    W = u * diag(1./sqrt(diag(s))) * u';
    zn = W * x;
end
