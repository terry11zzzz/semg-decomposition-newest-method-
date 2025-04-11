function dftFilterGUI
% dftFilterGUI - 一个GUI，用于配置ft_preproc_dftfilter的参数，对输入数据进行滤波，
% 并可生成滤波前后信号的时域图和频谱图，同时支持选择要显示的通道。
%
% 使用说明：
% 1. 在 MATLAB 工作区加载你的数据（例如变量名 dat），数据格式为：通道数 x 时间点。
% 2. 运行此函数，界面中会显示各参数设置项和几个按钮。
% 3. 根据需要修改采样率、噪声频率、滤波方法及相关参数，并指定显示通道（默认1）。
% 4. 点击“运行滤波”按钮滤波处理，然后点击“生成时域图”或“生成频谱图”
%    分别生成新窗口显示滤波前后的信号图形。
%
% 作者：你的朋友
% 日期：2025-02-16

%% 全局变量（在本函数内共享），存储滤波前后的数据及参数
origData = [];
filtData = [];
FsUsed   = [];
tAxis    = [];
channelDisplayed = 1; % 默认显示第1个通道

%% 创建主界面
hFig = figure('Name', 'DFT Filter GUI', 'NumberTitle', 'off', ...
    'Position', [100 100 900 600]);

%% 参数设置面板（左侧）
hPanel = uipanel('Parent', hFig, 'Title', '参数设置', ...
    'Position', [0.01 0.01 0.3 0.98]);

% 数据变量名
uicontrol('Parent', hPanel, 'Style', 'text', 'String', '数据变量名:', ...
    'HorizontalAlignment', 'left', 'Position', [10 560 120 20]);
hDataName = uicontrol('Parent', hPanel, 'Style', 'edit', 'String', 'dat', ...
    'Position', [140 560 120 25]);

% 采样率 Fs (Hz)
uicontrol('Parent', hPanel, 'Style', 'text', 'String', '采样率 Fs (Hz):', ...
    'HorizontalAlignment', 'left', 'Position', [10 520 120 20]);
hFs = uicontrol('Parent', hPanel, 'Style', 'edit', 'String', '1000', ...
    'Position', [140 520 120 25]);

% 电源噪声频率 Fl (Hz)
uicontrol('Parent', hPanel, 'Style', 'text', 'String', '电源噪声频率 Fl (Hz):', ...
    'HorizontalAlignment', 'left', 'Position', [10 480 120 20]);
hFl = uicontrol('Parent', hPanel, 'Style', 'edit', 'String', '50', ...
    'Position', [140 480 120 25]);

% 滤波方法 dftreplace
uicontrol('Parent', hPanel, 'Style', 'text', 'String', '滤波方法 dftreplace:', ...
    'HorizontalAlignment', 'left', 'Position', [10 440 120 20]);
hReplace = uicontrol('Parent', hPanel, 'Style', 'popupmenu', ...
    'String', {'zero','neighbour','neighbour_fft'}, ...
    'Position', [140 440 120 25]);

% dftbandwidth
uicontrol('Parent', hPanel, 'Style', 'text', 'String', 'dftbandwidth:', ...
    'HorizontalAlignment', 'left', 'Position', [10 400 120 20]);
hBandwidth = uicontrol('Parent', hPanel, 'Style', 'edit', 'String', '1', ...
    'Position', [140 400 120 25]);

% dftneighbourwidth
uicontrol('Parent', hPanel, 'Style', 'text', 'String', 'dftneighbourwidth:', ...
    'HorizontalAlignment', 'left', 'Position', [10 360 120 20]);
hNeighbour = uicontrol('Parent', hPanel, 'Style', 'edit', 'String', '2', ...
    'Position', [140 360 120 25]);

% 显示通道（默认1）
uicontrol('Parent', hPanel, 'Style', 'text', 'String', '显示通道:', ...
    'HorizontalAlignment', 'left', 'Position', [10 320 120 20]);
hChannel = uicontrol('Parent', hPanel, 'Style', 'edit', 'String', '1', ...
    'Position', [140 320 120 25]);

% 运行滤波按钮
hButton = uicontrol('Parent', hPanel, 'Style', 'pushbutton', 'String', '运行滤波', ...
    'FontSize', 12, 'Position', [40 260 220 40], 'Callback', @runFilter);

% 生成时域图按钮
hTimeButton = uicontrol('Parent', hPanel, 'Style', 'pushbutton', 'String', '生成时域图', ...
    'FontSize', 12, 'Position', [40 200 220 40], 'Callback', @plotTimeDomain);

% 生成频谱图按钮
hSpectrumButton = uicontrol('Parent', hPanel, 'Style', 'pushbutton', 'String', '生成频谱图', ...
    'FontSize', 12, 'Position', [40 140 220 40], 'Callback', @plotSpectrum);

%% 显示面板（右侧）：预留两个坐标区显示时域波形（滤波前与滤波后）
% 为避免重叠，这里分配了足够的空间
hAx1 = axes('Parent', hFig, 'Position', [0.35 0.55 0.6 0.35]);
title(hAx1, '原始信号（时域）');
xlabel(hAx1, '时间 (s)');
ylabel(hAx1, '幅值');

hAx2 = axes('Parent', hFig, 'Position', [0.35 0.1 0.6 0.35]);
title(hAx2, '滤波后信号（时域）');
xlabel(hAx2, '时间 (s)');
ylabel(hAx2, '幅值');

%% 回调函数：运行滤波处理
    function runFilter(~, ~)
        % 从界面获取参数
        dataName = get(hDataName, 'String');
        try
            % 尝试从基础工作区获取数据
            dat = evalin('base', dataName);
        catch
            errordlg(['无法加载数据变量：', dataName],'数据错误');
            return;
        end
        
        Fs_val = str2double(get(hFs, 'String'));
        if isnan(Fs_val) || Fs_val <= 0
            errordlg('采样率必须为正数','参数错误');
            return;
        end
        
        Fl_val = str2double(get(hFl, 'String'));
        if isnan(Fl_val) || Fl_val <= 0
            errordlg('电源噪声频率必须为正数','参数错误');
            return;
        end
        
        methods = get(hReplace, 'String');
        methodIdx = get(hReplace, 'Value');
        dftreplace_val = methods{methodIdx};
        
        dftbandwidth_val = str2double(get(hBandwidth, 'String'));
        if isnan(dftbandwidth_val)
            errordlg('dftbandwidth 输入无效','参数错误');
            return;
        end
        
        dftneighbourwidth_val = str2double(get(hNeighbour, 'String'));
        if isnan(dftneighbourwidth_val)
            errordlg('dftneighbourwidth 输入无效','参数错误');
            return;
        end
        
        % 获取显示通道号
        channelDisplayed = round(str2double(get(hChannel, 'String')));
        if isnan(channelDisplayed) || channelDisplayed < 1
            errordlg('显示通道必须为正整数','参数错误');
            return;
        end
        
        % 调用 ft_preproc_dftfilter 进行滤波处理
        try
            filt_data = ft_preproc_dftfilter(dat, Fs_val, Fl_val, ...
                'dftreplace', dftreplace_val, ...
                'dftbandwidth', dftbandwidth_val, ...
                'dftneighbourwidth', dftneighbourwidth_val);
        catch ME
            errordlg(['滤波过程中出错: ', ME.message],'运行错误');
            return;
        end
        
        % 保存原始数据和滤波后的数据，方便后续图形显示
        origData = dat;
        filtData = filt_data;
        FsUsed   = Fs_val;
        tAxis    = (0:size(dat,2)-1)/Fs_val;
        
        % 检查显示通道号是否超过数据通道数
        if channelDisplayed > size(origData, 1)
            errordlg('显示通道超出数据的通道数','参数错误');
            return;
        end
        
        % 更新右侧面板显示时域图（按选定的通道）
        axes(hAx1);
        plot(tAxis, origData(channelDisplayed, :), 'b', 'LineWidth', 1.5);
        title(['原始信号（时域）- 通道 ', num2str(channelDisplayed)]);
        xlabel('时间 (s)'); ylabel('幅值');
        
        axes(hAx2);
        plot(tAxis, filtData(channelDisplayed, :), 'r', 'LineWidth', 1.5);
        title(['滤波后信号（时域）- 通道 ', num2str(channelDisplayed)]);
        xlabel('时间 (s)'); ylabel('幅值');
        
        disp('滤波完成！');
    end

%% 回调函数：生成时域图（新窗口显示滤波前后信号）
    function plotTimeDomain(~, ~)
        if isempty(origData) || isempty(filtData)
            errordlg('请先运行滤波处理','提示');
            return;
        end
        % 检查显示通道号
        if channelDisplayed > size(origData,1)
            errordlg('显示通道超出数据的通道数','参数错误');
            return;
        end
        figure('Name','时域图','NumberTitle','off','Position',[200 200 800 600]);
        subplot(2,1,1);
        plot(tAxis, origData(channelDisplayed, :), 'b', 'LineWidth', 1.5);
        title(['原始信号（时域）- 通道 ', num2str(channelDisplayed)]);
        xlabel('时间 (s)'); ylabel('幅值');
        grid on;
        
        subplot(2,1,2);
        plot(tAxis, filtData(channelDisplayed, :), 'r', 'LineWidth', 1.5);
        title(['滤波后信号（时域）- 通道 ', num2str(channelDisplayed)]);
        xlabel('时间 (s)'); ylabel('幅值');
        grid on;
    end

%% 回调函数：生成频谱图（新窗口显示滤波前后信号的频谱）
    function plotSpectrum(~, ~)
        if isempty(origData) || isempty(filtData)
            errordlg('请先运行滤波处理','提示');
            return;
        end
        
        % 检查显示通道号
        if channelDisplayed > size(origData,1)
            errordlg('显示通道超出数据的通道数','参数错误');
            return;
        end
        
        % 计算FFT（以选定通道为例）
        N = length(tAxis);
        f = (0:N-1)*(FsUsed/N);  % 频率轴
        
        fft_orig = fft(origData(channelDisplayed, :));
        fft_filt = fft(filtData(channelDisplayed, :));
        
        % 归一化幅值
        mag_orig = abs(fft_orig)/N;
        mag_filt = abs(fft_filt)/N;
        
        % 新建窗口绘制频谱图
        figure('Name','频谱图','NumberTitle','off','Position',[250 250 800 600]);
        subplot(2,1,1);
        plot(f, mag_orig, 'b', 'LineWidth', 1.5);
        xlim([0, 100]);
        title(['原始信号频谱 - 通道 ', num2str(channelDisplayed)]);
        xlabel('频率 (Hz)'); ylabel('幅值');
        grid on;
        
        subplot(2,1,2);
        plot(f, mag_filt, 'r', 'LineWidth', 1.5);
        xlim([0, 100]);
        title(['滤波后信号频谱 - 通道 ', num2str(channelDisplayed)]);
        xlabel('频率 (Hz)'); ylabel('幅值');
        grid on;
    end

end
