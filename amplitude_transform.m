function transformed_signal = amplitude_transform(sEMG, method)
    % sEMG: 输入为 64 x 2000 的矩阵（通道 x 时间）
    % method: 'linear' 或 'sn'

    switch lower(method)
        case 'linear'
            transformed_signal = linear_amplitude_transform(sEMG);
        case 'sn'
            transformed_signal = segment_and_normalize(sEMG, 30); % 你可以调节分段数
        otherwise
            error('Unknown method. Use ''linear'' or ''sn''.');
    end
end

%% 方法1：线性振幅变换
function output = linear_amplitude_transform(sEMG)
    [nChannels, nSamples] = size(sEMG);
    output = zeros(size(sEMG));
    t = linspace(0, 1, nSamples); % 模拟时间轴

    for ch = 1:nChannels
        sig = sEMG(ch, :);

        m = round(nSamples * 0.1); % 取前10%点
        k = round(nSamples * 0.1); % 取后10%点

        Vm = mean(abs(sig(1:m)));
        Vk = mean(abs(sig(end-k+1:end)));

        y1 = Vk / Vm;
        y_end = 1;

        % 拟合线性函数 y = A*t + B
        A = (y_end - y1) / (t(end) - t(1));
        B = y1 - A * t(1);
        Y = A * t + B;

        output(ch, :) = sig .* Y;
    end
end

%% 方法2：分段归一化（S&N）
function output = segment_and_normalize(sEMG, numSegments)
    [nChannels, nSamples] = size(sEMG);
    output = zeros(size(sEMG));

    segmentLength = floor(nSamples / numSegments);

    for ch = 1:nChannels
        sig = sEMG(ch, :);
        normalized_sig = [];

        for i = 1:numSegments
            startIdx = (i-1)*segmentLength + 1;
            endIdx = min(i*segmentLength, nSamples);

            segment = sig(startIdx:endIdx);
            max_val = max(abs(segment));
            if max_val == 0
                norm_segment = segment; % 避免除以0
            else
                norm_segment = segment / max_val;
            end

            normalized_sig = [normalized_sig, norm_segment];
        end

        % 补上丢失的点（如果末尾没对齐）
        if length(normalized_sig) < nSamples
            normalized_sig = [normalized_sig, sig(length(normalized_sig)+1:end)];
        end

        output(ch, :) = normalized_sig;
    end
end
