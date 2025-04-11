function [xreal,xreal2,Pnr_con] = find_middle_2time(precise_init_MU,save_all_ctj,predict_gap,signal,large_mu_time,zn_exten15wgn,left_time,right_time,window_len)
%rough_time 使用准确初始时刻得到粗略的放电时刻
%   目的得到 all good time and some bad time
%   使用线性加权权向量
    amp_coeff=0.3;
    Pnr_con=0;
    k_gap=0.3;
    %% 左边时刻 left time
    ixvn1set=precise_init_MU;
    c=left_time+predict_gap;
    delta_a=min(c-1,window_len-2000);
    delta_b=min(precise_init_MU(end)-c,window_len-2000);
    hat_t_j=mean(save_all_ctj(ixvn1set,:),1);
    n_c = c - delta_a : c;
    f_values1 = hat_t_j(n_c) .* (0.7 - (c - delta_a - n_c) / (2*delta_a));
    n_c = c : c+delta_b ;
    f_values2 = hat_t_j(n_c) .* (0.7 + (c + delta_b - n_c) / (2*delta_b));
    f_values=[f_values1,f_values2(2:end)];
    hat_c_t_j_y = Eq10(c - delta_a :c+delta_b, f_values,zn_exten15wgn);
    Cstj=hat_c_t_j_y'*zn_exten15wgn;
    %plot(Cstj)
    xreal=cell(1,1);
    % 筛选时刻范围内幅值显著大于其他时刻的直接输出
    start_p=round(left_time+(1-k_gap)*predict_gap);
    end_p=round(left_time+(1+k_gap)*predict_gap);
    [c,ind]=max( Cstj(start_p:end_p));%从均值+3std到均值-3std找最大值
    ave_amp=mean( Cstj(start_p:end_p));%得到平均幅值
    std_amp=std( Cstj(start_p:end_p));%得到平均标准差
    if c>ave_amp+amp_coeff*std_amp
        xreal{1}=start_p+ind-1;
        xreal2=[];
        Pnr_con=100;
        return;
    end
    % IPT峰值时刻筛选  
    [PKS,LOCS]= findpeaks(Cstj,'MinPeakDistance',3, 'MinPeakHeight', 0);   
    for lp1=1:length(ixvn1set)
        LOCS(LOCS>ixvn1set(lp1)-round(0.5*predict_gap) & LOCS<ixvn1set(lp1)+round(predict_gap*0.5))=[];
    end
    LOCS=union(LOCS,ixvn1set);    
    xreal{1}=LOCS(LOCS>=left_time+(1-k_gap)*predict_gap & LOCS<=left_time+(1+k_gap)*predict_gap);
    % NEO 算子
    peak_detec=abs(Cstj(:,2:window_len-1)).*abs(Cstj(:,2:window_len-1))-abs(Cstj(:,1:window_len-2)).*abs(Cstj(:,3:window_len));
    max_values = max(peak_detec, [], 2);
    peak_detec=peak_detec./max_values;
    peak_detec=[0,peak_detec,0];
    [PKS,LOCS1]= findpeaks(peak_detec,'MinPeakDistance',3, 'MinPeakHeight', 0); 
    for lp1=1:length(ixvn1set)
        LOCS1(LOCS1>ixvn1set(lp1)-round(0.5*predict_gap) & LOCS1<ixvn1set(lp1)+round(predict_gap*0.5))=[];
    end
    LOCS1=union(LOCS1,ixvn1set);
    xreal{2}=LOCS1(LOCS1>=left_time+(1-k_gap)*predict_gap & LOCS1<=left_time+(1+k_gap)*predict_gap);
    binary_xreal=Cstj(xreal{2})>0;
    xreal{2}=xreal{2}(binary_xreal);
    xreal{1}=union(xreal{1},xreal{2});
    % 如果xreal{1}直接不存在，那么return
    if isempty(xreal{1})
        xreal=[];
        Pnr_con=0;
        xreal2=[];
        return;
    end

    % FCE筛选出绝对错误的时刻,要注意，FCE检测要求边界点不能小于25，大于4975,
    % 由于一般预留500左右时刻，因此不会到达4975边界点
    % error_pulse是加入新时刻的pulse
    FCE_RES=zeros(1,length(xreal{1}));
    COR_RES=zeros(1,length(xreal{1}));
    channel_num=round(size(signal,1)/2);
    for l1=1:length(xreal{1})
        error_pulse= xreal{1}(l1);
        temp_pulse1=ixvn1set(ixvn1set>25);%边界点不能小于25
        temp_pulse1=temp_pulse1(temp_pulse1<window_len-25);
        temp_pulse=union(temp_pulse1,error_pulse);
        temp_pulse(temp_pulse<1 & temp_pulse>window_len)=[];
        idx = find(temp_pulse == error_pulse);
        temp_pulse1(temp_pulse1<1 & temp_pulse1>window_len)=[];
        waveform_length=45;
        MU_waveforms = extractMUWaveforms(signal(channel_num,:), temp_pulse, waveform_length);
        MU_waveforms1 = extractMUWaveforms(signal(channel_num,:), temp_pulse1, waveform_length);
        % % e:噪声的均方根值 思路：先用neo算子去确定噪声时刻，再使用噪声时刻得到其均方根值 RMSE = sqrt(sum(noise_data.^2) / n)
        window_size = 15; 
        noise_time=find_noise(signal(channel_num,:),window_size,window_len);
        noise_data=signal(channel_num,noise_time);
        e = sqrt(sum(noise_data.^2) / length(noise_data));
        % 假设已经有了以下变量：
        % m: 选择的特征数
        % 初始化g[n]向量
        % 计算g[n]，即每个样本与其后一个样本的差值
        % yi 暂定为均值
        yi=mean(MU_waveforms1,1);
        xi=MU_waveforms(idx,:);
        g = abs(diff(xi));
        m=4;
        % 选择具有最大g[n]值的m个样本点的索引
        [~, sorted_idx] = sort(g, 'descend');
        selected_idx = sorted_idx(1:m);
        alpha=0.005;
        beta=0.05;
        % 获取选定的样本点作为特征
        features = xi(selected_idx);
        % 假设MUPT模板的MUP模板信号为yi(selected_idx)（大小为m的向量）
        % 计算形状不一致性百分比PSI_i
        diff1 = yi(selected_idx) - features ;
        term1=0;term2=0;
        for ii=1:length(diff1)
            term1 = term1+step_function(diff1(ii)- 3 * e);
            term2 = term2+step_function(-diff1(ii)- 3 * e);
        end
        PSI_i = (sum(term1) + sum(term2)) / m;
        % 计算d_i
        d_i = sum(diff1.^2) / e^2;
        % 显示结果
        %fprintf('MUP %d 的形状不一致性百分比 (PSI): %f\n', PSI_i);
        chi2_alpha = chi2inv(1 - alpha, m); % 这里假设 m=3
        if d_i > chi2_alpha && PSI_i > 0.7
            %fprintf('xi Definitely a FCE\n');
            FCE_RES(l1)=1;
        elseif d_i <chi2inv(1 - beta, m)
            COR_RES(l1)=1;
            %fprintf('Definitely a correct\n');
        end
        
    end
    % 筛选已知MU的时刻（用来排除）这个有点难搞，如果大于3个了，说明重叠时刻太多了，去掉一个
    is_largeMU=zeros(1,length(xreal{1}));
    if ~isempty(large_mu_time) 
        for l1=1:length(xreal{1})
            set1=union(ixvn1set,xreal{1}(l1));
            pu1 = zeros(window_len, 1);
            pu2 = zeros(window_len, 1);
            pu1(set1) = 1;
            for l2=1:length(large_mu_time)
                pulse=large_mu_time{l2};
                pulse(pulse<10 | pulse>window_len-100)=[];
                %pulse1 = bsxfun(@plus, pulse, -3:3);
                pulse1=[];
                for j = -3:3
                    pulse1 = [pulse1, j+pulse];
                end
                set2 = reshape(pulse1, 1, []);
                %set2 = repmat(pulse, 1, 7) + kron(-3:3, ones(1, length(pulse)));
                pu2(set2)=1;
                xcpul = xcorr(pu1, pu2, 15);
                [maxxc, ipmaxxc] = max(xcpul);
                % Remove the last occurrence of each common element if it appears more than 3 times
                if maxxc>=0.5* length(set1)
                    is_largeMU(l1)=1;
                    break;
                end
            end
        end
    end
    % 筛选那些 让 波形变化太大的时刻们
    MUAP_ratio=zeros(1,length(xreal{1}));
    for l1=1:length(xreal{1})
        finalres.pulse=union(ixvn1set,xreal{1}(l1));
        [MUAPSet]=CalMUAPNonShift(size(signal,2),signal,finalres,[],size(signal,1));
        close all
        for llll=1:length(finalres)
            tallmuap=[];
            for jjj=1:size(signal,1)
                tallmuap=[tallmuap,MUAPSet{jjj,llll}];
                %plot(MUAPSet{jjj,llll});hold on
            end
            gap=max(tallmuap)-min(tallmuap);
            [~, indices] = maxk(abs(tallmuap), 100);
            % 将这些数从数组中移除
            tallmuap(indices) = [];
            tstd=std(tallmuap);
            MUAP_ratio(l1)=gap/tstd;
        end
    end
    [~, indices] = maxk(MUAP_ratio, 15);
    % 创建一个与MUAP_ratio同样大小的全0矩阵
    top15Matrixl = ones(size(MUAP_ratio));
    % 将最大的5个数的位置置为1
    top15Matrixl(indices) = 0;
    %% 右边时刻 right time
    ixvn1set=precise_init_MU;
    c=right_time-predict_gap;
    delta_a=min(c-1,window_len-2000);
    delta_b=min(precise_init_MU(end)-c,window_len-2000);
    hat_t_j=mean(save_all_ctj(ixvn1set,:),1);
    n_c = c - delta_a : c;
    f_values1 = hat_t_j(n_c) .* (0.7 - (c - delta_a - n_c) / (2*delta_a));
    n_c = c : c+delta_b ;
    f_values2 = hat_t_j(n_c) .* (0.7 + (c + delta_b - n_c) / (2*delta_b));
    f_values=[f_values1,f_values2(2:end)];
    hat_c_t_j_y = Eq10(c - delta_a :c+delta_b, f_values,zn_exten15wgn);
    Cstj=hat_c_t_j_y'*zn_exten15wgn;
    %plot(Cstj)
    xreal2=cell(1,1);
    % 筛选时刻范围内幅值显著大于其他时刻的直接输出
    start_p=round(right_time-(1+k_gap)*predict_gap);
    end_p=round(right_time-(1-k_gap)*predict_gap);
    [c,ind]=max( Cstj(start_p:end_p));%从均值+3std到均值-3std找最大值
    ave_amp=mean( Cstj(start_p:end_p));%得到平均幅值
    std_amp=std( Cstj(start_p:end_p));%得到平均标准差
    if c>ave_amp+amp_coeff*std_amp% 如果xreal没有直接return，说明不太靠谱，没必要考虑
        xreal2{1}=start_p+ind-1;
        xreal=[];
        Pnr_con=100;
        return;
    end
    % IPT峰值时刻筛选  
    [PKS,LOCS]= findpeaks(Cstj,'MinPeakDistance',3,'MinPeakHeight', 0);   
    for lp1=1:length(ixvn1set)
        LOCS(LOCS>ixvn1set(lp1)-round(0.5*predict_gap) & LOCS<ixvn1set(lp1)+round(predict_gap*0.5))=[];
    end
    LOCS=union(LOCS,ixvn1set);    
    xreal2{1}=LOCS(LOCS>=right_time-(1+k_gap)*predict_gap & LOCS<=right_time-(1-k_gap)*predict_gap);
    % NEO 算子
    peak_detec=abs(Cstj(:,2:window_len-1)).*abs(Cstj(:,2:window_len-1))-abs(Cstj(:,1:window_len-2)).*abs(Cstj(:,3:window_len));
    max_values = max(peak_detec, [], 2);
    peak_detec=peak_detec./max_values;
    peak_detec=[0,peak_detec,0];
    [PKS,LOCS1]= findpeaks(peak_detec,'MinPeakDistance',3, 'MinPeakHeight', 0); 
    for lp1=1:length(ixvn1set)
        LOCS1(LOCS1>ixvn1set(lp1)-round(0.5*predict_gap) & LOCS1<ixvn1set(lp1)+round(predict_gap*0.5))=[];
    end
    LOCS1=union(LOCS1,ixvn1set);
    xreal2{2}=LOCS1(LOCS1>=right_time-(1+k_gap)*predict_gap & LOCS1<=right_time-(1-k_gap)*predict_gap);
    binary_xreal=Cstj(xreal2{2})>0;
    xreal2{2}=xreal2{2}(binary_xreal);
    xreal2{1}=union(xreal2{1},xreal2{2});
    % 如果xreal2{1}直接不存在，那么return
    if isempty(xreal2{1})
        %xreal=[];
        xreal2=[];
        Pnr_con=0;
        return;
    end

    % FCE筛选出绝对错误的时刻,要注意，FCE检测要求边界点不能小于25，大于4975,
    % 由于一般预留500左右时刻，因此不会到达4975边界点
    % error_pulse是加入新时刻的pulse
    FCE_RES2=zeros(1,length(xreal2{1}));
    COR_RES=zeros(1,length(xreal2{1}));
    channel_num=round(size(signal,1)/2);
    for l1=1:length(xreal2{1})
        error_pulse= xreal2{1}(l1);
        temp_pulse1=ixvn1set(ixvn1set>25);%边界点不能小于25
        temp_pulse1=temp_pulse1(temp_pulse1<window_len-25);
        temp_pulse=union(temp_pulse1,error_pulse);
        temp_pulse(temp_pulse<1 & temp_pulse>window_len)=[];
        idx = find(temp_pulse == error_pulse);
        temp_pulse1(temp_pulse1<1 & temp_pulse1>window_len)=[];
        waveform_length=45;
        MU_waveforms = extractMUWaveforms(signal(channel_num,:), temp_pulse, waveform_length);
        MU_waveforms1 = extractMUWaveforms(signal(channel_num,:), temp_pulse1, waveform_length);
        % % e:噪声的均方根值 思路：先用neo算子去确定噪声时刻，再使用噪声时刻得到其均方根值 RMSE = sqrt(sum(noise_data.^2) / n)
        window_size = 15; 
        noise_time=find_noise(signal(channel_num,:),window_size,window_len);
        noise_data=signal(channel_num,noise_time);
        e = sqrt(sum(noise_data.^2) / length(noise_data));
        % 假设已经有了以下变量：
        % m: 选择的特征数
        % 初始化g[n]向量
        % 计算g[n]，即每个样本与其后一个样本的差值
        % yi 暂定为均值
        yi=mean(MU_waveforms1,1);
        xi=MU_waveforms(idx,:);
        g = abs(diff(xi));
        m=4;
        % 选择具有最大g[n]值的m个样本点的索引
        [~, sorted_idx] = sort(g, 'descend');
        selected_idx = sorted_idx(1:m);
        alpha=0.005;
        beta=0.05;
        % 获取选定的样本点作为特征
        features = xi(selected_idx);
        % 假设MUPT模板的MUP模板信号为yi(selected_idx)（大小为m的向量）
        % 计算形状不一致性百分比PSI_i
        diff1 = yi(selected_idx) - features ;
        term1=0;term2=0;
        for ii=1:length(diff1)
            term1 = term1+step_function(diff1(ii)- 3 * e);
            term2 = term2+step_function(-diff1(ii)- 3 * e);
        end
        PSI_i = (sum(term1) + sum(term2)) / m;
        % 计算d_i
        d_i = sum(diff1.^2) / e^2;
        % 显示结果
        %fprintf('MUP %d 的形状不一致性百分比 (PSI): %f\n', PSI_i);
        chi2_alpha = chi2inv(1 - alpha, m); % 这里假设 m=3
        if d_i > chi2_alpha && PSI_i > 0.7
            %fprintf('xi Definitely a FCE\n');
            FCE_RES2(l1)=1;
        elseif d_i <chi2inv(1 - beta, m)
            COR_RES(l1)=1;
            %fprintf('Definitely a correct\n');
        end
        
    end
    % 筛选已知MU的时刻（用来排除）这个有点难搞，如果大于3个了，说明重叠时刻太多了，去掉一个
    is_largeMU2=zeros(1,length(xreal2{1}));
    if ~isempty(large_mu_time) 
        for l1=1:length(xreal2{1})
            set1=union(ixvn1set,xreal2{1}(l1));
            pu1 = zeros(window_len, 1);
            pu2 = zeros(window_len, 1);
            pu1(set1) = 1;
            for l2=1:length(large_mu_time)
                pulse=large_mu_time{l2};
                pulse(pulse<10 | pulse>window_len-100)=[];
                %pulse1 = bsxfun(@plus, pulse, -3:3);
                pulse1=[];
                for j = -3:3
                    pulse1 = [pulse1, j+pulse];
                end
                set2 = reshape(pulse1, 1, []);
                %set2 = repmat(pulse, 1, 7) + kron(-3:3, ones(1, length(pulse)));
                pu2(set2)=1;
                xcpul = xcorr(pu1, pu2, 15);
                [maxxc, ipmaxxc] = max(xcpul);
                % Remove the last occurrence of each common element if it appears more than 3 times
                if maxxc>=0.5* length(set1)
                    is_largeMU2(l1)=1;
                    break;
                end
            end
        end
    end
    % 筛选那些 让 波形变化太大的时刻们
    MUAP_ratio=zeros(1,length(xreal2{1}));
    for l1=1:length(xreal2{1})
        finalres.pulse=union(ixvn1set,xreal2{1}(l1));
        [MUAPSet]=CalMUAPNonShift(size(signal,2),signal,finalres,[],size(signal,1));
        close all
        for llll=1:length(finalres)
            tallmuap=[];
            for jjj=1:size(signal,1)
                tallmuap=[tallmuap,MUAPSet{jjj,llll}];
                %plot(MUAPSet{jjj,llll});hold on
            end
            gap=max(tallmuap)-min(tallmuap);
            [~, indices] = maxk(abs(tallmuap), 100);
            % 将这些数从数组中移除
            tallmuap(indices) = [];
            tstd=std(tallmuap);
            MUAP_ratio(l1)=gap/tstd;
        end
    end
    [~, indices] = maxk(MUAP_ratio, 15);
    % 创建一个与MUAP_ratio同样大小的全0矩阵
    top15Matrixr = ones(size(MUAP_ratio));
    % 将最大的5个数的位置置为1
    top15Matrixr(indices) = 0;
    %%  联合计算连续多个IPT的PNR,SIL与幅值
%     clear Pnr2 SIL
%     for k=1:length(xreal{1})
%     for i=1:length(xreal2{1})
%         pulse=union(ixvn1set,xreal2{1}(i));
%         %Cstj=mean(save_all_ctj(pulse,:),1);
%         pulse1 = repmat(pulse, 1, 11) + kron(-5:5, ones(1, length(pulse)));
%         for j=-2:1:2
%             idx=true(1,5000);
%             array=pulse1+j;
%             array(array < 1 | array > 5000) = [];
%             array1=pulse+j;
%             array1(array1 < 1 | array1 > 5000) = [];
%             idx(array)=false;
%             tsnj=mean(save_all_ctj(array1,idx),1);
%             tjn1=mean(save_all_ctj(array1,array1),1);
%             Pnr2(i,j+3) = pnr(tjn1,tsnj);
%             %SIL(i,j+3)= sil(tjn1,tsnj);
%         end
%     end
%     end
%     rowsToRemove = any(SIL < 0.87, 2);
%     xreal{1}(rowsToRemove)=[];
%     Pnr(rowsToRemove,:)=[];
    amp=Cstj(xreal{1});
    amp2=Cstj(xreal2{1});
%     Pnr2=sum(Pnr,2);
    %% 根据IPT的PNR,幅值,FCE,已知MU的时刻等条件共同筛选
    % 如果某一时刻PNR第1，IPT幅值第1，那么一定为真实的放电时刻?
    % 如果没有，依次往下，若都前3的只有1个，那么只输出1个，不然就说明这个幅值不明显，取PNR前7的时刻
    % 根据FCE和已知MU信息，筛选掉剩余时刻
    dele_timel=FCE_RES|is_largeMU|top15Matrixl;
    xreal{1}=xreal{1}(~dele_timel);
%     Pnr=Pnr(~dele_time);
    amp=amp(~dele_timel);
%     [~,idx_PNR]=sort(Pnr);
    [~,idx_amp]=sort(amp);
    len1=length(amp);
    len1=min(5,len1-1);
    dele_timer=FCE_RES2|is_largeMU2|top15Matrixr;
    xreal2{1}=xreal2{1}(~dele_timer);
%     Pnr2=Pnr2(~dele_time);
    amp2=amp2(~dele_timer);
%     [~,idx_PNR2]=sort(Pnr2);
    [~,idx_amp2]=sort(amp2);
    len2=length(amp2);
    len2=min(5,len2-1);
    % 若识别不到时刻，就说明要么是FCE,要么是largeMU，说明已经没有了放电时刻
    if isempty(amp) || isempty(amp2)
        xreal=[];
        xreal2=[];
        Pnr_con=0;
        return;
    end
    % 联合计算连续多个IPT的PNR,SIL与幅值
    for k=1:length(xreal{1})
        for i=1:length(xreal2{1})
            pulse=union(ixvn1set,xreal2{1}(i));
            pulse=union(pulse,xreal{1}(k));
            %Cstj=mean(save_all_ctj(pulse,:),1);
            idx1=[];
            for j=1:length(pulse)-1
                if pulse(j+1)-pulse(j)> predict_gap/1.5 
                    idx1=[idx1,pulse(j)+20:pulse(j)+round(predict_gap/3),pulse(j+1)-round(predict_gap/3):pulse(j+1)-20];
                else
                    idx1=[idx1,pulse(j)+10:pulse(j+1)-10];
                end
            end
            pulse1=[];
            for j = -5:5
                pulse1 = [pulse1, j+pulse];
            end
            for j=-2:1:2
                idx=idx1+j;
                array=pulse1+j;
                array(array < 1 | array > window_len) = [];
                array1=pulse+j;               
                array1(array1 < 1 | array1 > window_len) = [];
                tsnj=mean(save_all_ctj(array1,idx),1);
                tjn1=mean(save_all_ctj(array1,array1),1);
                tPnr(j+3) = pnr(tjn1,tsnj);
            end
            Pnr_con(k,i)=sum(tPnr);
        end
    end
    [topValues, topIndices] = maxk(Pnr_con(:), 3);
    % 将索引转换为行列坐标
    [rowIndices, colIndices] = ind2sub(size(Pnr_con), topIndices);
    %先考虑len1=0（xreal{1}只有一个时刻）或者len2=0（xreal2{1}只有一个时刻）的case
    if len1==0 && len2>0
        if colIndices(1)==idx_amp2(end)
            xreal{1}=xreal{1}(rowIndices(1));
            xreal2{1}=xreal2{1}(colIndices(1));
            Pnr_con=Pnr_con(rowIndices(1),colIndices(1));
            return
        elseif ismember(colIndices(1),idx_amp2(end-1:end))
            xreal{1}=xreal{1}(rowIndices(1));
            xreal2{1}=xreal2{1}(colIndices(1));
            Pnr_con=Pnr_con(rowIndices(1),colIndices(1));
            return
        elseif ismember(colIndices(2),idx_amp2(end-1:end))
            xreal{1}=xreal{1}(rowIndices(1));
            xreal2{1}=xreal2{1}(colIndices(2));
            Pnr_con=Pnr_con(rowIndices(1),colIndices(2));
            return
        elseif ismember(colIndices(2),idx_amp2(end-len2:end))
            xreal{1}=xreal{1}(rowIndices(1));
            xreal2{1}=xreal2{1}(colIndices(2));
            Pnr_con=Pnr_con(rowIndices(1),colIndices(2));
            return
        else
            xreal2=[];
            Pnr_con=0;
            return
        end
    end
    if len2==0 && len1>0
        if rowIndices(1)==idx_amp(end)
            xreal{1}=xreal{1}(rowIndices(1));
            xreal2{1}=xreal2{1}(colIndices(1));
            Pnr_con=Pnr_con(rowIndices(1),colIndices(1));
            return
        elseif ismember(rowIndices(1),idx_amp(end-1:end))
            xreal{1}=xreal{1}(rowIndices(1));
            xreal2{1}=xreal2{1}(colIndices(1));
            Pnr_con=Pnr_con(rowIndices(1),colIndices(1));
            return
        elseif ismember(rowIndices(2),idx_amp(end-1:end))
            xreal{1}=xreal{1}(rowIndices(2));
            xreal2{1}=xreal2{1}(colIndices(1));
            Pnr_con=Pnr_con(rowIndices(2),colIndices(1));
            return
        elseif ismember(rowIndices(2),idx_amp(end-len1:end))
            xreal{1}=xreal{1}(rowIndices(2));
            xreal2{1}=xreal2{1}(colIndices(1));
            Pnr_con=Pnr_con(rowIndices(2),colIndices(1));
            return
        else
            xreal=[];
            Pnr_con=0;
            return
        end
    end
    if len2==0 && len1==0
        return;
    end
    if rowIndices(1)==idx_amp(end) && colIndices(1)==idx_amp2(end)
        xreal{1}=xreal{1}(rowIndices(1));
        xreal2{1}=xreal2{1}(colIndices(1));
        Pnr_con=Pnr_con(rowIndices(1),colIndices(1));
        return
    elseif ismember(rowIndices(1),idx_amp(end-1:end)) && ismember(colIndices(1),idx_amp2(end-1:end))
        xreal{1}=xreal{1}(rowIndices(1));
        xreal2{1}=xreal2{1}(colIndices(1));
        Pnr_con=Pnr_con(rowIndices(1),colIndices(1));
        return
    elseif ismember(rowIndices(2),idx_amp(end-1:end)) && ismember(colIndices(2),idx_amp2(end-1:end))
        xreal{1}=xreal{1}(rowIndices(2));
        xreal2{1}=xreal2{1}(colIndices(2));
        Pnr_con=Pnr_con(rowIndices(2),colIndices(2));
        return
    elseif ismember(rowIndices(2),idx_amp(end-len1:end)) && ismember(colIndices(2),idx_amp2(end-len2:end))
        xreal{1}=xreal{1}(rowIndices(2));
        xreal2{1}=xreal2{1}(colIndices(2));
        Pnr_con=Pnr_con(rowIndices(2),colIndices(2));
        return
    else
        xreal=[];
        xreal2=[];
        Pnr_con=0;
    end
 end

function [Pnr] = pnr(tjn1,tsnj)
    Pnr=10*log(mean(tjn1.*tjn1)/mean(tsnj.*tsnj));
end
function y= sil(tjn1,tsnj)
    m1=mean(tsnj);
    m2=mean(tjn1);
    t1=abs(tjn1-m1);
    t2=abs(tjn1-m2);
    s = (t1 - t2) ./ max(t1, t2);
    y=mean(s);
end
function MU_waveforms = extractMUWaveforms(sequence, discharge_times, waveform_length)
    % sequence: 1x5000长度的序列
    % discharge_times: 包含MU放电时刻的向量
    % waveform_length: 放电波形的长度
    
    num_times = length(discharge_times); % MU放电时刻的数量
    MU_waveforms = zeros(num_times, waveform_length); % 存储MU波形的数组
    
    for i = 1:num_times
        discharge_time = discharge_times(i); % 当前MU的放电时刻
        start_idx = discharge_time-22; % 波形起始索引
        end_idx = start_idx + waveform_length - 1; % 波形结束索引
        left_pad=max(1-start_idx,0);
        true_start_idx=max(start_idx,1);
        true_end_idx=min(end_idx,length(sequence));
        valid_length=true_end_idx-true_start_idx+1;
        MU_waveform=zeros(1,waveform_length);
        MU_waveform(left_pad+(1:valid_length))=sequence(true_start_idx:true_end_idx);
        % 提取当前MU的放电波形，使用向量切片进行提取
        MU_waveforms(i,:) = MU_waveform;
    end
end
function noise_indices=find_noise(signal,window_size,window_len)
    %换个方法找噪声时刻不就好了
    NEO_sequence=signal(2:window_len-1).*signal(2:window_len-1)-signal(1:window_len-2).*signal(3:window_len);
    % window_size:窗口大小
    % step_size: 滑动窗口的步长
    %noise_indices = [];
    seq_mean=mean(NEO_sequence);
    seq_std=std(NEO_sequence);
    noise_indices=find(NEO_sequence<0.2*(seq_mean+seq_std));
%     for i = 1:5:length(NEO_sequence) - window_size + 1
%         window_data = NEO_sequence(i:i+window_size-1);
%         window_mean=mean(window_data);
%         window_std = std(window_data); % 计算窗口内数据的标准差
%         if window_mean+window_std < 0.1*(seq_mean+seq_std)
%             noise_indices = [noise_indices, i:i+window_size-1];
%         end
%     end
    noise_indices=noise_indices+1;
    noise_indices=unique(noise_indices);
end
function result = step_function(x)
    if x >= 0
        result = 1;
    else
        result = 0;
    end
end
function hat_c_t_j_y = Eq10(nc, f_values,y)
hat_c_t_j_y=sum(f_values.*y(:,nc),2);
end