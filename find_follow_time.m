function [xreal,Pnr] = find_follow_time(precise_init_MU,save_all_ctj,predict_gap,signal,large_mu_time,zn_exten15wgn,window_len,it_point)
%rough_time ʹ��׼ȷ��ʼʱ�̵õ����Եķŵ�ʱ��
%   Ŀ�ĵõ� all good time and some bad time
%   ʹ�����Լ�ȨȨ����
    k_gap=0.3;
    ixvn1set=precise_init_MU;
    if nargin==8
        init_pt=it_point(end);
    else
        init_pt=ixvn1set(end);
    end
    c=init_pt+1;
    delta_a=min(c-1,window_len-2000);
    hat_t_j=mean(save_all_ctj(ixvn1set,:),1);
    n_c = c - delta_a : c;
    f_values = hat_t_j(n_c) .* (0.7 - (c - delta_a - n_c) / (2*delta_a));
    hat_c_t_j_y = Eq10(n_c, f_values,zn_exten15wgn);
    Cstj=hat_c_t_j_y'*zn_exten15wgn;
    %plot(Cstj)
    xreal=cell(1,1);
    %% ɸѡʱ�̷�Χ�ڷ�ֵ������������ʱ�̵�ֱ�����
    start_p=round(init_pt+(1-k_gap)*predict_gap);
    end_p=round(init_pt+(1+k_gap)*predict_gap);
    [c,ind]=max( Cstj(start_p:end_p));%�Ӿ�ֵ+3std����ֵ-3std�����ֵ
    ave_amp=mean( Cstj(start_p:end_p));%�õ�ƽ����ֵ
    std_amp=std( Cstj(start_p:end_p));%�õ�ƽ����׼��
    if c>ave_amp+2.5*std_amp
        xreal{1}=start_p+ind-1;
        Pnr=100;
        return;
    end
    %% IPT��ֵʱ��ɸѡ  
    [PKS,LOCS]= findpeaks(Cstj,'MinPeakDistance',3,'MinPeakHeight', 0);   
    for lp1=1:length(ixvn1set)
        LOCS(LOCS>ixvn1set(lp1)-round(0.5*predict_gap) & LOCS<ixvn1set(lp1)+round(predict_gap*0.5))=[];
    end
    LOCS=union(LOCS,ixvn1set);    
    xreal{1}=LOCS(LOCS>=init_pt+(1-k_gap)*predict_gap & LOCS<=init_pt+(1+k_gap)*predict_gap);
    %% NEO ����
    peak_detec=abs(Cstj(:,2:window_len-1)).*abs(Cstj(:,2:window_len-1))-abs(Cstj(:,1:window_len-2)).*abs(Cstj(:,3:window_len));
    max_values = max(peak_detec, [], 2);
    peak_detec=peak_detec./max_values;
    peak_detec=[0,peak_detec,0];
    [PKS,LOCS1]= findpeaks(peak_detec,'MinPeakDistance',3,'MinPeakHeight', 0); 
    for lp1=1:length(ixvn1set)
        LOCS1(LOCS1>ixvn1set(lp1)-round(0.5*predict_gap) & LOCS1<ixvn1set(lp1)+round(predict_gap*0.5))=[];
    end
    LOCS1=union(LOCS1,ixvn1set);
    xreal{2}=LOCS1(LOCS1>=init_pt+(1-k_gap)*predict_gap & LOCS1<=init_pt+(1+k_gap)*predict_gap);
    binary_xreal=Cstj(xreal{2})>0;
    xreal{2}=xreal{2}(binary_xreal);
    xreal{1}=union(xreal{1},xreal{2});
    % ���xreal{1}ֱ�Ӳ����ڣ���ôreturn
    if  isempty(xreal{1})
        Pnr=0;
        return;
    end

    %% FCEɸѡ�����Դ����ʱ��,Ҫע�⣬FCE���Ҫ��߽�㲻��С��25������4975,
    % ����һ��Ԥ��500����ʱ�̣���˲��ᵽ��4975�߽��
    % error_pulse�Ǽ�����ʱ�̵�pulse
    FCE_RES=zeros(1,length(xreal{1}));
    COR_RES=zeros(1,length(xreal{1}));
    channel_num=round(size(signal,1)/2);
    for l1=1:length(xreal{1})
        error_pulse= xreal{1}(l1);
        temp_pulse1=ixvn1set(ixvn1set>25);%�߽�㲻��С��25
        temp_pulse1=temp_pulse1(temp_pulse1<window_len-25);
        temp_pulse=union(temp_pulse1,error_pulse);
        temp_pulse(temp_pulse<1 & temp_pulse>window_len)=[];
        idx = find(temp_pulse == error_pulse);
        temp_pulse1(temp_pulse1<1 & temp_pulse1>window_len)=[];
        waveform_length=45;
        MU_waveforms = extractMUWaveforms(signal(channel_num,:), temp_pulse, waveform_length);
        MU_waveforms1 = extractMUWaveforms(signal(channel_num,:), temp_pulse1, waveform_length);
        % % e:�����ľ�����ֵ ˼·������neo����ȥȷ������ʱ�̣���ʹ������ʱ�̵õ��������ֵ RMSE = sqrt(sum(noise_data.^2) / n)
        window_size = 15; 
        noise_time=find_noise(signal(channel_num,:),window_size,window_len);
        noise_data=signal(channel_num,noise_time);
        e = sqrt(sum(noise_data.^2) / length(noise_data));
        % �����Ѿ��������±�����
        % m: ѡ���������
        % ��ʼ��g[n]����
        % ����g[n]����ÿ�����������һ�������Ĳ�ֵ
        % yi �ݶ�Ϊ��ֵ
        yi=mean(MU_waveforms1,1);
        xi=MU_waveforms(idx,:);
        g = abs(diff(xi));
        m=4;
        % ѡ��������g[n]ֵ��m�������������
        [~, sorted_idx] = sort(g, 'descend');
        selected_idx = sorted_idx(1:m);
        alpha=0.01;
        beta=0.05;
        % ��ȡѡ������������Ϊ����
        features = xi(selected_idx);
        % ����MUPTģ���MUPģ���ź�Ϊyi(selected_idx)����СΪm��������
        % ������״��һ���԰ٷֱ�PSI_i
        diff1 = yi(selected_idx) - features ;
        term1=0;term2=0;
        for ii=1:length(diff1)
            term1 = term1+step_function(diff1(ii)- 3 * e);
            term2 = term2+step_function(-diff1(ii)- 3 * e);
        end
        PSI_i = (sum(term1) + sum(term2)) / m;
        % ����d_i
        d_i = sum(diff1.^2) / e^2;
        % ��ʾ���
        %fprintf('MUP %d ����״��һ���԰ٷֱ� (PSI): %f\n', PSI_i);
        chi2_alpha = chi2inv(1 - alpha, m); % ������� m=3
        if d_i > chi2_alpha && PSI_i > 0.7
            %fprintf('xi Definitely a FCE\n');
            FCE_RES(l1)=1;
        elseif d_i <chi2inv(1 - beta, m)
            COR_RES(l1)=1;
            %fprintf('Definitely a correct\n');
        end
        
    end
    %% ɸѡ��֪MU��ʱ�̣������ų�������е��Ѹ㣬�������3���ˣ�˵���ص�ʱ��̫���ˣ�ȥ��һ��
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
                pulse1 = bsxfun(@plus, pulse, -2:2);
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
    %%  ɸѡ��ʹ�� MUAP������õ�15��ʱ��
    MUAP_ratio=zeros(1,length(xreal{1}));
    for l1=1:length(xreal{1})
        finalres.pulse=union(ixvn1set,xreal{1}(l1));
        [MUAPSet]=CalMUAPNonShift1(size(signal,2),signal,finalres,[],size(signal,1));
        for llll=1:length(finalres)
            tallmuap=[];
            for jjj=1:size(signal,1)
                tallmuap=[tallmuap,MUAPSet(jjj,:)];
                %plot(MUAPSet{jjj,llll});hold on
            end
            gap=max(tallmuap)-min(tallmuap);
            [~, indices] = maxk(abs(tallmuap), 100);
            % ����Щ�����������Ƴ�
            tallmuap(indices) = [];
            tstd=std(tallmuap);
            MUAP_ratio(l1)=gap/tstd;
        end
    end
    [~, indices] = maxk(MUAP_ratio, 15);
    % ����һ����MUAP_ratioͬ����С��ȫ0����
    top15Matrix = ones(size(MUAP_ratio));
    % ������5������λ����Ϊ1
    top15Matrix(indices) = 0;
    %%  �������IPT��PNR,SIL���ֵ
    clear Pnr SIL
    for i=1:length(xreal{1})
        pulse=union(ixvn1set,xreal{1}(i));
        %Cstj=mean(save_all_ctj(pulse,:),1);
        %pulse1 = bsxfun(@plus, pulse, -5:5);
        pulse1=[];
        for j = -5:5
            pulse1 = [pulse1, j+pulse];
        end
        for j=-2:1:2
            idx=true(1,window_len);
            array=pulse1+j;
            array(array < 1 | array > window_len) = [];
            array1=pulse+j;
            array1(array1 < 1 | array1 > window_len) = [];
            idx(array)=false;
            tsnj=mean(save_all_ctj(array1,idx),1);
            tjn1=mean(save_all_ctj(array1,array1),1);
            Pnr(i,j+3) = pnr(tjn1,tsnj);
            %SIL(i,j+3)= sil(tjn1,tsnj);
        end
    end
%     rowsToRemove = any(SIL < 0.87, 2);
%     xreal{1}(rowsToRemove)=[];
%     Pnr(rowsToRemove,:)=[];
    amp=Cstj(xreal{1});
    Pnr=sum(Pnr,2);
    %% ����IPT��PNR,��ֵ,FCE,��֪MU��ʱ�̵�������ͬɸѡ
    % ���ĳһʱ��PNR��1��IPT��ֵ��1����ôһ��Ϊ��ʵ�ķŵ�ʱ��?
    % ���û�У��������£�����ǰ3��ֻ��1������ôֻ���1������Ȼ��˵�������ֵ�����ԣ�ȡPNRǰ7��ʱ��
    % ����FCE����֪MU��Ϣ��ɸѡ��ʣ��ʱ��
    dele_time=FCE_RES|is_largeMU|top15Matrix;
    xreal{1}=xreal{1}(~dele_time);
    len=length(xreal{1});
    len=min(5,len-1);
    Pnr=Pnr(~dele_time);
    amp=amp(~dele_time);
    [~,idx_PNR]=sort(Pnr);
    [~,idx_amp]=sort(amp);
    % ��ʶ�𲻵�ʱ�̣���˵��Ҫô��FCE,Ҫô��largeMU��˵���Ѿ�û���˷ŵ�ʱ��
    if isempty(Pnr)
        return;
    end
    % ���������ϸ�����
    if idx_PNR(end)==idx_amp(end)
        xreal{1}=xreal{1}(idx_PNR(end));
        Pnr=Pnr(idx_PNR(end));
        return
    elseif any(ismember(idx_PNR(end-1:end),idx_amp(end-1:end)))
        c=intersect(idx_PNR(end-1:end),idx_amp(end-1:end));
        xreal{1}=xreal{1}(c);
        Pnr=Pnr(c);
        return
    elseif any(ismember(idx_PNR(end-2:end),idx_amp(end-2:end)))
        c=intersect(idx_PNR(end-2:end),idx_amp(end-2:end));%%%%
        xreal{1}=xreal{1}(c);
        Pnr=Pnr(c);
        return
    elseif any(ismember(idx_PNR(end-len:end),idx_amp(end-len:end)))
        c=intersect(idx_PNR(end-len:end),idx_amp(end-len:end));%%%%
        xreal{1}=xreal{1}(c);
        Pnr=Pnr(c);
        return
    else
%         xxx=max(floor(length(Pnr)/2),length(Pnr)-5);
%         remove_idx=idx_PNR(1:xxx);
%         xreal{1}(remove_idx)=[];
%         Pnr(remove_idx)=[];
        
        xreal{1}=[];
        Pnr=[];
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
    % sequence: 1x5000���ȵ�����
    % discharge_times: ����MU�ŵ�ʱ�̵�����
    % waveform_length: �ŵ粨�εĳ���
    
    num_times = length(discharge_times); % MU�ŵ�ʱ�̵�����
    MU_waveforms = zeros(num_times, waveform_length); % �洢MU���ε�����
    
    for i = 1:num_times
        discharge_time = discharge_times(i); % ��ǰMU�ķŵ�ʱ��
        start_idx = discharge_time-22; % ������ʼ����
        end_idx = start_idx + waveform_length - 1; % ���ν�������
        start_idx=max(start_idx,1);
        left_pad=max(1-start_idx,0);
        true_start_idx=max(start_idx,1);
        true_end_idx=min(end_idx,length(sequence));
        valid_length=true_end_idx-true_start_idx+1;
        MU_waveform=zeros(1,waveform_length);
        MU_waveform(left_pad+(1:valid_length))=sequence(true_start_idx:true_end_idx);
        % ��ȡ��ǰMU�ķŵ粨�Σ�ʹ��������Ƭ������ȡ
        MU_waveforms(i,:) = MU_waveform;
    end
end
function noise_indices=find_noise(signal,window_size,window_len)
    NEO_sequence=signal(2:window_len-1).*signal(2:window_len-1)-signal(1:window_len-2).*signal(3:window_len);
    % window_size:���ڴ�С
    % step_size: �������ڵĲ���
    %noise_indices = [];
    seq_mean=mean(NEO_sequence);
    seq_std=std(NEO_sequence);
    noise_indices=find(NEO_sequence<0.2*(seq_mean+seq_std));
%     for i = 1:5:length(NEO_sequence) - window_size + 1
%         window_data = NEO_sequence(i:i+window_size-1);
%         window_mean=mean(window_data);
%         window_std = std(window_data); % ���㴰�������ݵı�׼��
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
function [MUAPAver]=CalMUAPNonShift1(NumSample,EMGFiltSeg,PulseFinal,Corrupted_Channel,CHANNEL_NUM)
    if (nargin<5)
        CHANNEL_NUM=128;
    end                
    if isnumeric(PulseFinal)
        % ����һ���µĽṹ�岢������洢�����е�һ���ֶ���
        PulseFinal = struct('pulse', PulseFinal);
    end
HalfWidth=100;
%MUAPSet=zeros(CHANNEL_NUM,2*HalfWidth+1);%Even there are corrupted channel that we have been deleted, we still need to draw the 128-channel MUAP
for loop2=1:length(PulseFinal)
   %Ϊ�˼�������ʱ�䣬����ֻ��ȡ���3000��ʱ�̵�MUAP��
   temp_pulse=PulseFinal(loop2).pulse;
   MUAP=zeros(CHANNEL_NUM,2*HalfWidth+1);%201 x ���г���

   for loop3=1:length(temp_pulse)%�õ���MU���еĳ���
    
        RightIPL=temp_pulse(loop3)-HalfWidth;%��loop3ʱ��ֵ��100
        RightIPH=temp_pulse(loop3)+HalfWidth;%��loop3ʱ��ֵ��100
        %����߽��ڷ�Χ�ڣ����¼loop3ʱ��loop1ͨ��������100��ʱ�̣��ܹ�201��ʱ�̣��ķ���ֵ
        if RightIPL>1 && RightIPH<NumSample                
           MUAP=MUAP+EMGFiltSeg(:,RightIPL:RightIPH);%64*201
        end                        
   end
   MUAPAver=MUAP/length(temp_pulse);     
end      
 

end