function finalres = Dynamic_for_real_stric_MUAP1_visual(sEMG, start_pt, end_pt, pulsefinal1, gap_save)
if nargin < 6 || isempty(gap_save)
    gap_save = [150,170,80,150,80,80,80,150,80,80,90,85,90,90,90,90,90,90,90,90,90];
end
init_time=start_pt+99;%start_point(sub,Task)+99;
end_time=end_pt;%init_time+sample_window;%end_point(sub,Task); 
gap_time=end_time-init_time;
window_len=10000; %这是每个可能有问题间隔，用于分解的时间窗长度，比如时刻2出现问题，就分解2-10002的时间窗
one_firing_gap_bound=200;
if end_time-init_time<5000
    end_time=init_time+2000;
end

% 提取信号及预处理
signal = sEMG(:, start_pt:end_pt);
zn_exten15wgn_all=get_zn(signal);
save_all_ctj_all=zn_exten15wgn_all'*zn_exten15wgn_all;
[y,~]=get_zn2(signal);
y=y(:,101:end);

% 更新每个MU的初始间隔gap_save
for MU=1:1:length(pulsefinal1)
    data=diff(pulsefinal1(MU).pulse);
    lowerBound = prctile(data, 20);
    upperBound = prctile( data, 80);
    % 筛选处于中间区间的数据
    filteredData = data(data >= lowerBound & data <= upperBound);
    % 计算中间数据的均值
    gap_save(MU) = mean(filteredData);
end
for MU=1:1:length(pulsefinal1)
    finalres(MU).gap=round(gap_save(MU));
end

for MU=1:1:length(pulsefinal1)
% 根据初始时刻确定时间窗范围
finalres(MU).pulse=pulsefinal1(MU).pulse+init_time;
last_time_num=length(finalres(MU).pulse);
while 1
%% 首先遍历整个时间窗内的good time，再根据这些good time 进一步筛选剩余时刻:global scan
%true_pulse=precise_init_MU(MU).res(iidx).pulse;
true_pulse=finalres(MU).pulse-init_time;
%signal=semg2sign(sEMG,init_time-99,end_time); 
a1 = true_pulse(1);           % a1 的初始值
a2 = true_pulse(end);         % a2 的初始值
delta_a = 2000;     % 初始的 delta_a 值,整个收缩的1/5，主要是为了把靠近后面时刻的权向量的权重加大
if a2+1+2000>=gap_time
    delta_a=100*ceil((gap_time-a2)/200);
end
spike_set=[];
t_j=zeros(1,size(y,2));% 初始化 t_j 数组
found_c_1=true;
found_c_2=true;
c_t_j_y=zeros(size(y,1),1);
Upsilon=[];% 创建一个空集合
na=a1:a2;
firing_t=true_pulse(true_pulse>=a1 & true_pulse<=a2);
Cxn_exten=y(:,na) * y(:,na)';
if isnan(det(Cxn_exten))
        Cxn_exten_inv=pinv(Cxn_exten);
elseif ~det(Cxn_exten)
    Cxn_exten_inv=pinv(Cxn_exten);
else
    Cxn_exten_inv=inv(Cxn_exten);
end
hat_t_j= mean(y(:,firing_t),2)'*Cxn_exten_inv*y(:,na);
    
% 使用 Eq. (8) 计算 c_t_j_y_n_a 在 [a1, a2] 上的值
% 注意，y是扩展之后的信号
c_t_j_y_n_a =sum(hat_t_j.*y(:,na),2) ;
% 循环直到 b1 + delta_a 不再超过信号长度
b1 = a2+1;
if finalres(MU).gap==0
    sequence=finalres(MU).pulse;
    subsets = {};
    currentSubset = sequence(1);
    % 循环遍历序列
    diff_se=diff(sequence);
    %med = median(diff_se);
    %diff_se(2*min(diff_se)<diff_se)=[];%去掉过大的间隔时刻
    med=1.3*min(diff_se);% 保证放电时刻的不太大也不太小
    for i = 2:length(sequence)
        if sequence(i) - sequence(i-1) <= 1.4*med % 保证不能让间隔太大的时刻进来，但是稍微大一点的要进来
            % 如果差值小于等于300，则将当前元素加入当前子集
            currentSubset = [currentSubset, sequence(i)];
        else
            % 如果差值大于300，则结束当前子集并开始一个新的子集
            subsets{end+1} = currentSubset;
            currentSubset = sequence(i);
        end
    end
    % 结束循环后，将最后一个子集添加到结果中
    subsets{end+1} = currentSubset;
    if length(subsets)==1
        break
    end
    if length(subsets)==length(sequence)%如果相等说明没法计算间隔，只能用初始值
        firing_gap=100;%初始值设置为100
    else
        t_diff=[];% 假设初始的放电时刻是100
        for i=1:length(subsets)
            if length(subsets{i})>=2
                t_diff=[t_diff,diff(subsets{i})];
            end
        end
        firing_gap=round(mean(t_diff));
    end
else
    firing_gap=finalres(MU).gap;
end
%firing_gap
%firing_gap=floor((a2-a1)/length(true_pulse));%%%%%%%%这里可能存在问题，因为这个gap不太准确
while b1 + delta_a <= size(y,2)
    % timeset存储这次循环的时间集合,将旧时刻与新时刻结合，检测新时刻中放电
    nb=b1:(b1+delta_a);
    timeset=[na, nb];
    % 这里假设计算结果存储在变量 hat_t_j_n_a 中(对应timeset)
    hat_t_j_n_a_n = CalculateHatT_j_n_a_n(c_t_j_y_n_a, y,timeset);           
    % 分割 hat_t_j_n_a 成尖峰和基线噪声，使用 Eq. (7) 计算 PNR
    [spikes,~] = Split_spike(hat_t_j_n_a_n,y,timeset,firing_gap,delta_a);%%%%%%%%%%%%%%%%%%这个也可能需要修改
    %plot(tjns1)
    % idx 是噪声时刻
    hat_t_j_n_a_nb=hat_t_j_n_a_n(end-delta_a:end);
    spike=spikes-length(na);
    spike(spike<1)=[];
    idx=true(1,length(hat_t_j_n_a_nb));
    idx(spike)=false;
    tjn1=hat_t_j_n_a_nb(spike);
    tsnj=hat_t_j_n_a_nb(idx);
    Pnr = pnr(tjn1,tsnj);
    % 如果在初始序列中，这些时刻的PNR都大于28，说明这些时刻非常靠谱，确实是这样
    % 判断 PNR 是否大于等于 28 dB,如果太苛刻，适当降低也无不可
    if Pnr >= 26
        spike_set=union(spike_set,nb(spike));
        % 更新 Upsilon
        Upsilon = union(Upsilon, nb);
    end
    % 更新 b1 的值
    b1 = b1 + delta_a;
end
if ~isempty(spike_set)
    finalres(MU).pulse=union(finalres(MU).pulse,spike_set+init_time);
    ixvn1set=finalres(MU).pulse-init_time;
    tjns1= mean(y(:,ixvn1set),2)'*Cxn_exten_inv*y;
    ixvn1set=dele_adjacent_time(ixvn1set,tjns1,firing_gap/2);
    finalres(MU).pulse=ixvn1set+init_time;
end
%finalres1(MU).pulse=finalres(MU).pulse-
% 观察这一步会不会出现问题,没出现问题
% [MUAPSet]=CalMUAPNonShift(9000,sEMG(:,1:9000),finalres,[],size(sEMG,1));
% close all
% for jjj=1:size(sEMG,1)
%     plot(MUAPSet{jjj,MU});hold on
% end

%% add missing time  对于间隔小于2000样本点，尝试补全时刻
% 首先判断间隔小于2000样本点
% 如果小于2000，每一次循环补全中间时刻
% 初始化子集
sequence=finalres(MU).pulse;
subsets = {};
currentSubset = sequence(1);
% 循环遍历序列
diff_se=diff(sequence);
med=1.3*min(diff_se);% 保证放电时刻的不太大也不太小
for i = 2:length(sequence)
    if sequence(i) - sequence(i-1) <= 1.4*med
        % 如果差值小于等于300，则将当前元素加入当前子集
        currentSubset = [currentSubset, sequence(i)];
    else
        % 如果差值大于300，则结束当前子集并开始一个新的子集
        subsets{end+1} = currentSubset;
        currentSubset = sequence(i);
    end
end
% 结束循环后，将最后一个子集添加到结果中
subsets{end+1} = currentSubset;
if length(subsets)==1
    break
end
t_diff=[];
for i=1:length(subsets)
    if length(subsets{i})>=2
        t_diff=[t_diff,diff(subsets{i})];
    end
end
predict_gap=round(mean(t_diff));
for i=2:length(subsets)
    while subsets{i}(1)-subsets{i-1}(end)>2.5*predict_gap && subsets{i}(1)-subsets{i-1}(end)<2000 
        %得到左边时刻与右边时刻
        %根据缺失段定义起始时刻与终止时刻
        it=subsets{i-1}(end)-2000;
        it=max(it,100); 
        it=max(start_pt+100,it);%防止加入还没放电时候的信息
        if it+window_len<end_time
            signal=semg2sign(sEMG,it-99,it+window_len);
        else
            it=end_time-window_len;
            signal=semg2sign(sEMG,it-99,it+window_len);
        end
        zn_exten15wgn=zn_exten15wgn_all(:,it-init_time:it+window_len-init_time);
        save_all_ctj=save_all_ctj_all(it-init_time:it+window_len-init_time,it-init_time:it+window_len-init_time);
        %zn_exten15wgn=get_zn(signal);
        %save_all_ctj=zn_exten15wgn'*zn_exten15wgn;
        % 根据MU确定已知时刻
        temp_pulse=finalres(MU).pulse-it;  
        temp_pulse(temp_pulse<1)=[];
        temp_pulse(temp_pulse>window_len)=[];
        left_time=subsets{i-1}(end)-it;
        right_time=subsets{i}(1)-it;
        large_mu_time=[];
        %补全中间的时刻
        [xreal,xreal2,Pnr_con] = find_middle_2time(temp_pulse,save_all_ctj,predict_gap,...
            signal(:,101:end),large_mu_time,zn_exten15wgn,left_time,right_time,window_len);
        %判断是否有新时刻，如果没有，直接结束
        if isempty(xreal) && isempty(xreal2)
            break;
        elseif ~isempty(xreal) && isempty(xreal2)
            subsets{i-1}=union(subsets{i-1},xreal{1}+it);
        elseif isempty(xreal) && ~isempty(xreal2)
            subsets{i}=union(subsets{i},xreal2{1}+it);
        else
            subsets{i-1}=union(subsets{i-1},xreal{1}+it);
            subsets{i}=union(subsets{i},xreal2{1}+it);
        end
    end
end
for i=1:length(subsets)
%补全后更新时刻
finalres(MU).pulse=union(finalres(MU).pulse,subsets{i});
end

if last_time_num == length(finalres(MU).pulse)
    break
end
last_time_num = length(finalres(MU).pulse);
end
%MU
end
%%
signal=semg2sign(sEMG,init_time-99,end_time); 
[y,~]=get_zn2(signal);
zn_exten15wgn=get_zn(signal);
save_all_ctj=zn_exten15wgn'*zn_exten15wgn;
for i=1:length(finalres)
    [finalres(i).mean_firing_gap,finalres(i).cov]=predictgap(finalres(i).pulse);
    finalres(i).mean_firing_rate=2000/finalres(i).mean_firing_gap;
    finalres(i).firing_time=finalres(i).pulse/2048;
    finalres(i).Cst=mean(y(:,finalres(i).pulse-(init_time-99-1)),2)';
    ctsj=mean(save_all_ctj(finalres(i).pulse-init_time,:),1);
    finalres(i).ctsj=ctsj;
end
% 存储最新的时刻
%从predict_time中人工确定一组时刻作为结果
%save_address=strcat('D:\dynamic_decompose\real data\',str{sub},'\',str{sub},'_Dynamic_Task',num2str(Task),'_Trial',num2str(trial),'_res_sys4.mat');
%save(save_address,'finalres');
%% 观察新生成的MUAP是否正常
% [MUAPSet]=CalMUAPNonShift(end_time,sEMG(:,1:end_time),finalres,[],size(sEMG,1));
% MUAP_ratio=zeros(1,length(finalres));
% for llll=length(finalres):-1:1
%     close all
%     tallmuap=[];
%     llll
%     figure_handle = figure('Position', [100, 100, 800, 1200]);
%     subplot(2, 1, 1);
%     for jjj=1:size(sEMG,1)
%         tallmuap=union(tallmuap,MUAPSet{jjj,llll});
%         plot(MUAPSet{jjj,llll});hold on
%     end
%     gap=max(tallmuap)-min(tallmuap);
%     [~, indices] = maxk(abs(tallmuap), round(length(tallmuap)/20));
%     subplot(2, 1, 2); % 创建一个新的图形窗口
%     plot(finalres(llll).ctsj, 'b'); % 'b' 表示蓝色线条
%     % 将这些数从数组中移除
%     tallmuap(indices) = [];
%     tstd=std(tallmuap);
%     gap/tstd
%     MUAP_ratio(llll)=gap/tstd;
% end
% MUAP_binary = MUAP_ratio > 15;%%%%%%%%%%%%%%%%%%%%%%%%%%
% finalres=finalres(MUAP_binary);
%save(save_address,'finalres');
%precision_pro_realdata
%pulsefinal1(bad_num)=[];

end


function zn_exten15wgn = Whiten(xn_exten15wgn)
    % Whiten - 对输入矩阵进行白化处理
    % 输入:
    %   xn_exten15wgn - 原始扩展矩阵，每行代表一个信号
    % 输出:
    %   zn_exten15wgn - 白化后的矩阵

    % 去均值
    for i = 1:size(xn_exten15wgn, 1)
        xn_exten15wgn(i, :) = xn_exten15wgn(i, :) - mean(xn_exten15wgn(i, :));
    end
    
    % 计算协方差矩阵
    Cxn_exten = xn_exten15wgn * xn_exten15wgn' / size(xn_exten15wgn, 2);
    
    % 通过 SVD 计算白化矩阵
    [u, s, v] = svd(Cxn_exten);
    W = v * diag(1./sqrt(diag(s))) * v';
    
    % 得到白化后的矩阵
    zn_exten15wgn = W * xn_exten15wgn;
end

function [xn_exten15wgn,Cxn_exten_inv]=get_zn2(signal)
    M=size(signal,1);
    K=15;
    nsamples=size(signal,2);
    xn_exten15wgn=zeros(M*K,nsamples);%M*K,n 
    for i=1:M
        for j=1:K
            xn_exten15wgn(K*i+j-K,:)=[zeros(1,j-1),signal(i,1:end-j+1)];
        end
    end
    for i =1:size(xn_exten15wgn,1)
        xn_exten15wgn(i,:)=xn_exten15wgn(i,:)-mean(xn_exten15wgn(i,:));
    end
    Cxn_exten=xn_exten15wgn*xn_exten15wgn'/nsamples;
    if isnan(det(Cxn_exten))
        Cxn_exten_inv=pinv(Cxn_exten);
    elseif ~det(Cxn_exten)
        Cxn_exten_inv=pinv(Cxn_exten);
    else
        Cxn_exten_inv=inv(Cxn_exten);
    end
    %zn_exten15wgn=Whiten(xn_exten15wgn(:,101:end));
end
function [zn_exten15wgn,xn_exten15wgn] =get_zn(signal)
    M=size(signal,1);
    K=15;
    nsamples=size(signal,2);
    xn_exten15wgn=zeros(M*K,nsamples);%M*K,n 
    for i=1:M
        for j=1:K
            xn_exten15wgn(K*i+j-K,:)=[zeros(1,j-1),signal(i,1:end-j+1)];
        end
    end
    for i =1:size(xn_exten15wgn,1)
        xn_exten15wgn(i,:)=xn_exten15wgn(i,:)-mean(xn_exten15wgn(i,:));
    end
    zn_exten15wgn=Whiten(xn_exten15wgn(:,101:end));
    xn_exten15wgn=xn_exten15wgn(:,101:end);
end

function signal=semg2sign(sEMG,start_point,end_point)
    if mod(size(sEMG,1),2)==0
        l=size(sEMG,1)/2;
    else
        l=(size(sEMG,1)-1)/2;
    end
    if size(sEMG,1)>100
        for i=1:l
            signal(i,:)=sEMG(2*i,start_point:end_point);
        end
    else
        signal=sEMG(:,start_point:end_point);
    end
end

function hat_t_j_n_a_n = CalculateHatT_j_n_a_n(c_t_j_y_n_a, y,timeset)
    % c_t_j_y_n_a: 已知的 c_t_j_y_n_a 向量
    % y: 信号向量
    % n: 要计算的时间点 n
    % a1, a2: 时间范围 [a1, a2]
    
    y_matrix= y(:,timeset) * y(:,timeset)';
    if isnan(det(y_matrix))
            Cxn_exten_inv=pinv(y_matrix);
    elseif ~det(y_matrix)
        Cxn_exten_inv=pinv(y_matrix);
    else
        Cxn_exten_inv=pinv(y_matrix);
    end
    hat_t_j_n_a_n = c_t_j_y_n_a' * Cxn_exten_inv*y(:,timeset); 
end

function [ixvn1set,tjns1] = Split_spike(tjns1,y,timeset,firing_gap,delta_a)
    % 这个函数要怎么写，才能达到一个比较O,K的结果呢
    % 直接用CKC的好了
    Nnloop=round(delta_a/firing_gap)-5;
    Nnloop=min(Nnloop,10);
    ixvn1set=[];
    Nfind=round((length(timeset)-delta_a)/firing_gap)+Nnloop-1;
    vn_1=tjns1;
    ixvn1set=[];
    %这里每次都从头开始迭代实在太慢了，我得加快迭代速度
    while 1
        [~,ixvn1]=max(vn_1);
        if ~any(ixvn1set==ixvn1)%若不在集合里面 或者 一开始ixvn1set为空
            nvr=ixvn1;
            vn_1(ixvn1)=-inf;
        else
            vn_1(ixvn1)=-inf;
        end
        ixvn1set=[ixvn1set,nvr];%ixvn1set就是这一次最大点的时刻集合
        ixvn1set=sort(ixvn1set);
        ixvn1set=dele_adjacent_time(ixvn1set,tjns1,firing_gap/1.5);
        if length(ixvn1set)>=Nfind 
                break; 
        end
    end
    %tjns1 = mean(y(:,timeset(ixvn1set)),2)' * (tmp); 
    %plot(tjns1)
end

function [Pnr] = pnr(tjn1,tsnj)
    Pnr=10*log(mean(tjn1.*tjn1)/mean(tsnj.*tsnj));
end

function [predict_gap1,predict_cov] =predictgap(pulse)
    sequence=pulse;
    subsets = {};
    currentSubset = sequence(1);
    % 循环遍历序列
    diff_se=diff(sequence);
    med=1.3*min(diff_se);% 保证放电时刻的不太大也不太小
    for i = 2:length(sequence)
        if sequence(i) - sequence(i-1) <= 1.5*med
            % 如果差值小于等于300，则将当前元素加入当前子集
            currentSubset = [currentSubset, sequence(i)];
        else
            % 如果差值大于300，则结束当前子集并开始一个新的子集
            subsets{end+1} = currentSubset;
            currentSubset = sequence(i);
        end
    end
    % 结束循环后，将最后一个子集添加到结果中
    subsets{end+1} = currentSubset;
    if length(subsets)==1
        predict_gap1=med;
        predict_cov=std(sequence)/predict_gap1;
        return
    end
    t_diff=[];
    for i=1:length(subsets)
        if length(subsets{i})>=2
            t_diff=[t_diff,diff(subsets{i})];
        end
    end
    predict_gap1=round(mean(t_diff));
    predict_cov=std(t_diff)/predict_gap1;
end