
function finalres = Dynamic_for_real_stric_MUAP2_visual(sEMG, start_pt, end_pt, finalres)
init_time=start_pt+99;
end_time=end_pt;
window_len=10000;
large_mu_time=[];
large_mu_num=0; 
large_mu_idx=[];
small_mu_idx=[];
one_firing_gap_bound=200;
temp_large_mu_time=[];
signal=semg2sign(sEMG,init_time-99,end_time); 
zn_exten15wgn_all=get_zn(signal);
save_all_ctj_all=zn_exten15wgn_all'*zn_exten15wgn_all;
mid_time=(init_time-99+end_time)/2;
for MU=1:1:length(finalres)
    predict_gap=precise_gap(finalres(MU).pulse);
    if length(finalres(MU).pulse)>=1.8*(mid_time-finalres(MU).pulse(1))/predict_gap
        large_mu_num=large_mu_num+1;
        large_mu_time{large_mu_num}=finalres(MU).pulse;
        large_mu_idx=[large_mu_idx,MU];
    else
        small_mu_idx=[small_mu_idx,MU];
    end
end
% 对于那些没识别到的MU，利用大MU判别时刻
% 如果可以补全，尽量补全
%%
for MU=small_mu_idx
    if  ~any(diff(finalres(MU).pulse)>2.5*finalres(MU).mean_firing_gap) 
        %如果目前的结果识别到比较少的连续时刻，直接往后迭代就行了
        while finalres(MU).pulse(end)< end_time-2000
            it=finalres(MU).pulse(end)-(window_len-2000);
            is_finish=false;
            it=max(it,start_pt+100);
            clear predict_time predict_pnr
            signal=semg2sign(sEMG,it-99,it+window_len); 
            zn_exten15wgn=zn_exten15wgn_all(:,it-init_time:it+window_len-init_time);
            save_all_ctj=save_all_ctj_all(it-init_time:it+window_len-init_time,it-init_time:it+window_len-init_time);
            % 根据MU确定已知时刻
            temp_pulse=finalres(MU).pulse-it;  
            temp_pulse(temp_pulse<1)=[];
            if ~isempty(large_mu_time)
                for l2=1:length(large_mu_time)
                    pulse=large_mu_time{l2}-it;
                    pulse(pulse<10 | pulse>window_len-100)=[];
                    temp_large_mu_time{l2}=pulse;
                end
                temp_large_mu_time = temp_large_mu_time(~cellfun('isempty', temp_large_mu_time));
            end
            [xreal,Pnr] = find_follow_time(temp_pulse,save_all_ctj,...
                finalres(MU).mean_firing_gap,signal(:,101:end),temp_large_mu_time,zn_exten15wgn,window_len);
            if isempty(xreal{1})
                is_finish=1;
                break;
            end
            predict_time1=xreal{1};
            predict_time=num2cell(predict_time1);
            ttt=union(temp_pulse,predict_time{1});
            UU_iter=0;
            while predict_time{1}(end)<window_len-600 && is_finish==0
                UU_iter=UU_iter+1;
                while length(predict_time)<50 && predict_time{1}(end)<window_len-600
                    clear t_p t_pnr
                    for i=1:length(predict_time)
                        if predict_time{i}(end)<window_len-600
                            temp_pulse1=union(temp_pulse,predict_time{i});
                            [xreal,Pnr] = find_follow_time(temp_pulse1,save_all_ctj,...
                                finalres(MU).mean_firing_gap,signal(:,101:end),temp_large_mu_time,zn_exten15wgn,window_len);
                %             toRemove = cellfun(@(x) length(x) == 1, xreal);
                %             xreal(toRemove) = [];
                            for j=1:length(xreal{1})
                                t_p{i}{j}=union(predict_time{i},xreal{1}(j));
                                t_pnr{i}{j}=Pnr(j);
                            end
                        end
                    end
                    % 如果t_p不存在，则说明没有新时刻了
                    if ~exist('t_p','var')
                        is_finish=1;
                        break;
                    end
                    % 用predict_time存储铺平后的t_p：这一轮的所有时刻组合
                    flattenedArray = cellfun(@(x) x(:), t_p, 'UniformOutput', false);
                    predict_time = cat(1, flattenedArray{:});
                    flattenedArray = cellfun(@(x) x(:), t_pnr, 'UniformOutput', false);
                    predict_pnr = cat(1, flattenedArray{:});
                    %find(cellfun(@(x) isequal(x, [2386,2594,2828,3036,3240]), predict_time))
                    %find(cellfun(@(x) isequal(x, [2386,2594,2828,3036,3240,3441]), predict_time))
                end
                if is_finish %如果结束了，就直接退出，不跟新predict_time，因为没有pnr
                    break;
                end
                % 根据PNR筛选predict_time，保留TOP5时刻
                [~,idx_PNR]=sort(cell2mat(predict_pnr),'descend');
                if length(predict_pnr)>6
                    predict_time(idx_PNR(6:end))=[];
                    predict_pnr(idx_PNR(6:end))=[];
                end
            end
            % 加入TOP1PNR组合
            if is_finish % 如果结束了，可能在最后一个窗口的中间就结束了，就更换下一个MU
                break;
            end
            if ~exist('predict_pnr')%如果不存在predict_pnr，但是存在predict_time，说明只有一个时刻
                finalres(MU).pulse=union(finalres(MU).pulse,predict_time{1}+it);
                is_finish=1;
                break;
            end
            [~,idx_PNR]=sort(cell2mat(predict_pnr),'descend');
            finalres(MU).pulse=union(finalres(MU).pulse,predict_time{idx_PNR(1)}+it);
        end
    else %如果目前的结果存在较多的缺失时刻，需要尽可能补全
        last_time_num=length(finalres(MU).pulse);
        while 1 %finalres(MU).pulse(end)< end_point(sub,Task)-2000 %加这个while的目的是为了让大于2k的gap小于2k后得到更准确的迭代
            sequence=finalres(MU).pulse;
            subsets = {};
            currentSubset = sequence(1);
            % 循环遍历序列
            diff_se=diff(sequence);
            med=1.3*min(diff_se);% 保证放电时刻的不太大也不太小
            for i = 2:length(sequence)
                if sequence(i) - sequence(i-1) <= 1.5*med
                    % 如果差值小于等于350，则将当前元素加入当前子集
                    currentSubset = [currentSubset, sequence(i)];
                else
                    % 如果差值大于350，则结束当前子集并开始一个新的子集
                    subsets{end+1} = currentSubset;
                    currentSubset = sequence(i);
                end
            end
            % 结束循环后，将最后一个子集添加到结果中
            subsets{end+1} = currentSubset;
            t_diff=[];
            for i=1:length(subsets)
                if length(subsets{i})>=2
                    t_diff=[t_diff,diff(subsets{i})];
                end
            end
            predict_gap=round(mean(t_diff));
            
            for i=2:length(subsets)
                %如果识别到了几个时刻，但是不多，gap小于2.5k的先迭代中间的
                if subsets{i}(1)-subsets{i-1}(end)<2500
                    while subsets{i}(1)-subsets{i-1}(end)>2.5*predict_gap && subsets{i}(1)-subsets{i-1}(end)<2000%如果间隔没那么远，小于2.5k，先从2头往中间找时刻，直到2头的时刻比较接近或者找不到时刻了
                        it=subsets{i-1}(end)-2500;
                        it=max(it,start_pt+100); 
                        if it+window_len<end_time
                            signal=semg2sign(sEMG,it-99,it+window_len);
                        else
                            it=end_time-window_len;
                            signal=semg2sign(sEMG,it-99,it+window_len);
                        end
                        zn_exten15wgn=zn_exten15wgn_all(:,it-init_time:it+window_len-init_time);
                        save_all_ctj=save_all_ctj_all(it-init_time:it+window_len-init_time,it-init_time:it+window_len-init_time);
                        % 根据MU确定已知时刻
                        temp_pulse=finalres(MU).pulse-it;  
                        temp_pulse(temp_pulse<1)=[];
                        temp_pulse(temp_pulse>window_len)=[];
                        %tjns1= sum(save_all_ctj(temp_pulse,:),1)/length(temp_pulse);
                        %plot(tjns1)
                        left_time=subsets{i-1}(end)-it;
                        right_time=subsets{i}(1)-it;
                        if ~isempty(large_mu_time)
                            for l2=1:length(large_mu_time)
                                pulse=large_mu_time{l2}-it;
                                pulse(pulse<10 | pulse>window_len-100)=[];
                                temp_large_mu_time{l2}=pulse;
                            end
                            temp_large_mu_time = temp_large_mu_time(~cellfun('isempty', temp_large_mu_time));
                        end
                        %补全中间的时刻
                        [xreal,xreal2,~] = find_middle_2time(temp_pulse,save_all_ctj,predict_gap,...
                            signal(:,101:end),temp_large_mu_time,zn_exten15wgn,left_time,right_time,window_len);
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
                elseif subsets{i}(1)-subsets{i-1}(end)>=2500%如果间隔太远了，就只能从一头慢慢往后找时刻
                    %如果先往后迭代，一旦小于2.5k，结束迭代，准备往中间迭代
                    %要更新subsets
                    while subsets{i}(1)-subsets{i-1}(end)>=2500
                        it=subsets{i-1}(end)-(window_len-2500);
                        it=max(it,start_pt+100); 
                        is_finish=false;
                        clear predict_time predict_pnr
                        signal=semg2sign(sEMG,it-99,it+window_len); 
                        zn_exten15wgn=zn_exten15wgn_all(:,it-init_time:it+window_len-init_time);
                        save_all_ctj=save_all_ctj_all(it-init_time:it+window_len-init_time,it-init_time:it+window_len-init_time);
                        % 根据MU确定已知时刻
                        temp_pulse=finalres(MU).pulse-it;  
                        temp_pulse(temp_pulse<1)=[];
                        temp_pulse(temp_pulse>window_len-100)=[];
                        %tjns1= sum(save_all_ctj(temp_pulse,:),1)/length(temp_pulse);
                        %plot(tjns1)
                        %true_firing=firing_time{MU, 1}(1:end)-it;  
                        if ~isempty(large_mu_time)
                            for l2=1:length(large_mu_time)
                                pulse=large_mu_time{l2}-it;
                                pulse(pulse<10 | pulse>window_len-100)=[];
                                temp_large_mu_time{l2}=pulse;
                            end
                            temp_large_mu_time = temp_large_mu_time(~cellfun('isempty', temp_large_mu_time));
                        end
                        [xreal,~] = find_follow_time(temp_pulse,save_all_ctj,...
                            predict_gap,signal(:,101:end),temp_large_mu_time,zn_exten15wgn,window_len,subsets{i-1}(end)-it+25);
                        %%%% 这个点也是2525，之所以要从这个点开始搜索，因为temp_pulse=finalres(MU).pulse-it，其最后一个不在中间
                        if isempty(xreal{1})
                            is_finish=1;
                            break;
                        end
                        predict_time1=xreal{1};
                        predict_time=num2cell(predict_time1);
                        ttt=union(temp_pulse,predict_time{1});
                        UU_iter=0;
                        while predict_time{1}(end)<window_len-600 && is_finish==0
                            UU_iter=UU_iter+1;
                            while length(predict_time)<50 && predict_time{1}(end)<window_len-600
                                clear t_p t_pnr
                                for lp1=1:length(predict_time)
                                    if predict_time{lp1}(end)<window_len-600
                                        temp_pulse1=union(temp_pulse,predict_time{lp1});
                                        [xreal,Pnr] = find_follow_time(temp_pulse1,save_all_ctj,...
                                            predict_gap,signal(:,101:end),large_mu_time,zn_exten15wgn,window_len,predict_time{lp1});
                            %             toRemove = cellfun(@(x) length(x) == 1, xreal);
                            %             xreal(toRemove) = [];
                                        for j=1:length(xreal{1})
                                            t_p{lp1}{j}=[predict_time{lp1},xreal{1}(j)];
                                            t_pnr{lp1}{j}=Pnr(j);
                                        end
                                    end
                                end
                                % 如果t_p不存在，则说明没有新时刻了
                                if ~exist('t_p','var')
                                    is_finish=1;
                                    break;
                                end
                                % 用predict_time存储铺平后的t_p：这一轮的所有时刻组合
                                flattenedArray = cellfun(@(x) x(:), t_p, 'UniformOutput', false);
                                predict_time = cat(1, flattenedArray{:});
                                flattenedArray = cellfun(@(x) x(:), t_pnr, 'UniformOutput', false);
                                predict_pnr = cat(1, flattenedArray{:});
                                %find(cellfun(@(x) isequal(x, [2386,2594,2828,3036,3240]), predict_time))
                                %find(cellfun(@(x) isequal(x, [2386,2594,2828,3036,3240,3441]), predict_time))
                            end
                            if is_finish %如果结束了，就直接退出，不跟新predict_time，因为没有pnr
                                break;
                            end
                            % 根据PNR筛选predict_time，保留TOP5时刻
                            [~,idx_PNR]=sort(cell2mat(predict_pnr),'descend');
                            if length(predict_pnr)>6
                                predict_time(idx_PNR(6:end))=[];
                                predict_pnr(idx_PNR(6:end))=[];
                            end
                        end
                        % 加入TOP1PNR组合
                        if is_finish % 如果结束了，可能在最后一个窗口的中间就结束了，就更换下一个MU
                            break;
                        end
                        if ~exist('predict_pnr')%如果不存在predict_pnr，但是存在predict_time，说明只有一个时刻
                            subsets{i-1}=union(subsets{i-1},predict_time{1}-25+it);
                            finalres(MU).pulse=union(finalres(MU).pulse,predict_time{1}+it);
                            is_finish=1;
                            break;
                        end
                        [~,idx_PNR]=sort(cell2mat(predict_pnr),'descend');
                        subsets{i-1}=union(subsets{i-1},predict_time{idx_PNR(1)}+it);
                        finalres(MU).pulse=union(finalres(MU).pulse,predict_time{idx_PNR(1)}+it);
                    end
                end
            end
            if last_time_num == length(finalres(MU).pulse)
                break
            end
            last_time_num = length(finalres(MU).pulse);
        end
    end
    MU
end
%save_add=strcat('D:\dynamic_decompose\real data\',str{sub},'\',str{sub},'_Dynamic_Task',num2str(Task),'_Trial',num2str(trial),'_res_sys4_up1.mat');
%%
signal=semg2sign(sEMG,init_time-99,end_time); 
[y,Cxn_exten_inv]=get_zn2(signal);
zn_exten15wgn=get_zn(signal);
save_all_ctj=zn_exten15wgn'*zn_exten15wgn;
for i=1:length(finalres)
    [predict_gap,predict_cov]=precise_gap(finalres(i).pulse);
    finalres(i).mean_firing_gap=predict_gap;
    finalres(i).cov=predict_cov;
    finalres(i).mean_firing_rate=2048/finalres(i).mean_firing_gap;
    finalres(i).firing_time=finalres(i).pulse/2048;
    finalres(i).Cst=mean(y(:,finalres(i).pulse-(init_time-99-1)),2)';
    ctsj=mean(save_all_ctj(finalres(i).pulse-init_time,:),1);
    finalres(i).ctsj=ctsj;
end
%save(save_add,'finalres');
%% 若MUAP形状正常，根据MUAP形状判断MU质量
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
%     % 将这些数从数组中移除
%     subplot(2, 1, 2); % 创建一个新的图形窗口
%     plot(finalres(llll).ctsj, 'b'); % 'b' 表示蓝色线条
%     tallmuap(indices) = [];
%     tstd=std(tallmuap);
%     gap/tstd
%     MUAP_ratio(llll)=gap/tstd;
% end
% binary_jd=MUAP_ratio>13;
% finalres=finalres(binary_jd);
end

%% 分解+验证
function zn_exten15wgn=Whiten(xn_exten15wgn)
    %用来白化矩阵
    for i =1:size(xn_exten15wgn,1)
        xn_exten15wgn(i,:)=xn_exten15wgn(i,:)-mean(xn_exten15wgn(i,:));
    end
    Cxn_exten=xn_exten15wgn*xn_exten15wgn'/size(xn_exten15wgn,2);
%     [U, D] = eig(Cxn_exten);
%     W=U*diag(1./sqrt(diag(D)))*U';
%     zn_exten15wgn=W*xn_exten15wgn;
    [u,s,v]=svd(Cxn_exten);
    W=v*diag(1./sqrt(diag(s)))*v';
    zn_exten15wgn=W*xn_exten15wgn;
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
function zn_exten15wgn=get_zn(signal)
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
    zn_exten15wgn=Whiten(xn_exten15wgn(:,101:end));
end
function hat_t_j_n_a_n = CalculateHatT_j_n_a_n(c_t_j_y_n_a, y,timeset)
    % c_t_j_y_n_a: 已知的 c_t_j_y_n_a 向量
    % y: 信号向量
    % n: 要计算的时间点 n
    % a1, a2: 时间范围 [a1, a2]
    
    y_matrix= y(:,timeset) * y(:,timeset)';
%     for p = timeset
%         y_matrix= y_matrix+y(:,p) * y(:,p)';
%     end
    % 计算 \hat{t}_{j,n_a}(n)
    hat_t_j_n_a_n = c_t_j_y_n_a' * (y_matrix\y(:,timeset)); 
end
function [ixvn1set,tjns1] = Split_spike(tjns1,y,timeset)
% 这个函数要怎么写，才能达到一个比较O,K的结果呢
% 直接用CKC的好了
Nnloop=7;
y_matrix= y(:,timeset) * y(:,timeset)';
tmp=y_matrix\y(:,timeset);
for j=1:Nnloop
  if j==1
         Nfind=8;%6000-14 4000-11 12000-24
  else
         Nfind=Nfind+1;%每步增加1个点*****************************************
  end
vn_1=tjns1;
clear ixvn1set;
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
    ixvn1set=dele_adjacent_time(ixvn1set,tjns1,120);
    if length(ixvn1set)>=Nfind
            break; 
    end
end
tjns1 = mean(y(:,timeset(ixvn1set)),2)' * (tmp); 
%plot(tjns1)
end
end
function [Pnr] = pnr(tjn1,tsnj)
    Pnr=10*log(mean(tjn1.*tjn1)/mean(tsnj.*tsnj));
end
function [predict_gap,predict_cov] = precise_gap(pulse)
sequence=pulse;
subsets = {};
currentSubset = sequence(1);
% 循环遍历序列
diff_se=diff(sequence);
med=1.3*min(diff_se);% 保证放电时刻的不太大也不太小
for i = 2:length(sequence)
if sequence(i) - sequence(i-1) <=  1.5*med
    % 如果差值小于等于350，则将当前元素加入当前子集
    currentSubset = [currentSubset, sequence(i)];
else
    % 如果差值大于350，则结束当前子集并开始一个新的子集
    subsets{end+1} = currentSubset;
    currentSubset = sequence(i);
end
end
% 结束循环后，将最后一个子集添加到结果中
subsets{end+1} = currentSubset;
t_diff=[];
for i=1:length(subsets)
if length(subsets{i})>=2
    t_diff=[t_diff,diff(subsets{i})];
end
end
predict_gap=round(mean(t_diff));
predict_cov=std(t_diff)/predict_gap;
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