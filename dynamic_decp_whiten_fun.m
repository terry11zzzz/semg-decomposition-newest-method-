function [pulsefinal,union_pulsetime,union_set,tjns1set] = dynamic_decp_whiten_fun(xn_exten15wgn,gap,para)
nsamples=size(xn_exten15wgn,2);
iteration_num = gap+1;%迭代次数
numn1=1;
% par=4;
clear n0set;
tjns1set=cell(iteration_num-1,1);
Rset=cell(iteration_num-1,1);
min_gap=200;
gap_n0=nsamples/gap;%n0的间隔
tic;
for g=1:nsamples
    lamuga(g)=xn_exten15wgn(:,g)'*xn_exten15wgn(:,g);%选出的加入15db的高斯白噪声，必须加入噪声，不然无法算出n0
    %     lamuga(g)=xn_exten15wgn(:,g)'* Cxn_exten_trun_inv*xn_exten15wgn(:,g);
end
save_all_ctj=xn_exten15wgn'*xn_exten15wgn;
[~,index2]=sort(lamuga);
shared_ctj = parallel.pool.Constant(save_all_ctj);  % 每个 worker 共享一份
parfor UU=1:gap
    %n0=6400-UU;
    %[~,index2]=sort(lamuga);%排序后，sA是排序好的向量，index 是 向量sA 中对 A 的索引
    %n0=index2(round(nsamples/2));%取一个序列的中位数
    n0=index2(round(gap_n0*UU));
    %n0=2*UU;%index2(nsamples-UU+1);%取一个序列的最大值
    %n0=2*UU;%由于10个延时，事实上连续4个时刻对应的MU大致相同，并且能够极大减少计算复杂度
    %n0set(UU) =n0;
    vn0=shared_ctj.Value(n0,:);%这个序列代表着，通过第n0时刻，算出该时刻MU发放的序列
    vn1set=vn0;
    ele_num=100;
    [~,idx]=sort(vn0);
    Gn0=idx(end-ele_num:end-1);
    for loop1=ele_num:-1:1
        vn1=shared_ctj.Value(Gn0(loop1),:);
        [~,idx_n1]=sort(vn1);
        Gn1=idx_n1(end-ele_num:end-1);
        if any(Gn1==n0)
            n1=Gn0(loop1);
            vn1set=vn1.*vn0;
            break;
        end
    end
    Gn0Gn1=intersect(Gn1,Gn0);
    for loop1=length(Gn0Gn1):-1:1
        vn2=shared_ctj.Value(Gn0Gn1(loop1),:);
        [~,idx_n2]=sort(vn2);
        Gn2=idx_n2(end-ele_num:end-1);
        if any(Gn2==n0)&&any(Gn2==n1)
            vn1set=(vn1+vn0+vn2)/3;
            break;
        end
    end
    %根据初始时刻n0，找到峰值时刻n1
    %STEP1 OVER
    %================================================
    
    
    
    %================================================
    %STEP2 迭代n0时刻对应的MU的发放时刻，并且不delete，增加改MU发放信息与能量
    for i=1:numn1
     for j=1:para(1)  %6000-15 4000-10 12000-28
        Nfind=1*j;%每步增加4个点*****************************************
        vn_1(i,:)=vn1set(i,:);
        %clear ixvn1set;
        ixvn1set=[];
        MUAP_ratio=20;
        
        while 1
            [~,ixvn1]=max(vn_1);
%             if length(ixvn1set)>5
%                 pulse=union(ixvn1set,ixvn1);
%                 [MUAPSet]=CalMUAPNonShift(size(signal,2),signal,pulse,[],size(signal,1));
%                 tallmuap=[];
%                 for jjj=1:size(signal,1)
%                     tallmuap=union(tallmuap,MUAPSet{jjj});
%                     %plot(MUAPSet{jjj,llll});hold on
%                 end
%                 gap=max(tallmuap)-min(tallmuap);
%                 [~, indices] = maxk(abs(tallmuap), 100);
%                 % 将这些数从数组中移除
%                 tallmuap(indices) = [];
%                 tstd=std(tallmuap);
%                 MUAP_ratio=gap/tstd;
%             end
            if ~any(ixvn1set==ixvn1) && MUAP_ratio>15%若不在集合里面 或者 一开始ixvn1set为空
                nvr=ixvn1;
                vn_1(ixvn1)=-inf;
                %vn_1(ixvn1-100:ixvn1+100)=-inf;
                ixvn1set=[ixvn1set,nvr];%ixvn1set就是这一次最大点的时刻集合
                ixvn1set=sort(ixvn1set);
                if length(dele_adjacent_time(ixvn1set,vn1set,min_gap))>=Nfind
                        break; 
                end
            else
                vn_1(ixvn1)=-inf;
            end
        end
        
            %如果时刻的个数在10到15之间，那就结束对    r个时刻产生的序列共同大于阈值的时刻集合的寻找
        vn1set(i,:)=sum(shared_ctj.Value(ixvn1set,:),1)/Nfind;   
%         if length(ixvn1set) >200 %如果连续时刻太多了，超过300个，说实话这个时刻质量肯定不是很好，直接跳过就好
%             break;
%         end
       % plot(vn1set)
     end
    end
     %STEP2 OVER
    %==================================================
    R=dele_adjacent_time(ixvn1set,vn1set(i,:),min_gap);
    tjns1= sum(shared_ctj.Value(R,:),1)/length(R) ;
    %共同发放时刻寻找完毕，结果是R'
    %STEP3 OVER
%     if length(ixvn1set) <200 %如果小于300，说明质量还行，那就11个
%         eq=1;
%         Nnloop=11;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     else
%         eq=0;
         Nnloop=para(2);%%6000-7 4000-8 12000-14
%     end
    %%%仿真数据70，真实数据10(序列长度的1/10)+2*23-12547/100 的70%到1/3 
    %%% 真实数据5(序列长度的1/10)+2*13-10831/200 的70%到1/3 
    %%6000-7 4000-5 12000-14
    %STEP5 OVER,得到的一次大循环后的第jMU的发放序列TJNS1
    %============================
    %%%通过
    %%%之前得到的MU时刻以及发放序列，进行delete迭代，这里会担心得到混杂的时刻，如果有噪声，其他比该MU能量大的MU的信息没有排除干净
    %%%就会导致混杂，没有delete的话，这样混杂的时刻会少一点，因为临近时刻会把该MU的信息相比于deletetime更凸显

     for j=1:Nnloop
          if j==1
              %if eq==1
                 Nfind=para(3);%6000-14 4000-8 12000-24
%               else %如果序列质量不咋地就别折磨他了
%                   Nfind=14;
%               end
          else
                 Nfind=Nfind+1;%每步增加1个点*****************************************
          end
        vn_1=tjns1;
        %clear ixvn1set;
        ixvn1set=[];
        MUAP_ratio=20;
        %这里每次都从头开始迭代实在太慢了，我得加快迭代速度
        while 1

            [~,ixvn1]=max(vn_1);
            % 如果去掉邻近时刻领导MUAP-ratio满足要求
            if ~any(ixvn1set==ixvn1) && MUAP_ratio>15%若不在集合里面 或者 一开始ixvn1set为空
                nvr=ixvn1;
                %将周围100的时刻都变成-inf
                %vn_1(ixvn1-100:ixvn1+100)=-inf;
                vn_1(ixvn1)=-inf;
                ixvn1set=[ixvn1set,nvr];%ixvn1set就是这一次最大点的时刻集合
                ixvn1set=sort(ixvn1set);
                ixvn1set=dele_adjacent_time(ixvn1set,tjns1,min_gap);
                if length(ixvn1set)>=Nfind
                        break; 
                end
            else
                vn_1(ixvn1)=-inf;
            end
            
        end
        tjns1=sum(shared_ctj.Value(ixvn1set,:),1)/length(ixvn1set) ;
        %plot(tjns1)
     end
    %这一行是用来
    tjns1=sum(shared_ctj.Value(ixvn1set,:),1)/length(ixvn1set);
    tjns1set{UU,1}=tjns1;%这个是得到的一次大循环后的第jMU的发放序列
    %STEP7 OVER   
    UU
end
%%
tspend=toc;
tic;
clear pulsefinal;
clear pulsetimef;
clear pulsetimefset;
pulsetimefset=struct;%定义结构体
segthr=[2,5];%设置3个值，有利于多求，准确求,限制了找最小标准差的阈值的个数
n=1;%不初始化，三次结果pulsetimefset放一起,一长条
        %======================================
union_set=[];
union_set1=[];
union_pulsetime=[];
%for runtime=1:length(segthr)
    for m=1:size(tjns1set,1)
         %step1:得到一个序列tjns1set的最大值以及幅度平均值，作为判断峰值时刻的阈值
        stjn=tjns1set{m,1};
        nummaxs=para(4);%找15次最大值
        
        stjns=stjn;
        maxsset=zeros(1,nummaxs);
        for j=1:nummaxs%每个信号找15次最大值
            [maxs,maxaips]=max(stjns);
            maxsset(1,j)=maxs;%存放每行信号最大值
            stjns(maxaips)=-inf;
        end
        avermaxa=sum(maxsset)/nummaxs;%求幅度平均值
        %step1 OVER 最大值=maxsset，幅度平均值=avermaxa
        %======================================
        
        %======================================
        %step2:从大到小设置一系列阈值，分别求每一个阈值找到所有时刻间隔的平均值和标准差
        sigtimdifstdset=zeros(1,61);
        for i=100:-1:40
            thrtim=(i/100)*avermaxa;
            sigpulset=find(stjn>=thrtim);%sigtim是比均值大的（峰值）序列时刻
            
            if length(sigpulset)<7%如果大于阈值的长度小于10
                sigtimdifstdset(1,101-i)=inf;
                continue;
            else
                 %%%%%%%%%感觉这里还需要加一个STEP，需要序列间隔大于100才行
                 %如果 sigtim 存在间隔为小于10的，就记录选其中数值小的那个序号，目的是在下一次循环的时候直接去掉
                while 1
                    sigtimdiff=diff(sigpulset);
                    temp_record=find(sigtimdiff<=100);%%%这个根据实际信号的发放间隔更改，原则是间隔越大，阈值越大
                    if ~isempty(temp_record)%如果存在间隔在3个以内的时刻，去掉小的那一个
                       for j=1:length(temp_record)

                           if stjn(sigpulset(temp_record(j)))>=stjn(sigpulset(temp_record(j)+1))
                               %说明前面一个时刻的发放值大,去掉后面那个
                               sigpulset(temp_record(j)+1)=-1;
                               sigpulset=sigpulset(sigpulset>0);
                               break;
                           else%说明后面一个时刻的发放值大,去掉前面那个
                               sigpulset(temp_record(j))=-1;
                               sigpulset=sigpulset(sigpulset>0);
                               break;
                           end
                       end
                    else
                        break;
                    end
                end
                if i<50 %这个限制条件主要是防止：因为标准差过小导致，阈值过高，导致序列过少
                    sigtimdiff=diff(sigpulset);
                    %%%%%%%%%%标准变成COV
                    sigtimdifstd=std(sigtimdiff);%求峰值间隔的标准差
                    sigtimdifcov=sigtimdifstd/mean(sigtimdiff);
                    sigtimdifstdset(1,101-i)=sigtimdifcov;%记录得到的峰值间隔的标准差
                else
                    sigtimdifstdset(1,101-i)=inf;
                end

                
             end
        end
       %step2 OVER 所有时刻间隔的平均值=sigtimdifmeanset，所有时刻间隔的标准差=sigtimdifstdset
        %======================================
       
        %======================================
        %step3:找到合适的 时刻间隔的标准差，来确定对应的序列发放值的阈值，这个阈值就是序列tjns1set的峰值时刻的阈值
        %比如通过阈值thrtim=(i/100)*avermaxa;，得到了一些大于这个阈值的峰值时刻，而当这些峰值时刻间隔的标准差最小时候
        %这意味着峰值时刻间隔的分布比较均匀，然后调整阈值，导致了阈值在一定范围内变化的时候，标准差不变，
        %也就说明了在这段阈值内，得到的峰值时刻是不变的，当调整多次阈值都不变，说明阈值thrtim=(i/100)*avermaxa;就应该时我们想要的阈值了
        sigtimdifstdset1=sigtimdifstdset;
        [stdmin,i_best]=min(sigtimdifstdset1(sigtimdifstdset1~=0));%找到标准差的最小值
        thrtim_best=(1.01-i_best/100)*avermaxa;

        sigpulset=find(stjn>=thrtim_best);%比如pulsetime0：1x34，就说明有34个大于阈值的峰值时刻
                while 1
                sigtimdiff=diff(sigpulset);
                temp_record=find(sigtimdiff<=100);
                    if ~isempty(temp_record)%如果存在间隔在3个以内的时刻，去掉小的那一个
                       for j=1:length(temp_record)

                           if stjn(sigpulset(temp_record(j)))>=stjn(sigpulset(temp_record(j)+1))
                               %说明前面一个时刻的发放值大,去掉后面那个
                               sigpulset(temp_record(j)+1)=-1;
                               sigpulset=sigpulset(sigpulset>0);
                               break;
                           else%说明后面一个时刻的发放值大,去掉前面那个
                               sigpulset(temp_record(j))=-1;
                               sigpulset=sigpulset(sigpulset>0);
                               break;
                           end
                       end
                    else
                        break;
                    end
                end
        pulsetime0= sigpulset;
        %STEP3 OVER
        %已经找到了最小标准差对应的阈值，且这个阈值在一段比较宽的区间内，标准差都是最小的，得到的峰值时刻是不变，且应该是要求的时刻
        %阈值thrtim_best，最小标准差stdmin，最小标准差对应的个数secstdmin
        %求出基于最佳阈值的峰值时刻pulsetime0
        %================================================
        %STEP4 去掉标准差很大的序列，调整标准差较大的序列（因为还能抢救一下）
        %如果求出的标准差【stdmin】太大，则需要去掉一些发放间距不合理的脉冲**********，注意可能需要调节参数
%          if stdmin>140 %去掉标准差太大的一个IPT
%             continue;
%          end
        %减小峰值间隔的标准差，目的是，让这个峰值间隔都在一个比较小的范围里面波动
        %每连续k个峰值间隔，计算他的最小间隔以及均值，
        %找到连续k个间隔中，标准差特别大的，应该就是那个不太对的峰值时刻了，把这个不太对的峰值时刻去掉，剩下的峰值时刻的标准差就会比较小
        if stdmin<=0.8%******************如果信号不好，就用0.4
            pulsetimef=pulsetime0;
            union_set1=union(union_set1,Rset{m,1});
            pulsetimefdiff=diff(pulsetimef);%记录得到最终序列的时刻间隔序列
            pulsetimefstd=std(pulsetimefdiff);%%记录得到最终序列的时刻间隔的标准差
            %if pulsetimefstd<70 %*********最终合乎要求的是标准差【stdmin】不能太大的，标准差大的可以手工操作
            pulsetimefset(1,n).pulse=pulsetimef;
            pulsetimefset(1,n).std=pulsetimefstd/mean(pulsetimefdiff);
            pulsetimefset(1,n).m=m;%找到所求原始脉冲tjns1set的对应位置
            union_set=union(union_set,Rset{m,1});
            %union_pulsetime=union(union_pulsetime,pulsetimef);
            n=n+1;
        end
           
    end
    toc
%end       
%%
clear pulseclass
pulseclass=struct;
pulsefinal=struct;
i=1;
iw=1;
IDallset=(1:length(pulsetimefset));%初始化为全部的ID
tic
if ~isfield(pulsetimefset,'pulse')
    return;
end
if length(IDallset)==1
    pulsefinal(1).pulse=pulsetimefset(1,1).pulse;
    return;
end
pulsetimefset1=pulsetimefset;
        %step1归类，认为发放序列重叠的发放时刻大于0.5，就是一个MU发放的
    %  每一个序列类别的 pulseclass，pulseclass 的第1类元素就是  pulsetimefset 中的属于第一类mu的下标：pulsetimefset （pulseclass ）
    while 1  %执行若干次归类操作
        iq=1;
        clear pulseclass(iw).ID;
        pulseclass(iw).ID=[];%存放每一类的ID（号码）
        array=zeros(1,IDallset(end)-1);
        parfor iq=1:length(IDallset)-1   %执行每一类的归类操作
            %寻找IDallset(i)的发放序列（时刻在+-50范围内）与之后的发放序列的共同时刻，，并且将找到的时刻加入这个聚类中
            %
                pu1 = zeros(nsamples,1);
                pu2 = zeros(nsamples,1);
                pu1(pulsetimefset(1,IDallset(1)).pulse)=1;% 因为已经被分类的id会被del，所以一直是1
                temp1=pulsetimefset(1,IDallset(1+iq)).pulse+1;
                temp1(temp1>nsamples)=[];
                pu2(temp1)=1;
                temp2=pulsetimefset(1,IDallset(i+iq)).pulse-1;
                temp2(temp2<1)=[];
                pu2(temp2)=1;
                pu2(pulsetimefset(1,IDallset(i+iq)).pulse)=1;
                xcpul=xcorr(pu1,pu2,50);
                [maxxc,ipmaxxc]=max(xcpul);

                %pulseinterleng=maxxc;
                tp1=maxxc/length(pulsetimefset(1,IDallset(i)).pulse);
                tp2=maxxc/length(pulsetimefset(1,IDallset(i+iq)).pulse);
                if tp1>=0.3&&tp2>=0.30
                    %大MU的IPT近似度，很容易大于0.5.但是其中包含了小MU信息，如果设置太低，那些包含不少小MU信息的IPT就被舍弃，但是速度很快。
                    %但是被覆盖的小MU之后也很难搞出来，所以还是舍弃吧。如果需要提取这部分MU了，那就降低速度，把这部分阈值设置高一点。
                    %当前还是设置低一点，把垃圾的IPT舍弃算了。这些垃圾IPT反正在后面也会被舍弃（重复次数小于2的会被cut）
                                %pulseclass(iw).ID=[pulseclass(iw).ID,IDallset(i+iq)];
                                array(iq)=IDallset(i+iq);%将找到的时刻加入这个聚类中  pulseclass(iw).ID
                                %pulsetimefset1(1,IDallset(i+iq)).pulse=pulsetimefset(1,IDallset(i+iq)).pulse+(ipmaxxc-50-1);%ipmaxxc-50-1

                end
           % end
            %pulseclass(iw).ID就是同一mu的聚类
%             if IDallset(i+iq)==IDallset(end)
%                 break;
%             end
          %  iq=iq+1;
        end
        if ~any(array)
             pulseclass(iw).ID=IDallset(i);
        else
            pulseclass(iw).ID=[nonzeros(array)',IDallset(i)];
        end
        NewIDallset=setdiff(IDallset,pulseclass(iw).ID);%去掉已经归类的ID,最后可以检查IDallset,看是否全部归类
        iw=iw+1;
        IDallset=NewIDallset;
        if isempty(IDallset)||length(IDallset)<3%防止单独，没有同组情况的存在
            break;
        end
    end
    toc
    %step1归类 over pulseclass 是归类的结果，并将每个时刻对齐后的序列记录为pulsetimefset.pulse1
    %====================================================
    %%
    %====================================================
    %step2  对时刻较多的类，选择一个标准差小的序列作为 pulsefinal
    %%对时刻较少的类，因为这种MU对应的序列就算是标准差很小，只能说明序列包含的时刻比较准确
    %但是由于个数比较少，很可能序列的发放是有漏的,所以想把是吧这些序列对齐后，
    %再union所有对齐后的序列，这样尽可能的将此MU的发放时刻都结合起来，再排除这个类中与其他类（union或者finalpulse） 重合度大于0.3的交集时刻，将这个步骤叫做去重
    %最后对union去重后的MU处理，加入 pulsefinal


    %区分类中元素个数是否大于10，
    for i=1:length(pulseclass)
       if length(pulseclass(i).ID)<2
           pulseclass(i).judge =1;
       else
           union_pulsetime=[];
           pulseclass(i).judge =0;
           for j=1:length( pulseclass(i).ID)
             union_pulsetime=union(union_pulsetime,pulsetimefset(pulseclass(i).ID(j)).pulse);
             pulseclass(i).union=union_pulsetime;
           end
       end
    end
    % pulseclass(i).judge 如果是1，那么就说明序列的个数小于10

    %如果小于5就union（上一步已经对齐过）
    for i=1:length(pulseclass)
        all_union=[];
        if pulseclass(i).judge==0
            %pulseclass(i).union=all_union;
            continue;
        end 
        for j=1:length( pulseclass(i).ID)
           all_union=union(all_union,pulsetimefset(pulseclass(i).ID(j)).pulse) ;
        end
        pulseclass(i).union=all_union;
    end


    %如果大于30就选择一个标准差小的序列作为 pulsefinal（上一步已经对齐过）
    %选每类里面标准差最小/最多的一个（最多比最小好）,并且数目要多余一个的
    classstdset=struct;
    pulsefinal=struct;
    stdsetunique=struct;
    classstdip=struct;
    % unifindlen=struct;
    selstd=zeros(2,length(pulseclass));
    for k=1:length(pulseclass)
        if pulseclass(k).judge==1
            continue;
        end
        clear classstdset(k).std;
        classstdset(k).std=[];
        classstdset(k).std=[classstdset(k).std,pulsetimefset(1,pulseclass(k).ID).std];%收集每一类的标准差，结合空集，相当于多次循环添加元素
        stdsetunique(k).std=unique(classstdset(k).std);
        unifindlen=zeros(1,length(stdsetunique(k).std));%一定要初始化，不然下次运行可能出现错误
        %如果每一类里面的标准差有重复的，则选择多次重复的标准差，这样准确率高，如果上面运行四次，重复的要多于四次为好
        if length(stdsetunique(k).std)<length(classstdset(k).std)
            for ie=1:length(stdsetunique(k).std)
                unifind=find(classstdset(k).std==stdsetunique(k).std(ie));
                unifindlen(ie)=length(unifind);
            end
            [~,unilenmaxip]=max(unifindlen);%得到标准差最多的那个标准差的序列
            classstdip(k).ip=find(classstdset(k).std==stdsetunique(k).std(unilenmaxip));
            selstd(1,k)=stdsetunique(k).std(unilenmaxip);
            selstd(2,k)= classstdip(k).ip(1);%可随机选，都是一样的，因为这样多的保留小数位精度，间隔需要一点不差才能相等
        else  %没有重复的标准差就选择最小的
            [selstd(1,k),selstd(2,k)]=min(classstdset(k).std);%存放最小标准差的位置*****最小的有问题，用最多的试试看
            classstdip(k).ip=selstd(2,k);
        end
        pulsefinal(k).pulse=pulsetimefset(1,pulseclass(k).ID(selstd(2,k))).pulse;%找出每一类最小标准差所对应的发放，作为最终要分解的MU
        pulsefinal(k).unionpulsetime=pulseclass(k).union;
    end

    

timespend=toc;
timespend=timespend/59

end

