function [pulsefinal,union_pulsetime,union_set,tjns1set] = dynamic_decp_whiten_fun(xn_exten15wgn,gap,para)
nsamples=size(xn_exten15wgn,2);
iteration_num = gap+1;%��������
numn1=1;
% par=4;
clear n0set;
tjns1set=cell(iteration_num-1,1);
Rset=cell(iteration_num-1,1);
min_gap=200;
gap_n0=nsamples/gap;%n0�ļ��
tic;
for g=1:nsamples
    lamuga(g)=xn_exten15wgn(:,g)'*xn_exten15wgn(:,g);%ѡ���ļ���15db�ĸ�˹�����������������������Ȼ�޷����n0
    %     lamuga(g)=xn_exten15wgn(:,g)'* Cxn_exten_trun_inv*xn_exten15wgn(:,g);
end
save_all_ctj=xn_exten15wgn'*xn_exten15wgn;
[~,index2]=sort(lamuga);
shared_ctj = parallel.pool.Constant(save_all_ctj);  % ÿ�� worker ����һ��
parfor UU=1:gap
    %n0=6400-UU;
    %[~,index2]=sort(lamuga);%�����sA������õ�������index �� ����sA �ж� A ������
    %n0=index2(round(nsamples/2));%ȡһ�����е���λ��
    n0=index2(round(gap_n0*UU));
    %n0=2*UU;%index2(nsamples-UU+1);%ȡһ�����е����ֵ
    %n0=2*UU;%����10����ʱ����ʵ������4��ʱ�̶�Ӧ��MU������ͬ�������ܹ�������ټ��㸴�Ӷ�
    %n0set(UU) =n0;
    vn0=shared_ctj.Value(n0,:);%������д����ţ�ͨ����n0ʱ�̣������ʱ��MU���ŵ�����
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
    %���ݳ�ʼʱ��n0���ҵ���ֵʱ��n1
    %STEP1 OVER
    %================================================
    
    
    
    %================================================
    %STEP2 ����n0ʱ�̶�Ӧ��MU�ķ���ʱ�̣����Ҳ�delete�����Ӹ�MU������Ϣ������
    for i=1:numn1
     for j=1:para(1)  %6000-15 4000-10 12000-28
        Nfind=1*j;%ÿ������4����*****************************************
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
%                 % ����Щ�����������Ƴ�
%                 tallmuap(indices) = [];
%                 tstd=std(tallmuap);
%                 MUAP_ratio=gap/tstd;
%             end
            if ~any(ixvn1set==ixvn1) && MUAP_ratio>15%�����ڼ������� ���� һ��ʼixvn1setΪ��
                nvr=ixvn1;
                vn_1(ixvn1)=-inf;
                %vn_1(ixvn1-100:ixvn1+100)=-inf;
                ixvn1set=[ixvn1set,nvr];%ixvn1set������һ�������ʱ�̼���
                ixvn1set=sort(ixvn1set);
                if length(dele_adjacent_time(ixvn1set,vn1set,min_gap))>=Nfind
                        break; 
                end
            else
                vn_1(ixvn1)=-inf;
            end
        end
        
            %���ʱ�̵ĸ�����10��15֮�䣬�Ǿͽ�����    r��ʱ�̲��������й�ͬ������ֵ��ʱ�̼��ϵ�Ѱ��
        vn1set(i,:)=sum(shared_ctj.Value(ixvn1set,:),1)/Nfind;   
%         if length(ixvn1set) >200 %�������ʱ��̫���ˣ�����300����˵ʵ�����ʱ�������϶����Ǻܺã�ֱ�������ͺ�
%             break;
%         end
       % plot(vn1set)
     end
    end
     %STEP2 OVER
    %==================================================
    R=dele_adjacent_time(ixvn1set,vn1set(i,:),min_gap);
    tjns1= sum(shared_ctj.Value(R,:),1)/length(R) ;
    %��ͬ����ʱ��Ѱ����ϣ������R'
    %STEP3 OVER
%     if length(ixvn1set) <200 %���С��300��˵���������У��Ǿ�11��
%         eq=1;
%         Nnloop=11;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     else
%         eq=0;
         Nnloop=para(2);%%6000-7 4000-8 12000-14
%     end
    %%%��������70����ʵ����10(���г��ȵ�1/10)+2*23-12547/100 ��70%��1/3 
    %%% ��ʵ����5(���г��ȵ�1/10)+2*13-10831/200 ��70%��1/3 
    %%6000-7 4000-5 12000-14
    %STEP5 OVER,�õ���һ�δ�ѭ����ĵ�jMU�ķ�������TJNS1
    %============================
    %%%ͨ��
    %%%֮ǰ�õ���MUʱ���Լ��������У�����delete����������ᵣ�ĵõ����ӵ�ʱ�̣�����������������ȸ�MU�������MU����Ϣû���ų��ɾ�
    %%%�ͻᵼ�»��ӣ�û��delete�Ļ����������ӵ�ʱ�̻���һ�㣬��Ϊ�ٽ�ʱ�̻�Ѹ�MU����Ϣ�����deletetime��͹��

     for j=1:Nnloop
          if j==1
              %if eq==1
                 Nfind=para(3);%6000-14 4000-8 12000-24
%               else %�������������զ�ؾͱ���ĥ����
%                   Nfind=14;
%               end
          else
                 Nfind=Nfind+1;%ÿ������1����*****************************************
          end
        vn_1=tjns1;
        %clear ixvn1set;
        ixvn1set=[];
        MUAP_ratio=20;
        %����ÿ�ζ���ͷ��ʼ����ʵ��̫���ˣ��ҵüӿ�����ٶ�
        while 1

            [~,ixvn1]=max(vn_1);
            % ���ȥ���ڽ�ʱ���쵼MUAP-ratio����Ҫ��
            if ~any(ixvn1set==ixvn1) && MUAP_ratio>15%�����ڼ������� ���� һ��ʼixvn1setΪ��
                nvr=ixvn1;
                %����Χ100��ʱ�̶����-inf
                %vn_1(ixvn1-100:ixvn1+100)=-inf;
                vn_1(ixvn1)=-inf;
                ixvn1set=[ixvn1set,nvr];%ixvn1set������һ�������ʱ�̼���
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
    %��һ��������
    tjns1=sum(shared_ctj.Value(ixvn1set,:),1)/length(ixvn1set);
    tjns1set{UU,1}=tjns1;%����ǵõ���һ�δ�ѭ����ĵ�jMU�ķ�������
    %STEP7 OVER   
    UU
end
%%
tspend=toc;
tic;
clear pulsefinal;
clear pulsetimef;
clear pulsetimefset;
pulsetimefset=struct;%����ṹ��
segthr=[2,5];%����3��ֵ�������ڶ���׼ȷ��,����������С��׼�����ֵ�ĸ���
n=1;%����ʼ�������ν��pulsetimefset��һ��,һ����
        %======================================
union_set=[];
union_set1=[];
union_pulsetime=[];
%for runtime=1:length(segthr)
    for m=1:size(tjns1set,1)
         %step1:�õ�һ������tjns1set�����ֵ�Լ�����ƽ��ֵ����Ϊ�жϷ�ֵʱ�̵���ֵ
        stjn=tjns1set{m,1};
        nummaxs=para(4);%��15�����ֵ
        
        stjns=stjn;
        maxsset=zeros(1,nummaxs);
        for j=1:nummaxs%ÿ���ź���15�����ֵ
            [maxs,maxaips]=max(stjns);
            maxsset(1,j)=maxs;%���ÿ���ź����ֵ
            stjns(maxaips)=-inf;
        end
        avermaxa=sum(maxsset)/nummaxs;%�����ƽ��ֵ
        %step1 OVER ���ֵ=maxsset������ƽ��ֵ=avermaxa
        %======================================
        
        %======================================
        %step2:�Ӵ�С����һϵ����ֵ���ֱ���ÿһ����ֵ�ҵ�����ʱ�̼����ƽ��ֵ�ͱ�׼��
        sigtimdifstdset=zeros(1,61);
        for i=100:-1:40
            thrtim=(i/100)*avermaxa;
            sigpulset=find(stjn>=thrtim);%sigtim�ǱȾ�ֵ��ģ���ֵ������ʱ��
            
            if length(sigpulset)<7%���������ֵ�ĳ���С��10
                sigtimdifstdset(1,101-i)=inf;
                continue;
            else
                 %%%%%%%%%�о����ﻹ��Ҫ��һ��STEP����Ҫ���м������100����
                 %��� sigtim ���ڼ��ΪС��10�ģ��ͼ�¼ѡ������ֵС���Ǹ���ţ�Ŀ��������һ��ѭ����ʱ��ֱ��ȥ��
                while 1
                    sigtimdiff=diff(sigpulset);
                    temp_record=find(sigtimdiff<=100);%%%�������ʵ���źŵķ��ż�����ģ�ԭ���Ǽ��Խ����ֵԽ��
                    if ~isempty(temp_record)%������ڼ����3�����ڵ�ʱ�̣�ȥ��С����һ��
                       for j=1:length(temp_record)

                           if stjn(sigpulset(temp_record(j)))>=stjn(sigpulset(temp_record(j)+1))
                               %˵��ǰ��һ��ʱ�̵ķ���ֵ��,ȥ�������Ǹ�
                               sigpulset(temp_record(j)+1)=-1;
                               sigpulset=sigpulset(sigpulset>0);
                               break;
                           else%˵������һ��ʱ�̵ķ���ֵ��,ȥ��ǰ���Ǹ�
                               sigpulset(temp_record(j))=-1;
                               sigpulset=sigpulset(sigpulset>0);
                               break;
                           end
                       end
                    else
                        break;
                    end
                end
                if i<50 %�������������Ҫ�Ƿ�ֹ����Ϊ��׼���С���£���ֵ���ߣ��������й���
                    sigtimdiff=diff(sigpulset);
                    %%%%%%%%%%��׼���COV
                    sigtimdifstd=std(sigtimdiff);%���ֵ����ı�׼��
                    sigtimdifcov=sigtimdifstd/mean(sigtimdiff);
                    sigtimdifstdset(1,101-i)=sigtimdifcov;%��¼�õ��ķ�ֵ����ı�׼��
                else
                    sigtimdifstdset(1,101-i)=inf;
                end

                
             end
        end
       %step2 OVER ����ʱ�̼����ƽ��ֵ=sigtimdifmeanset������ʱ�̼���ı�׼��=sigtimdifstdset
        %======================================
       
        %======================================
        %step3:�ҵ����ʵ� ʱ�̼���ı�׼���ȷ����Ӧ�����з���ֵ����ֵ�������ֵ��������tjns1set�ķ�ֵʱ�̵���ֵ
        %����ͨ����ֵthrtim=(i/100)*avermaxa;���õ���һЩ���������ֵ�ķ�ֵʱ�̣�������Щ��ֵʱ�̼���ı�׼����Сʱ��
        %����ζ�ŷ�ֵʱ�̼���ķֲ��ȽϾ��ȣ�Ȼ�������ֵ����������ֵ��һ����Χ�ڱ仯��ʱ�򣬱�׼��䣬
        %Ҳ��˵�����������ֵ�ڣ��õ��ķ�ֵʱ���ǲ���ģ������������ֵ�����䣬˵����ֵthrtim=(i/100)*avermaxa;��Ӧ��ʱ������Ҫ����ֵ��
        sigtimdifstdset1=sigtimdifstdset;
        [stdmin,i_best]=min(sigtimdifstdset1(sigtimdifstdset1~=0));%�ҵ���׼�����Сֵ
        thrtim_best=(1.01-i_best/100)*avermaxa;

        sigpulset=find(stjn>=thrtim_best);%����pulsetime0��1x34����˵����34��������ֵ�ķ�ֵʱ��
                while 1
                sigtimdiff=diff(sigpulset);
                temp_record=find(sigtimdiff<=100);
                    if ~isempty(temp_record)%������ڼ����3�����ڵ�ʱ�̣�ȥ��С����һ��
                       for j=1:length(temp_record)

                           if stjn(sigpulset(temp_record(j)))>=stjn(sigpulset(temp_record(j)+1))
                               %˵��ǰ��һ��ʱ�̵ķ���ֵ��,ȥ�������Ǹ�
                               sigpulset(temp_record(j)+1)=-1;
                               sigpulset=sigpulset(sigpulset>0);
                               break;
                           else%˵������һ��ʱ�̵ķ���ֵ��,ȥ��ǰ���Ǹ�
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
        %�Ѿ��ҵ�����С��׼���Ӧ����ֵ���������ֵ��һ�αȽϿ�������ڣ���׼�����С�ģ��õ��ķ�ֵʱ���ǲ��䣬��Ӧ����Ҫ���ʱ��
        %��ֵthrtim_best����С��׼��stdmin����С��׼���Ӧ�ĸ���secstdmin
        %������������ֵ�ķ�ֵʱ��pulsetime0
        %================================================
        %STEP4 ȥ����׼��ܴ�����У�������׼��ϴ�����У���Ϊ��������һ�£�
        %�������ı�׼�stdmin��̫������Ҫȥ��һЩ���ż�಻���������**********��ע�������Ҫ���ڲ���
%          if stdmin>140 %ȥ����׼��̫���һ��IPT
%             continue;
%          end
        %��С��ֵ����ı�׼�Ŀ���ǣ��������ֵ�������һ���Ƚ�С�ķ�Χ���沨��
        %ÿ����k����ֵ���������������С����Լ���ֵ��
        %�ҵ�����k������У���׼���ر��ģ�Ӧ�þ����Ǹ���̫�Եķ�ֵʱ���ˣ��������̫�Եķ�ֵʱ��ȥ����ʣ�µķ�ֵʱ�̵ı�׼��ͻ�Ƚ�С
        if stdmin<=0.8%******************����źŲ��ã�����0.4
            pulsetimef=pulsetime0;
            union_set1=union(union_set1,Rset{m,1});
            pulsetimefdiff=diff(pulsetimef);%��¼�õ��������е�ʱ�̼������
            pulsetimefstd=std(pulsetimefdiff);%%��¼�õ��������е�ʱ�̼���ı�׼��
            %if pulsetimefstd<70 %*********���պϺ�Ҫ����Ǳ�׼�stdmin������̫��ģ���׼���Ŀ����ֹ�����
            pulsetimefset(1,n).pulse=pulsetimef;
            pulsetimefset(1,n).std=pulsetimefstd/mean(pulsetimefdiff);
            pulsetimefset(1,n).m=m;%�ҵ�����ԭʼ����tjns1set�Ķ�Ӧλ��
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
IDallset=(1:length(pulsetimefset));%��ʼ��Ϊȫ����ID
tic
if ~isfield(pulsetimefset,'pulse')
    return;
end
if length(IDallset)==1
    pulsefinal(1).pulse=pulsetimefset(1,1).pulse;
    return;
end
pulsetimefset1=pulsetimefset;
        %step1���࣬��Ϊ���������ص��ķ���ʱ�̴���0.5������һ��MU���ŵ�
    %  ÿһ���������� pulseclass��pulseclass �ĵ�1��Ԫ�ؾ���  pulsetimefset �е����ڵ�һ��mu���±꣺pulsetimefset ��pulseclass ��
    while 1  %ִ�����ɴι������
        iq=1;
        clear pulseclass(iw).ID;
        pulseclass(iw).ID=[];%���ÿһ���ID�����룩
        array=zeros(1,IDallset(end)-1);
        parfor iq=1:length(IDallset)-1   %ִ��ÿһ��Ĺ������
            %Ѱ��IDallset(i)�ķ������У�ʱ����+-50��Χ�ڣ���֮��ķ������еĹ�ͬʱ�̣������ҽ��ҵ���ʱ�̼������������
            %
                pu1 = zeros(nsamples,1);
                pu2 = zeros(nsamples,1);
                pu1(pulsetimefset(1,IDallset(1)).pulse)=1;% ��Ϊ�Ѿ��������id�ᱻdel������һֱ��1
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
                    %��MU��IPT���ƶȣ������״���0.5.�������а�����СMU��Ϣ���������̫�ͣ���Щ��������СMU��Ϣ��IPT�ͱ������������ٶȺܿ졣
                    %���Ǳ����ǵ�СMU֮��Ҳ���Ѹ���������Ի��������ɡ������Ҫ��ȡ�ⲿ��MU�ˣ��Ǿͽ����ٶȣ����ⲿ����ֵ���ø�һ�㡣
                    %��ǰ�������õ�һ�㣬��������IPT�������ˡ���Щ����IPT�����ں���Ҳ�ᱻ�������ظ�����С��2�Ļᱻcut��
                                %pulseclass(iw).ID=[pulseclass(iw).ID,IDallset(i+iq)];
                                array(iq)=IDallset(i+iq);%���ҵ���ʱ�̼������������  pulseclass(iw).ID
                                %pulsetimefset1(1,IDallset(i+iq)).pulse=pulsetimefset(1,IDallset(i+iq)).pulse+(ipmaxxc-50-1);%ipmaxxc-50-1

                end
           % end
            %pulseclass(iw).ID����ͬһmu�ľ���
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
        NewIDallset=setdiff(IDallset,pulseclass(iw).ID);%ȥ���Ѿ������ID,�����Լ��IDallset,���Ƿ�ȫ������
        iw=iw+1;
        IDallset=NewIDallset;
        if isempty(IDallset)||length(IDallset)<3%��ֹ������û��ͬ������Ĵ���
            break;
        end
    end
    toc
    %step1���� over pulseclass �ǹ���Ľ��������ÿ��ʱ�̶��������м�¼Ϊpulsetimefset.pulse1
    %====================================================
    %%
    %====================================================
    %step2  ��ʱ�̽϶���࣬ѡ��һ����׼��С��������Ϊ pulsefinal
    %%��ʱ�̽��ٵ��࣬��Ϊ����MU��Ӧ�����о����Ǳ�׼���С��ֻ��˵�����а�����ʱ�̱Ƚ�׼ȷ
    %�������ڸ����Ƚ��٣��ܿ������еķ�������©��,��������ǰ���Щ���ж����
    %��union���ж��������У����������ܵĽ���MU�ķ���ʱ�̶�������������ų���������������ࣨunion����finalpulse�� �غ϶ȴ���0.3�Ľ���ʱ�̣�������������ȥ��
    %����unionȥ�غ��MU�������� pulsefinal


    %��������Ԫ�ظ����Ƿ����10��
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
    % pulseclass(i).judge �����1����ô��˵�����еĸ���С��10

    %���С��5��union����һ���Ѿ��������
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


    %�������30��ѡ��һ����׼��С��������Ϊ pulsefinal����һ���Ѿ��������
    %ѡÿ�������׼����С/����һ����������С�ã�,������ĿҪ����һ����
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
        classstdset(k).std=[classstdset(k).std,pulsetimefset(1,pulseclass(k).ID).std];%�ռ�ÿһ��ı�׼���Ͽռ����൱�ڶ��ѭ�����Ԫ��
        stdsetunique(k).std=unique(classstdset(k).std);
        unifindlen=zeros(1,length(stdsetunique(k).std));%һ��Ҫ��ʼ������Ȼ�´����п��ܳ��ִ���
        %���ÿһ������ı�׼�����ظ��ģ���ѡ�����ظ��ı�׼�����׼ȷ�ʸߣ�������������ĴΣ��ظ���Ҫ�����Ĵ�Ϊ��
        if length(stdsetunique(k).std)<length(classstdset(k).std)
            for ie=1:length(stdsetunique(k).std)
                unifind=find(classstdset(k).std==stdsetunique(k).std(ie));
                unifindlen(ie)=length(unifind);
            end
            [~,unilenmaxip]=max(unifindlen);%�õ���׼�������Ǹ���׼�������
            classstdip(k).ip=find(classstdset(k).std==stdsetunique(k).std(unilenmaxip));
            selstd(1,k)=stdsetunique(k).std(unilenmaxip);
            selstd(2,k)= classstdip(k).ip(1);%�����ѡ������һ���ģ���Ϊ������ı���С��λ���ȣ������Ҫһ�㲻��������
        else  %û���ظ��ı�׼���ѡ����С��
            [selstd(1,k),selstd(2,k)]=min(classstdset(k).std);%�����С��׼���λ��*****��С�������⣬���������Կ�
            classstdip(k).ip=selstd(2,k);
        end
        pulsefinal(k).pulse=pulsetimefset(1,pulseclass(k).ID(selstd(2,k))).pulse;%�ҳ�ÿһ����С��׼������Ӧ�ķ��ţ���Ϊ����Ҫ�ֽ��MU
        pulsefinal(k).unionpulsetime=pulseclass(k).union;
    end

    

timespend=toc;
timespend=timespend/59

end

