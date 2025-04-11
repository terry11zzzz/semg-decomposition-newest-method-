% load('D:\dynamic_decompose\real_time_decomposition\Raw_S1_to_S6\S1\CaThKT_30_GM.mat')
% temp=reshape(SIG,13*5,1);
% bad_channel=reshape(discardChannelsVec,13*5,1);
% for i=2:13*5
%     sEMG(i,:)=temp{i};%(1:fs*period+1);
% end
% sEMG(bad_channel==1,:)=[];
% sEMG(1,:)=[];
% start_pt=50001-100;end_pt=80000;
%clc;clear;
% mission='dynamic';trial=2;
% sub=3;Task=1;str={'Yang','Nhi','Micheal','Chen','Bryan','Dorothy','Tien','Nick','Test1'};%动态
% start_point(1,1)=39901;end_point(1,1)=48000;start_point(1,2)=20901;end_point(1,2)=27000;start_point(1,3)=34901;end_point(1,3)=43000;
% start_point(2,1)=29901;end_point(2,1)=38000;start_point(2,2)=27901;end_point(2,2)=34000;
% start_point(3,1)=19901;end_point(3,1)=29000;start_point(3,2)=37901;end_point(3,2)=44000;start_point(3,3)=38901;end_point(3,3)=45000;
% start_point(4,1)=29901;end_point(4,1)=40000;start_point(4,2)=53901;end_point(4,2)=63000;start_point(4,3)=27901;end_point(4,3)=36000;
% start_point(5,1)=39901;end_point(5,1)=50000;start_point(5,2)=43901;end_point(5,2)=54000;start_point(5,3)=39901;end_point(5,3)=49000;
% start_point(6,1)=32901;end_point(6,1)=43000;start_point(6,2)=29901;end_point(6,2)=40000;start_point(6,3)=39901;end_point(6,3)=49000;
% start_point(7,1)=28901;end_point(7,1)=40000;start_point(7,2)=22901;end_point(7,2)=32000;start_point(7,3)=43901;end_point(7,3)=54000;
% start_point(8,1)=26901;end_point(8,1)=38000;start_point(8,2)=53901;end_point(8,2)=63000;start_point(8,3)=15901;end_point(8,3)=24000;
% start_point(9,1)=12901;end_point(9,1)=24000;start_point(9,2)=45901;end_point(9,2)=55000;start_point(9,3)=63901;end_point(9,3)=73000;                
% load_address1=strcat('D:\dynamic_decompose\real data\',str{sub},'\',str{sub},'_Dynamic_Task',num2str(Task),'_Trial',num2str(trial),'.mat');
% load_address2=strcat('D:\dynamic_decompose\real data\',str{sub},'\',str{sub},'_Dynamic_Task',num2str(Task),'_Trail',num2str(trial),'.mat');
% if exist(load_address1, 'file')
%     load_address = load_address1;
%     fprintf('文件存在，使用路径：%s\n', load_address);
% elseif exist(load_address2, 'file')
%     load_address = load_address2;
%     fprintf('第一个文件不存在，使用备用路径：%s\n', load_address);
% else
%     error('两个文件均不存在，程序无法继续。');
% end
% load(load_address)
% load('D:\dynamic_decompose\real data\bad_chan_dynamic.mat')
% sEMG(bad_chan{sub,Task},:)=[];
% if size(sEMG,1)>110
%     if mod(size(sEMG,1),2)==0
%         l=size(sEMG,1)/2;
%     else
%         l=(size(sEMG,1)-1)/2;
%     end
%     for i=1:l
%         signal(i,:)=sEMG(2*i,:);%静态 81000:93000 动态 66000 75000
%     end
% end
%%
sEMG=trials{1,3}';
% if mod(size(sEMG,2),2)==0
%     l=size(sEMG,1)/2;
% else
%     l=(size(sEMG,2)-1)/2;
% end
% for i=1:l
%     signal(:,i)=sEMG(:,2*i);
% end
para = [19, 10, 14, 40];
[finalres, union_set, tjns1set] = semgTechVisualization(sEMG, 3901, 24000, 10000, [], para,[]);
%%
L=16;
th=0.30;
record=[];
badMU_idx2=[];
ep=20100;
for w=1:length(finalres)-1 
    for  v=w+1:length(finalres)
         pu1 = zeros(ep,1);
        pu2 = zeros(ep,1);
        pu1(finalres(v).pulse)=1;
        for loop1=-L:1:L
            temp1=finalres(w).pulse+loop1;
            temp1(temp1>ep)=[];
            temp1(temp1<1)=[];
            pu2(temp1)=1;
        end
        xcpul=xcorr(pu1,pu2,50);
        [maxxc,ipmaxxc]=max(xcpul);
        tp1=maxxc/length(finalres(w).pulse);
        tp2=maxxc/length(finalres(v).pulse);
        if tp1>=th && tp2>=th && w~=v %%目的达到，及早跳出循环，节省时间参数选择决定结果准确性
            record=[record;w v tp1 tp2];
            badMU_idx2=[badMU_idx2,v];
        end
    end
end
badMU_idx1=unique(badMU_idx2)
finalres(badMU_idx1)=[];
%%
for i=1:length(finalres)
    finalres(i).mean_firing_rate=fs/finalres(i).mean_firing_gap;
    finalres(i).firing_time=finalres(i).pulse/fs;
end