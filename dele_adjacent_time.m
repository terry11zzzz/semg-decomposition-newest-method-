function [sigpulset] = dele_adjacent_time(sigpulset,stjn,interval)
    if (nargin<3)
        interval=80;
    end
        
    sigpulset=sort(sigpulset);
    while 1
    sigtimdiff=diff(sigpulset);
    temp_record=find(sigtimdiff<=interval);
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
end

