function [sigpulset] = dele_adjacent_time(sigpulset,stjn,interval)
    if (nargin<3)
        interval=80;
    end
        
    sigpulset=sort(sigpulset);
    while 1
    sigtimdiff=diff(sigpulset);
    temp_record=find(sigtimdiff<=interval);
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
end

