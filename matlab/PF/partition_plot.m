function [le]=partition_plot(proc_Index,f1,f2,numprocs)


%point_style={'yo','m+','c*','gx','bd','cs','k^','rh'};
point_style={'m+','c*','gx','bd','cs','k^','rh','yo'};
legend_set={};

iterator_head=0;
iterator_tail=0;
for procs=1:numprocs
    %find next head and tail
    iterator_head=iterator_tail+1;
    iterator_tail=iterator_head;
    while iterator_tail<length(proc_Index)&&proc_Index(iterator_tail) == procs-1
        iterator_tail=iterator_tail+1;
    %while end
    end
    iterator_tail=iterator_tail-1;
    if procs<=8
        plot(f1(iterator_head:iterator_tail),f2(iterator_head:iterator_tail),point_style{procs},'MarkerSize',9);
    else 
        plot(f1(iterator_head:iterator_tail),f2(iterator_head:iterator_tail),point_style{8},'MarkerSize',9);
    end
    hold on;

    str=sprintf('partition %d',procs-1);
    legend_set=[legend_set str];
 %for end
end

le=legend(legend_set);
%funtion end
end
