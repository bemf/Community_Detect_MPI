addpath('../');
graphics_toolkit("gnuplot");
file_name='../../labData/param.txt';
[problem_str,nobj,alg,gen,pop,numprocs,run,paraTopology,paraType]=variation_generation_engine(file_name);

for problem_i=1:length(problem_str)
    procs=[];
    speedups_1=[];
    serial_time1=1.0;
    for proc_num=1:numprocs
        filename1 = sprintf('LOG_%s_%s(%d)_%d_%d_%dnp_%s_%s',alg{1},problem_str{problem_i},nobj,pop,gen,proc_num,paraType{3},paraTopology{1});

        procs=[procs proc_num];
        
        %proc_num

        [run_index,time]=textread(['../../labData/LOG/',filename1,'.dat'],'%d %f');
        mean_time=sum(time)/length(run_index)
        if proc_num==1
            serial_time1=mean_time;
        end
        speedup=serial_time1/mean_time;
        speedups_1=[speedups_1 speedup];



    %for proc end
    end

    %fprintf('testinstance:%s%d(%d)\n',problem_str,problem_index(problem_i),nobj);
    %fprintf('speedup of peCAEA:\n');
    %for i=1:length(speedups_1)
        %fprintf('\twhen the num_procs is %d,the speedup is %f\n',i,speedups_1(i));
    %end

    %fprintf('speedup of CAEA-IDEAL:\n');
    %for i=1:length(speedups_2)
        %fprintf('\twhen the num_procs is %d,the speedup is %f\n',i,speedups_2(i));
    %end
    %fprintf('speedup of CAEA-H:\n');
    %for i=1:length(speedups_3)
        %fprintf('\twhen the num_procs is %d,the speedup is %f\n',i,speedups_3(i));
    %end
    %fprintf('speedup of peCAEA-B:\n');
    %for i=1:length(speedups_4)
        %fprintf('\twhen the num_procs is %d,the speedup is %f\n',i,speedups_4(i));
    %end
    %fprintf('speedup of peCAEA-C:\n');
    %for i=1:length(speedups_5)
        %fprintf('\twhen the num_procs is %d,the speedup is %f\n',i,speedups_5(i));
    %end

    header=sprintf('speedup-%s',problem_str{problem_i})
    gcf=figure('Name',header,'NumberTitle','off','DefaultAxesFontSize',8);
    %set(gcf,'position',[200 100 1000 750]);
    handle1=plot(procs,speedups_1,'-bo');






    xl=xlabel('The number of processors');
    yl=ylabel('Speedup');
    set(xl,'FontSize',14,'FontWeight','bold');
    set(yl,'FontSize',14,'FontWeight','bold');
    ylim([0.00,7]);
    %grid on;
    %output=sprintf('speedup_%s%d.eps',problem_str,problem_index(problem_i));
 %   print ("-color","-deps", output);
%for problem end
end
figure;

handle1=plot(1,1,'-bo','lineWidth',lineSize,'markerSize',markerSize);
hold on;
handle2=plot(1,1,'-mh','lineWidth',lineSize,'markerSize',markerSize);
hold on;
handle3= plot(1,1,'-c*','lineWidth',lineSize,'markerSize',markerSize);
hold on;
handle4=plot(1,1,'-kx','lineWidth',lineSize,'markerSize',markerSize);
hold on;
handle5=plot(1,1,'-g^','lineWidth',lineSize,'markerSize',markerSize);
axis([10,11,10,11]);
le=legend([handle1 handle2 handle3 handle4,handle5],'peCAEA','peCAEA^p','peCAEA^h','peCAEA^B','peCAEA^C');
set(le,'Location','northwest');
set(le,'FontSize',32);
legend('boxoff');
axis off;
%set(gcf,'position',[500,300,100,400]);
%print("-color","-deps",'speedup_ZDT_legend.eps');
