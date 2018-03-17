addpath('../');
graphics_toolkit("gnuplot");
file_name='../../labData/param.txt';
[problem_str,nobj,alg,gen,pop,numprocs,run,paraTopology,paraType,dis_style]=variation_generation_engine(file_name);
lineSize=4;
markerSize=5;
%set (0, "defaultaxesfontname", "Helvetica");
for problem_i= 1:length(problem_str)
    header=sprintf('%s',problem_str{problem_i})
    figure('Name',header,'NumberTitle','off','DefaultAxesFontSize',28);
    %set(gca,'FontSize',20);
%</modifiable
    %filename1 = sprintf('Q_%s_%s%d(%d)_%d_%d_1np_%s_%s',alg{1},problem_str{problem_i},nobj,pop,gen,paraType{1},paraTopology{1});
    filename1 = sprintf('Q_%s_%s(%d)_%d_%d_%dnp_%s_%s',alg{1},problem_str{problem_i},nobj,pop,gen,numprocs,paraType{3},paraTopology{1})
 %   filename3 = sprintf('Q_%s_%s%d(%d)_%d_%d_%dnp_%s_%s',alg{1},problem_str,problem_index(problem_i),nobj,pop,gen,numprocs,paraType{2},paraTopology{1});
    %filename4 = sprintf('Q_%s_%s%d(%d)_%d_%d_%dnp_%s_%s',alg{1},problem_str,problem_index(problem_i),nobj,pop,gen,numprocs,paraType{3},paraTopology{1});
    %filename5 = sprintf('Q_%s_%s%d(%d)_%d_%d_%dnp_%s_%s',alg{1},problem_str,problem_index(problem_i),nobj,pop,gen,numprocs,paraType{1},paraTopology{2});
    %filename6 = sprintf('Q_%s_%s%d(%d)_%d_%d_%dnp_%s_%s',alg{1},problem_str,problem_index(problem_i),nobj,pop,gen,numprocs,paraType{1},paraTopology{3})
%modifiable/>
   [Index,Q] = textread(['../../labData/Q/',filename1,'.dat'],'%d %f');
    min_Q = 1e3;
    max_Q = -1;
    record_one = [];
    record_min = [];
    record_max = [];
    record_whole = [];
    Index_gen = [];
    
    ItemNum =length(Index)/run
    
    
    for j= 1:ItemNum
       Index_gen(j) = Index(j);
       record_whole(j) = 0;
    end
    
    for i = 1:run
        for j = 1:ItemNum
            record_one(j) = Q((i-1)*ItemNum + j);
            record_whole(j) = record_whole(j) + record_one(j);
        end
        if record_one(ItemNum) > max_Q 
            max_Q = record_one(ItemNum);
            for j = 1:ItemNum
                record_max(j) = record_one(j);
            end
        end
        if record_one(ItemNum) < min_Q 
            min_Q = record_one(ItemNum);
            for j = 1:ItemNum
                record_min(j) = record_one(j);
            end
        end
    end 
    
    avg_whole  = [];
    for j = 1:ItemNum
        avg_whole(j) = record_whole(j)/run;
    end
    
    for j = 1:ItemNum
        record_whole(j) = record_whole(j) - record_min(j);
        record_whole(j) = record_whole(j) - record_max(j);
    end
    
    avg_ex = [];
    for j = 1 :ItemNum
       avg_ex(j) = record_whole(j)/(run-2);
    end
    
    fprintf('== %s =====\n',filename1);
    fprintf('min_Q :\t%.6f\n',min_Q);
    fprintf('avg_whole :\t%.6f\n',avg_whole(ItemNum));
    fprintf('max_Q :\t%.6f\n',max_Q);
    fprintf('avg_ex :\t%.6f\n',avg_ex(ItemNum));
    handle1 = plot(Index_gen,avg_ex,dis_style{2},'MarkerSize',markerSize,'lineWidth',3);
    hold on;
    
    %%{
    %fprintf('*****avg_whole******\n')
    %for i = 1:ItemNum
       %fprintf('%.8f\n',avg_whole(i)); 
    %end
    %%}
    
    %[Index,Q] = textread(['../../labData/Q/',filename2,'.dat'],'%d %f');
    %min_Q = 1e3;
    %max_Q = -1;
    
    %%ItemNum =length(Index)/run
    %record_one = [];
    %record_min = [];
    %record_max = [];
    %record_whole = [];
    %for j= 1:ItemNum
       %record_whole(j) = 0;
    %end
    
    %for i = 1:run
        %for j = 1:ItemNum
            %record_one(j) = Q((i-1)*ItemNum + j);
            %record_whole(j) = record_whole(j) + record_one(j);
        %end
        %if record_one(ItemNum) > max_Q 
            %max_Q = record_one(ItemNum);
            %for j = 1:ItemNum
                %record_max(j) = record_one(j);
            %end
        %end
        %if record_one(ItemNum) < min_Q 
            %min_Q = record_one(ItemNum);
            %for j = 1:ItemNum
                %record_min(j) = record_one(j);
            %end
        %end
    %end

    %avg_whole  = [];
    %for j = 1:ItemNum
        %avg_whole(j) = record_whole(j)/run;
    %end

    %for j = 1:ItemNum
        %record_whole(j) = record_whole(j) - record_min(j);
        %record_whole(j) = record_whole(j) - record_max(j);
    %end

    %avg_ex = [];
    %for j = 1 :ItemNum
       %avg_ex(j) = record_whole(j)/(run-2);
    %end

    %fprintf('== %s =====\n',filename2);
    %fprintf('min_Q :\t%.6f\n',min_Q);
    %fprintf('avg_whole :\t%.6f\n',avg_whole(ItemNum));
    %fprintf('max_Q :\t%.6f\n',max_Q);
    %fprintf('avg_ex :\t%.6f\n',avg_ex(ItemNum));
    %handle2 = plot(Index_gen,avg_ex,dis_style{2},'lineWidth',lineSize);
    %hold on;





    %for i = 1:run
        %for j = 1:ItemNum
            %record_one(j) = Q((i-1)*ItemNum + j);
            %record_whole(j) = record_whole(j) + record_one(j);
        %end
        %if record_one(ItemNum) > max_Q 
            %max_Q = record_one(ItemNum);
            %for j = 1:ItemNum
                %record_max(j) = record_one(j);
            %end
        %end
        %if record_one(ItemNum) < min_Q 
            %min_Q = record_one(ItemNum);
            %for j = 1:ItemNum
                %record_min(j) = record_one(j);
            %end
        %end
    %end 
    
    %avg_whole  = [];
    %for j = 1:ItemNum
        %avg_whole(j) = record_whole(j)/run;
    %end
    
    %for j = 1:ItemNum
        %record_whole(j) = record_whole(j) - record_min(j);
        %record_whole(j) = record_whole(j) - record_max(j);
    %end
    
    %avg_ex = [];
    %for j = 1 :ItemNum
       %avg_ex(j) = record_whole(j)/(run-2);
    %end
    

    
    
    %[Index,Q] = textread(['../../labData/Q/',filename4,'.dat'],'%d %f');
    %min_Q = 1e3;
    %max_Q = -1;
    
    %record_one = [];
    %record_min = [];
    %record_max = [];
    %record_whole = [];
    %for j= 1:ItemNum
       %record_whole(j) = 0;
    %end
    
    %for i = 1:run
        %for j = 1:ItemNum
            %record_one(j) = Q((i-1)*ItemNum + j);
            %record_whole(j) = record_whole(j) + record_one(j);
        %end
        %if record_one(ItemNum) > max_Q 
            %max_Q = record_one(ItemNum);
            %for j = 1:ItemNum
                %record_max(j) = record_one(j);
            %end
        %end
        %if record_one(ItemNum) < min_Q 
            %min_Q = record_one(ItemNum);
            %for j = 1:ItemNum
                %record_min(j) = record_one(j);
            %end
        %end
    %end 
    
    %avg_whole  = [];
    %for j = 1:ItemNum
        %avg_whole(j) = record_whole(j)/run;
    %end
    
    %for j = 1:ItemNum
        %record_whole(j) = record_whole(j) - record_min(j);
        %record_whole(j) = record_whole(j) - record_max(j);
    %end
    
    %avg_ex = [];
    %for j = 1 :ItemNum
       %avg_ex(j) = record_whole(j)/(run-2);
    %end
    
    %fprintf('== %s =====\n',filename4);
    %fprintf('min_Q :\t%.6f\n',min_Q);
    %fprintf('avg_whole :\t%.6f\n',avg_whole(ItemNum));
    %fprintf('max_Q :\t%.6f\n',max_Q);
    %fprintf('avg_ex :\t%.6f\n',avg_ex(ItemNum));
    %handle4 = semilogy(Index_gen,avg_ex,dis_style{4},'MarkerSize',markerSize);
    %hold on;
    
    %[Index,Q] = textread(['../../labData/Q/',filename5,'.dat'],'%d %f');
    %min_Q = 1e3;
    %max_Q = -1;
    
    %record_one = [];
    %record_min = [];
    %record_max = [];
    %record_whole = [];
    %for j= 1:ItemNum
       %record_whole(j) = 0;
    %end
    
    %for i = 1:run
        %for j = 1:ItemNum
            %record_one(j) = Q((i-1)*ItemNum + j);
            %record_whole(j) = record_whole(j) + record_one(j);
        %end
        %if record_one(ItemNum) > max_Q 
            %max_Q = record_one(ItemNum);
            %for j = 1:ItemNum
                %record_max(j) = record_one(j);
            %end
        %end
        %if record_one(ItemNum) < min_Q 
            %min_Q = record_one(ItemNum);
            %for j = 1:ItemNum
                %record_min(j) = record_one(j);
            %end
        %end
    %end 
    
    %avg_whole  = [];
    %for j = 1:ItemNum
        %avg_whole(j) = record_whole(j)/run;
    %end
    
    %for j = 1:ItemNum
        %record_whole(j) = record_whole(j) - record_min(j);
        %record_whole(j) = record_whole(j) - record_max(j);
    %end
    
    %avg_ex = [];
    %for j = 1 :ItemNum
       %avg_ex(j) = record_whole(j)/(run-2);
    %end
    
    %fprintf('== %s =====\n',filename5);
    %fprintf('min_Q :\t%.6f\n',min_Q);
    %fprintf('avg_whole :\t%.6f\n',avg_whole(ItemNum));
    %fprintf('max_Q :\t%.6f\n',max_Q);
    %fprintf('avg_ex :\t%.6f\n',avg_ex(ItemNum));
    %handle5 = semilogy(Index_gen,avg_ex,dis_style{5},'lineWidth',lineSize);
    %hold on;
    
    %[Index,Q] = textread(['../../labData/Q/',filename6,'.dat'],'%d %f');
    %min_Q = 1e3;
    %max_Q = -1;
    
    %record_one = [];
    %record_min = [];
    %record_max = [];
    %record_whole = [];
    %for j= 1:ItemNum
       %record_whole(j) = 0;
    %end
    
    %for i = 1:run
        %for j = 1:ItemNum
            %record_one(j) = Q((i-1)*ItemNum + j);
            %record_whole(j) = record_whole(j) + record_one(j);
        %end
        %if record_one(ItemNum) > max_Q 
            %max_Q = record_one(ItemNum);
            %for j = 1:ItemNum
                %record_max(j) = record_one(j);
            %end
        %end
        %if record_one(ItemNum) < min_Q 
            %min_Q = record_one(ItemNum);
            %for j = 1:ItemNum
                %record_min(j) = record_one(j);
            %end
        %end
    %end

    %avg_whole  = [];
    %for j = 1:ItemNum
        %avg_whole(j) = record_whole(j)/run;
    %end

    %for j = 1:ItemNum
        %record_whole(j) = record_whole(j) - record_min(j);
        %record_whole(j) = record_whole(j) - record_max(j);
    %end

    %avg_ex = [];
    %for j = 1 :ItemNum
       %avg_ex(j) = record_whole(j)/(run-2);
    %end

    %fprintf('== %s =====\n',filename6);
    %fprintf('min_Q :\t%.6f\n',min_Q);
    %fprintf('avg_whole :\t%.6f\n',avg_whole(ItemNum));
    %fprintf('max_Q :\t%.6f\n',max_Q);
    %fprintf('avg_ex :\t%.6f\n',avg_ex(ItemNum));
    %handle6 = semilogy(Index_gen,avg_ex,dis_style{6},'MarkerSize',markerSize);
    %handle6 = semilogy(Index_gen,avg_ex,dis_style{6},'MarkerSize',markerSize,'markerfacecolor','g');
    %hold on;
    %[Index,Q] = textread(['../../labData/Q/',filename7,'.dat'],'%d %f');
    %min_Q = 1e3;
    %max_Q = -1;

    %record_one = [];
    %record_min = [];
    %record_max = [];
    %record_whole = [];
    %for j= 1:ItemNum
       %record_whole(j) = 0;
    %end

    %for i = 1:run
        %for j = 1:ItemNum
            %record_one(j) = Q((i-1)*ItemNum + j);
            %record_whole(j) = record_whole(j) + record_one(j);
        %end
        %if record_one(ItemNum) > max_Q 
            %max_Q = record_one(ItemNum);
            %for j = 1:ItemNum
                %record_max(j) = record_one(j);
            %end
        %end
        %if record_one(ItemNum) < min_Q 
            %min_Q = record_one(ItemNum);
            %for j = 1:ItemNum
                %record_min(j) = record_one(j);
            %end
        %end
    %end
    
    %avg_whole  = [];
    %for j = 1:ItemNum
        %avg_whole(j) = record_whole(j)/run;
    %end
    
    %for j = 1:ItemNum
        %record_whole(j) = record_whole(j) - record_min(j);
        %record_whole(j) = record_whole(j) - record_max(j);
    %end
    
    %avg_ex = [];
    %for j = 1 :ItemNum
       %avg_ex(j) = record_whole(j)/(run-2);
    %end
    
    %fprintf('== %s =====\n',filename7);
    %fprintf('min_Q :\t%.6f\n',min_Q);
    %fprintf('avg_whole :\t%.6f\n',avg_whole(ItemNum));
    %fprintf('max_Q :\t%.6f\n',max_Q);
    %fprintf('avg_ex :\t%.6f\n',avg_ex(ItemNum));
    %handle7 = semilogy(Index_gen,avg_ex,'y');hold on;
%</modifiable
%hold on;
    yl=ylabel('Q');
    xl=xlabel('No. of generations');
    set(xl,'FontSize',12,'FontWeight','bold');
    set(yl,'FontSize',12,'FontWeight','bold');
    %//karate///
    %ylim([0.3,0.45]);
    %set(gca,'YTick',[0.3:0.01:0.45]);
    %text(550,0.43,"0.4198");
    %dolphin
    %ylim([0.43,0.54]);
    %set(gca,'YTick',[0.43:0.01:0.54]);
    %text(550,0.53,"0.5285");
    %polbooks
    %ylim([0.22,0.54]);
    %set(gca,'YTick',[0.22:0.02:0.54]);
    %text(550,0.51,"0.5272");
    %adjnoun
    %ylim([0.16,0.32]);
    %set(gca,'YTick',[0.16:0.01:0.32]);
    %text(550,0.31,"0.3014");
    %football
    %ylim([0.18,0.62]);
    %set(gca,'YTick',[0.18:0.03:0.62]);
    %text(550,0.58,"0.6044");
    %JAZZ
    %ylim([0.05,0.46]);
    %set(gca,'YTick',[0.05:0.03:0.46]);
    %text(550,0.45,"0.4420");
    %netscience
    %ylim([0.85,0.97]);
    %set(gca,'YTick',[0.85:0.03:0.97]);
    %text(550,0.95,"0.9566");
    %email
    ylim([0.27,0.57]);
    set(gca,'YTick',[0.27:0.02:0.57]);
    text(1750,0.55,"0.5606");
    %annotation('arrow',[450 0.42],[500 0.43]);
    %annotation('line',550,0.43);
    %le=legend('CAEA','peCAEA','peCAEA^u','peCAEA^h','peCAEA^B','peCAEA^C');
    %set(le,'Location','northoutside');
    %set(le,'Orientation','horizontal');
    %print -color -dpdf  "-S560,420" 1.pdf
    %output=sprintf('Q_%s%d.eps',problem_str,problem_index(problem_i));
  %  print ("-color","-deps", output);
%modifiable/>
end
%figure;
%semilogy(Index_gen,avg_ex,dis_style{1},'MarkerSize',markerSize);hold on;
%semilogy(Index_gen,avg_ex,dis_style{2},'MarkerSize',markerSize,'lineWidth',lineSize);hold on;
%%semilogy(Index_gen,avg_ex,dis_style{3},'MarkerSize',markerSize,'lineWidth',lineSize);hold on;
%semilogy(Index_gen,avg_ex,dis_style{4},'MarkerSize',markerSize);hold on;
%semilogy(Index_gen,avg_ex,dis_style{5},'MarkerSize',markerSize,'lineWidth',lineSize);hold on;
%semilogy(Index_gen,avg_ex,dis_style{6},'MarkerSize',markerSize);hold on;
%axis([10,11,10,11]);
%le=legend('CAEA','peCAEA','peCAEA^p','peCAEA^B','peCAEA^C');
%set(le,'Location','west');
%set(le,'FontSize',12);
%legend('boxoff');
%axis off;
%set(gcf,'position',[500,300,100,400]);
%print("-color","-deps",'Q_ZDT_legend.eps');
