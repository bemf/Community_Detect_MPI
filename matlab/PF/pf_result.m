%</modifiable
addpath('../');
graphics_toolkit("gnuplot");
file_name='../../labData/param.txt';
[problem_str,nobj,alg,gen,pop,numprocs,run,paraTopology,paraType]=variation_generation_engine(file_name);
alg={'CAEA','MOEAD_TCH'};
%modifiable/>

for problem_i=1:length(problem_str)
    filename1 = sprintf('PF_%s_%s(%d)_%d_%dR%d_%dnp_%s_%s',alg{1},problem_str{problem_i},nobj,pop,gen,run,numprocs,paraType{3},paraTopology{1});
    %filename2 = sprintf('PF_%s_%s(%d)_%d_%dR%d_%dnp_%s_%s',alg{2},problem_str{problem_i},nobj,pop,gen,run_index,numprocs,paraType{1},paraTopology{2});
    
    
    header=sprintf('pf-%s',problem_str{problem_i})
    figure('Name',header,'NumberTitle','off','DefaultAxesFontSize',8);
    
    
    [proc_Index,f1,f2]=textread(['../../labData/POF/',filename1,'.dat'],'%d %f %f');
    le=partition_plot(proc_Index,f1,f2,numprocs);
    xl=xlabel('InterQ');
    yl=ylabel('IntraQ');
    %grid on;
    set(xl,'FontSize',14,'FontWeight','bold');
    set(yl,'FontSize',14,'FontWeight','bold');
    set(le,'FontSize',14);
%print  -deps 'peCAEA-ZDT2-front.eps';
end
