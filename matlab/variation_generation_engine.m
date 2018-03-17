function [problem_str,nobj,alg,gen,pop,numprocs,run,paraTopology,paraType,dis_style]=variation_generation_engine(file_name);

%</modifiable
    %note!!
    testInstance=fopen(file_name);
    for i=1:8 %note
        line{i}=fgets(testInstance);
    end
    run =str2num(line{5}); %note
    line{8}=deblank(line{8}); %去掉首尾多余空格
    splits=regexp(line{8},'\s+','split');
    nobj=str2num(splits{2});
    pop=str2num(splits{4});
    gen=str2num(splits{5});

    %problem_str='ZDT';
    %problem_index=[1,2,3,4,6];
    %problem_str={'karate','dolphin','polbooks','adjnoun','football','jazz','netscience','email','power'};
    problem_str={'email'};
    %problem_index=[1];
    %problem_str='DTLZ';
    %problem_index=[1,2,3,4];
    %problem_str='UF';
    %problem_index=[1,2,3,4,5,6,7];
    alg={'CAEA','MOEAD_TCH'};
    numprocs=1;
    %paraTopology={'FULL','BI_RING'};
    %paraTopology={'DYNAMIC','BI_RING','FULL'};
    paraTopology={'BI_RING','FULL'};
    paraType={'INDIV','IDEAL','IDEAL_INDIV'};
    dis_style={'--ok','-b','-m','--*r','-y','--hg','-ob','-.m','-^b'};
%modifiable/>
end

