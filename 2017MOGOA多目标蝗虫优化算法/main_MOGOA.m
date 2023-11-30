%_________________________________________________________________________________
%  Multi-objective Grasshopper Optimization Algorithm (MOGOA) source codes version 1.0
%
%  Developed in MATLAB R2016a
%
%  Author and programmer: Seyedali Mirjalili
%
%         e-Mail: ali.mirjalili@gmail.com
%                 seyedali.mirjalili@griffithuni.edu.au
%
%       Homepage: http://www.alimirjalili.com
%
%   Main paper:
%   S. Z. Mirjalili, S. Mirjalili, S. Saremi, H. Fatis, H. Aljarah, 
%   Grasshopper optimization algorithm for multi-objective optimization problems, 
%   Applied Intelligence, 2017, DOI: http://dx.doi.org/10.1007/s10489-017-1019-8
%____________________________________________________________________________________

% clc;
% clear;
% close all;

% % Change these details with respect to your problem%%%%%%%%%%%%%%
% ObjectiveFunction=@ZDT1;
% dim=5;
% lb=0;
% ub=1;
% obj_no=2;
% 
% if size(ub,2)==1
%     ub=ones(1,dim)*ub;
%     lb=ones(1,dim)*lb;
% end

clear all;
close all;
clc;

%%%% ����ʵ�������Χ������ʵ�������ʵ�麯����Χ
Num_Test=5;   %%%% ÿ��������������Num_Test��?
Num_Experiment=30;   %%%% �����Ǵ�F1-FNum_Functions
AlgorithmName='MOGOA'; %%% ���ƺ�����?


% gmax = 100;    %����������
% % FEE = 5000; %%%���Ŀ�꺯�����۴���
% n = 50;       %��Ⱥ��ģ
max_iter=200;  % Maximum Number of Iterations
N=200;    % Population Size (Number of Sub-Problems)
ArchiveMaxSize=100;  %%%Archive Size(number of rep)

m = 5;   %Ŀ��ά��

ALLFunction_AllTest=[];

% for ff=[1:21];
% for ff=[10:21];
% for ff=[1:4,6:9];
for ff=[7];
    clearvars -except Num_Test Num_Experiment AlgorithmName ALLFunction_AllTest ff max_iter N ArchiveMaxSize m AllTest_Results  problem_name 
    %%%%% �����ļ���?
    string_0ALL=['000\',AlgorithmName,'_5άĿ��800�ε���100��Ⱥʵ��20210923\'];
    dirname00=[string_0ALL,'\F',num2str(ff),'\'];
   display(['**********  ',AlgorithmName,'�㷨�Ż�F',num2str(ff),'�� ', 'M',num2str(m), ' άʵ��   **********']);
   for testi=1:Num_Test   %%%% ����ÿ��ʵ����Դ���
       dirname0=[dirname00,'test',num2str(testi),'_F',num2str(ff)];
       system(['mkdir ' dirname0]) %�������ļ���
       dirname1=[dirname0,'\F',num2str(ff),'_fig'];
       system(['mkdir ' dirname1]) %�����ļ���  �ȴ�����ʵ��ͼ��
       dirname2=[dirname0,'\F',num2str(ff),'_data'];
       system(['mkdir ' dirname2]) %�����ļ���  �ȴ�����ʵ��ͼ��
       for kk=1:30 %%%% ����ʵ�������ѭ��
           display(['**********  ',AlgorithmName,'�㷨�Ż�F',num2str(ff),'��  ��  ', num2str(kk), ' ��ʵ��   **********']);
            rand('state',sum(100*clock));
            problem_name=['F',num2str(ff)];
   
           [ ub,lb,dim ] = generate_boundary1( problem_name,m );%Upper and Lower Bound of Decision Variables  %%%���ɾ��߿ռ��б����Ͻ硢�½��ά��
tic;  % CPU time measure
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
flag=0;
if (rem(dim,2)~=0)
    dim = dim+1;
    ub = [ub, 1];
    lb = [lb, 0];
    flag=1;
end


% max_iter=100;
% N=200;
% ArchiveMaxSize=100;

Archive_X=zeros(100,dim);
Archive_F=ones(100,m)*inf;

Archive_member_no=0;

%Initialize the positions of artificial whales
GrassHopperPositions=initialization(N,dim,ub,lb);

TargetPosition=zeros(dim,1);
TargetFitness=inf*ones(1,m);

cMax=1;
cMin=0.00004;
%calculate the fitness of initial grasshoppers

for iter=1:max_iter
    for i=1:N
        
        Flag4ub=GrassHopperPositions(:,i)>ub';
        Flag4lb=GrassHopperPositions(:,i)<lb';
        GrassHopperPositions(:,i)=(GrassHopperPositions(:,i).*(~(Flag4ub+Flag4lb)))+ub'.*Flag4ub+lb'.*Flag4lb;
        
%         GrassHopperFitness(i,:)=ObjectiveFunction(GrassHopperPositions(:,i)');
        GrassHopperFitness(i,:)=test_function(GrassHopperPositions(:,i)',dim,m,problem_name);
        if dominates(GrassHopperFitness(i,:),TargetFitness)
            TargetFitness=GrassHopperFitness(i,:);
            TargetPosition=GrassHopperPositions(:,i);
        end
        
    end
    
    [Archive_X, Archive_F, Archive_member_no]=UpdateArchive(Archive_X, Archive_F, GrassHopperPositions, GrassHopperFitness, Archive_member_no);
    
    if Archive_member_no>ArchiveMaxSize
        Archive_mem_ranks=RankingProcess(Archive_F, ArchiveMaxSize, m);
        [Archive_X, Archive_F, Archive_mem_ranks, Archive_member_no]=HandleFullArchive(Archive_X, Archive_F, Archive_member_no, Archive_mem_ranks, ArchiveMaxSize);
    else
        Archive_mem_ranks=RankingProcess(Archive_F, ArchiveMaxSize, m);
    end
    
    Archive_mem_ranks=RankingProcess(Archive_F, ArchiveMaxSize, m);
    index=RouletteWheelSelection(1./Archive_mem_ranks);
    if index==-1
        index=1;
    end
    TargetFitness=Archive_F(index,:);
    TargetPosition=Archive_X(index,:)';
    
    c=cMax-iter*((cMax-cMin)/max_iter); % Eq. (3.8) in the paper
    
    for i=1:N
        
        temp= GrassHopperPositions;
        
        for k=1:2:dim
            S_i=zeros(2,1);
            for j=1:N
                if i~=j
                    Dist=distance(temp(k:k+1,j), temp(k:k+1,i));
                    r_ij_vec=(temp(k:k+1,j)-temp(k:k+1,i))/(Dist+eps);
                    xj_xi=2+rem(Dist,2);
                       
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% % Eq. (3.2) in the paper 
                    s_ij=((ub(k:k+1)' - lb(k:k+1)') .*c/2)*S_func(xj_xi).*r_ij_vec;
                    S_i=S_i+s_ij;
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                end
            end
            S_i_total(k:k+1, :) = S_i;
            
        end
        
        X_new=c*S_i_total'+(TargetPosition)'; % Eq. (3.7) in the paper
        GrassHopperPositions_temp(i,:)=X_new';
    end
    % GrassHopperPositions
    GrassHopperPositions=GrassHopperPositions_temp';
    
%     display(['At the iteration ', num2str(iter), ' there are ', num2str(Archive_member_no), ' non-dominated solutions in the archive']);
    HisPF{iter} = Archive_F;
end


if (flag==1)
    TargetPosition = TargetPosition(1:dim-1);
end

% figure
% 
% Draw_ZDT1();
% 
% hold on
% 
% plot(Archive_F(:,1),Archive_F(:,2),'ro','MarkerSize',8,'markerfacecolor','k');
% 
% legend('True PF','Obtained PF');
% title('MOGOA');
% 
% set(gcf, 'pos', [403   466   230   200])


            time=toc;
            PF = Archive_F;
            cg_curve=HisPF; %%% ��ʷĿ�꺯��ֵ
            Time(kk)=time;
            cc=strcat(dirname2,'\',AlgorithmName,'�Ż�����_',num2str(kk),'.mat');
            result.time=time;
           
            true_PF=TPF(m,ArchiveMaxSize, problem_name);

         %%% using the matlab codes for calculating metric values
         

            hv = HV(PF,true_PF);   %�����?
            gd = GD(PF, true_PF);                       %�������릻

            sp = Spacing(PF, true_PF);                  %�ռ�ֲ� 
            igd = IGD(PF, true_PF);            %�����������릻

            hvd(kk)=hv;
            gdd(kk)=gd;
            ssp(kk)=sp;
            igdd(kk)=igd;

             save(cc)
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       end
       mean_IGD = mean(igdd);
       display([AlgorithmName,'_Functions_F',num2str(ff),'����',num2str(kk),'��ʵ���IGDƽ��ֵ : ', num2str(mean_IGD)]);
       std_IGD=std(igdd);
       display([AlgorithmName,'_Functions_F',num2str(ff),'����',num2str(kk),'��ʵ���IGD��׼�� : ', num2str(std_IGD)]);
       max_IGD=max(igdd);
       display([AlgorithmName,'_Functions_F',num2str(ff),'����',num2str(kk),'��ʵ���IGD���ֵ : ', num2str(max_IGD)]);
       min_IGD=min(igdd);
       display([AlgorithmName,'_Functions_F',num2str(ff),'����',num2str(kk),'��ʵ���IGD��Сֵ : ', num2str(min_IGD)]);
       display('******************************** ');
       
       
       mean_GD = mean(gdd);
       display([AlgorithmName,'_Functions_F',num2str(ff),'����',num2str(kk),'��ʵ���GDƽ��ֵ : ', num2str(mean_GD)]);
       std_GD=std(gdd);
       display([AlgorithmName,'_Functions_F',num2str(ff),'����',num2str(kk),'��ʵ���GD��׼�� : ', num2str(std_GD)]);
       max_GD=max(gdd);
       display([AlgorithmName,'_Functions_F',num2str(ff),'����',num2str(kk),'��ʵ���GD���ֵ : ', num2str(max_GD)]);
       min_GD=min(gdd);
       display([AlgorithmName,'_Functions_F',num2str(ff),'����',num2str(kk),'��ʵ���GD��Сֵ : ', num2str(min_GD)]);
       display('******************************** ');
      
       
       mean_HV = mean(hvd);
       display([AlgorithmName,'_Functions_F',num2str(ff),'����',num2str(kk),'��ʵ���HVƽ��ֵ : ', num2str(mean_HV)]);
       std_HV=std(hvd);
       display([AlgorithmName,'_Functions_F',num2str(ff),'����',num2str(kk),'��ʵ���HV��׼�� : ', num2str(std_HV)]);
       max_HV=max(hvd);
       display([AlgorithmName,'_Functions_F',num2str(ff),'����',num2str(kk),'��ʵ���HV���ֵ : ', num2str(max_HV)]);
       min_HV=min(hvd);
       display([AlgorithmName,'_Functions_F',num2str(ff),'����',num2str(kk),'��ʵ���HV��Сֵ : ', num2str(min_HV)]);
       display('******************************** ');
       
       mean_SP = mean(ssp);
       display([AlgorithmName,'_Functions_F',num2str(ff),'����',num2str(kk),'��ʵ���SPƽ��ֵ : ', num2str(mean_SP)]);
       std_SP=std(ssp);
       display([AlgorithmName,'_Functions_F',num2str(ff),'����',num2str(kk),'��ʵ���SP��׼�� : ', num2str(std_SP)]);
       max_SP=max(ssp);
       display([AlgorithmName,'_Functions_F',num2str(ff),'����',num2str(kk),'��ʵ���SP���ֵ : ', num2str(max_SP)]);
       min_SP=min(ssp);
       display([AlgorithmName,'_Functions_F',num2str(ff),'����',num2str(kk),'��ʵ���SP��Сֵ : ', num2str(min_SP)]);
       display('******************************** ');
       
       mean_time=mean(Time);
       display([AlgorithmName,'_Functions_F',num2str(ff),'����',num2str(kk),'��ʵ�������ʱ��ƽ��ֵ : ', num2str(mean_time)]);
       std_time=std(Time);
       display([AlgorithmName,'_Functions_F',num2str(ff),'����',num2str(kk),'��ʵ�������ʱ���׼�� : ', num2str(std_time)]);
       display('******************************** ');
%        mean_X=mean(Best_X);
%         display([AlgorithmName,'_Functions_F',num2str(ff),'����',num2str(kk),'��ʵ������Ž�ƽ��ֵ : ', num2str(mean_X)]);
%         std_X=std(Best_X);
%         display([AlgorithmName,'_Functions_F',num2str(ff),'����',num2str(kk),'��ʵ������Ž��׼�� : ', num2str(std_X)]);
%         %%%%%%%%%%%%%%%%%%
        cd=strcat(dirname0,'\Result���ܽ��.mat');
        Result.IGDmean=mean_IGD;
        Result.IGDstd=std_IGD;
        Result.IGDmax=max_IGD;
        Result.IGDmin=min_IGD;
      
        
        Result.GDmean=mean_GD;
        Result.GDstd=std_GD;
        Result.GDmax=max_GD;
        Result.GDmin=min_GD;
        
        Result.HVmean=mean_HV;
        Result.HVstd=std_HV;
        Result.HVmax=max_HV;
        Result.HVmin=min_HV;
        
        Result.SPmean=mean_SP;
        Result.SPstd=std_SP;
        Result.SPmax=max_SP;
        Result.SPmin=min_SP;
        
        Result.Tmean=mean_time;
        Result.Tstd=std_time;
        %         Result.Xmean=mean_X;
        %         Result.Xstd=std_X;
        %         Result.Best_Y=Best_Y;
        %         Result.Best_X=Best_X;
        Result.Time=Time;
        %         Result.ResultVector=[mean_IGD,std_IGD,max_IGD,min_IGD,mean_GD,std_GD,max_GD,min_GD,mean_time,std_time];
        Result.ResultVector=[mean_IGD,std_IGD,max_IGD,min_IGD,mean_GD,std_GD,max_GD,min_GD,mean_HV,std_HV,max_HV,min_HV,mean_SP,std_SP,max_SP,min_SP,mean_time,std_time];
        %         Result.Best_History_Y=History_Y;
        save(cd,'Result')
        
        
        %         AllTest_Results(testi,:)=[mean_IGD,std_IGD,max_IGD,min_IGD,mean_time,std_time];
        % AllTest_Results(testi,:)=[mean_IGD,std_IGD,max_IGD,min_IGD,mean_GD,std_GD,max_GD,min_GD,mean_time,std_time];
        AllTest_Results(testi,:)=[mean_IGD,std_IGD,mean_GD,std_GD,mean_HV,std_HV,mean_SP,std_SP,mean_time,std_time];
    end
    cd=strcat(dirname00,'Result_AllTest.mat');
    save(cd,'AllTest_Results')
    ALLFunction_AllTest=[ALLFunction_AllTest;AllTest_Results];
end


cd=strcat(string_0ALL,'ALLFunction_AllTest.mat');
save(cd,'ALLFunction_AllTest')