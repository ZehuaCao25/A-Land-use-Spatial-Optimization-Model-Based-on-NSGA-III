clc,clear all
close all
tic
global N M D  PopCon name gen Cons const_num 
 %约束逻辑矩阵
N = 20;                        % 种群个数
M = 3;                          % 目标个数
name = '27';                 % 测试函数选择，目前只有：DTLZ1、DTLZ2、DTLZ3
gen = 300;                      %迭代次数
const_num = 3;

%% Generate the reference points and random population
EEE = xlsread('yulin.xlsx');
EEE = reshape(EEE,1,22500);
% for i = 1:15
%     Population(i,:) = EEE;
% end

[Z,N] = UniformPoint(N,M);        % 生成一致性参考解
[res,Population,PF] = funfun(); % 生成初始种群与目标值     [PopObj,PopDec,P] =  funfun()
% for i = 1:15
%     Population(i,:) = [fliplr(EEE),EEE,EEE,fliplr(EEE)];
% end
% hhh = randi([1,7],1,10000);
for i = 1:1:15
    Population(i,:) = EEE(1:2500);
end
Pop_objs = -CalObj(Population); % 计算适应度函数值，加了负号
Zmin  = min(Pop_objs(all(PopCon<=0,2),:),[],1); %求理想点，其实PopCon是处理有约束问题的，这里并没有用到

%% Optimization
HHH = [];
for i = 1:gen
    MatingPool = TournamentSelection(2,N,sum(max(0,PopCon),2));
    Offspring  = GA(Population(MatingPool,:));
    %% 加入整数 
    Offspring = round(Offspring);
    %%
    Offspring_objs = -CalObj(Offspring); % 函数值
    Zmin       = min([Zmin;Offspring_objs],[],1);%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Zmin
    Population = EnvironmentalSelection([Population;Offspring;],N,Z,Zmin); 
    Popobj = -CalObj(Population);
    %%
    [FrontNo,MaxFNo] = NDSort(Population,PopCon,round(N));
    A = find(FrontNo == 1);
    Pop_Front = Population(A',:);
    Obj_front = -CalObj(Pop_Front);
    
    HHH = [HHH;Obj_front];
    %%
    if(M<=3)
        %figure(1);
%         if(1<=i)
%             plot3(Popobj(:,1),Popobj(:,2),Popobj(:,3),'ro');
%             hold on;
%         end
        size(Popobj)
        if (1<=i)
            plot3(Obj_front(:,1),Obj_front(:,2),Obj_front(:,3),'ro');
            hold on;       
        end
%         if (200<=i)
%             plot3(Obj_front(:,1),Obj_front(:,2),Obj_front(:,3),'go');
%             hold on;       
%         end
        grid on;
        title(num2str(i));
        drawnow;  % 更新图窗并处理回调
    end
end
% if(M<=3)
%     hold on
%     plot3(PF(:,1),PF(:,2),PF(:,3),'g*')
%     grid on;
% else
%     for i=1:size(Popobj,1)
%         plot(Popobj(i,:))
%         grid on;
%         hold on
%     end
% end
%%IGD
score = IGD(Popobj,PF);
t_III = toc