function [PopObj,x,P] =  funfun()

global M D lower upper encoding N PopCon cons cons_flag name xl xu IntVar const_num err c

    %% Initialization
    % D = M + 4;
    D = 2500;
    lower    = 1.*ones(1,D);
    upper    = 7.*ones(1,D);
    xl = lower(1);
    xu = upper(1);
    
    encoding = 'real';
    switch encoding
        case 'binary'
            x = randi([0,1],N,D);
        case 'permutation'
            [~,x] = sort(rand(N,D),2);
        otherwise
%             x = unifrnd(repmat(lower,N,1),repmat(upper,N,1));
%             x = round(x);
%             xl_temp=repmat(xl, pop_size,1);
%             xu_temp=repmat(xu, pop_size,1);
             IntVar = [1:1:D];
             x = lower+((upper - lower).*rand(N,D));
             x(:, IntVar) = round(x(:,IntVar)); % Round Integer variables.
    end
    cons = zeros(size(x,1),1);
    cons_flag = 1;
    PopCon = cons;
%                                     ff = zeros(pop_size, M);
%                                     err = zeros(pop_size, Const);
    %% Calculate objective values
    switch name
        case 'DTLZ1'
            PopObj      = 100*(D-M+1+sum((x(:,M:end)-0.5).^2-cos(20.*pi.*(x(:,M:end)-0.5)),2));
            % 378*1
            PopObj = 0.5*repmat(1+PopObj,1,M).*fliplr(cumprod([ones(N,1),...  % 目标值
                x(:,1:M-1)],2)).*[ones(N,1),1-x(:,M-1:-1:1)];
            P = UniformPoint(N,M)/2;  % 产生均匀分布的点
        case 'DTLZ2'
            PopObj      = sum((x(:,M:end)-0.5).^2,2);
            PopObj = repmat(1+PopObj,1,M).*fliplr(cumprod([ones(size(PopObj,1),1),cos(x(:,1:M-1)*pi/2)],2)).*[ones(size(PopObj,1),1),sin(x(:,M-1:-1:1)*pi/2)];
            P = UniformPoint(N,M);
            P = P./repmat(sqrt(sum(P.^2,2)),1,M);
        case 'DTLZ3'
            PopObj      = 100*(D-M+1+sum((x(:,M:end)-0.5).^2-cos(20.*pi.*(x(:,M:end)-0.5)),2));
            PopObj = repmat(1+PopObj,1,M).*fliplr(cumprod([ones(N,1),cos(x(:,1:M-1)*pi/2)],2)).*[ones(N,1),sin(x(:,M-1:-1:1)*pi/2)];
            P = UniformPoint(N,M);
            P = P./repmat(sqrt(sum(P.^2,2)),1,M);
        case '27'
    v1 = [4.77 1599.43  17.85   29.04  7.14   5.46  0];
    %v1 = repmat(v1,[N,1]);
    v2 = [1.02 -0.66    0.69    0.86   3.97   2.44  0.10];
    %v2 = repmat(v2,[N,1]);
    PopObj = 0.*ones(N,M);
    % P = UniformPoint(N,M);
    % P = P./repmat(sqrt(sum(P.^2,2)),1,M);
    
c = 0.*ones(N,const_num);
c = repmat([-3 -2 -1],N,1);
PopCon = 0.*ones(N,1);
    for k = 1:N
        z = x(k,:);     f1 = 0;     f2 = 0;     f3 = 0;
%%
        for i = 1:length(lower)          
            f1 = f1 + v2(z(i));
            f2 = f2 + v1(z(i));
        end
            PopObj(k,1) =  f1;
            PopObj(k,2) =  f2;
%%          
       len = sqrt(D);
        for i = (len + 2):length(D)-(len+2) % 聚集度
            if z(i) == z(i+1) || z(i) == z(i-1) || z(i) == z(i-len) || z(i) == z(i+len) ||...
                z(i) == z(i+len-1) || z(i) == z(i+len+5) || z(i) == z(i-len-1) || z(i) == z(i-len+1)
                f3 = f3 + 1;
            end
        end
        PopObj(k,3) = f3;
        
        %% 约束
        
    for i = 1:1:length(upper)
        if x(k,i) == 1
            c(k,1) = c(k,1) + 1;
        end
        if x(k,i) == 2
            c(k,2) = c(k,2) + 1;
        end
        if x(k,i) == 3
            c(k,3) = c(k,3) + 1;
        end
    end
    
    err = (c<0);  % 小于零是有问题的，
    err = sum(err,2);
    PopCon = err;
     %%   
    end

    P = UniformPoint(N,M)/2;        
    % f4 = x.*x;
    c(N,1:3) = [0];
%     for i = 1:1:length(upper)
%         if x(i) == 1
%             c(k,1) = c(k,1) + 1;
%         end
%         if x(i) == 2
%             c(k,2) = c(k,2) + 1;
%         end
%         if x(i) == 3
%             c(1,3) = c(1,3) + 1;
%         end
%     end
%     
%     err = (c>0).*c;
    case '28'
        PopObj = x.^2-exp(x);

        P = UniformPoint(N,M);  % 产生均匀分布的点   
    end
   %% Sample reference points on Pareto front
    %P = UniformPoint(N,obj.Global.M)/2;
end