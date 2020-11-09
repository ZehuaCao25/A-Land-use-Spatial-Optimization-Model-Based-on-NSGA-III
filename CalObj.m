function PopObj = CalObj(x)
global M D name upper lower N PopCon const_num c
    NN = size(x,1);  % NN = 496
    switch name
        case 'DTLZ1'
            PopObj      = 100*(D-M+1+sum((x(:,M:end)-0.5).^2-cos(20.*pi.*(x(:,M:end)-0.5)),2));
            PopObj = 0.5*repmat(1+PopObj,1,M).*fliplr(cumprod([ones(NN,1),x(:,1:M-1)],2)).*[ones(NN,1),1-x(:,M-1:-1:1)];
        case 'DTLZ2'
            PopObj      = sum((x(:,M:end)-0.5).^2,2);
            PopObj = repmat(1+PopObj,1,M).*fliplr(cumprod([ones(size(PopObj,1),1),cos(x(:,1:M-1)*pi/2)],2)).*[ones(size(PopObj,1),1),sin(x(:,M-1:-1:1)*pi/2)];
        case 'DTLZ3'
            PopObj      = 100*(D-M+1+sum((x(:,M:end)-0.5).^2-cos(20.*pi.*(x(:,M:end)-0.5)),2));
            PopObj = repmat(1+PopObj,1,M).*fliplr(cumprod([ones(NN,1),cos(x(:,1:M-1)*pi/2)],2)).*[ones(NN,1),sin(x(:,M-1:-1:1)*pi/2)];
        
            
            %%
        case '27'
            % v1 = [2.51 8.31 22.63 16.59 0 0];
            % v1 = [17.85 29.04 5.46 4.77 2899.68 219.18 579.96 7.14 0.00] 
            % v2 = [0.69  0.86  2.44 1.02 -0.98   -0.33   -0.26 3.97 0.10]
            %       耕地   园地 林地 草地  城镇    农村    风景 水域 未利用
            
            v1 = [4.77 1599.43  17.85   29.04  7.14   5.46  0 ];
            v2 = [1.02 -0.66    0.69    0.86   3.97   2.44  0.10];
            %      草地  聚落     耕地   园地    水域   森林   荒地
            
            %      草地 聚落     耕地    园地   水域   森林   荒地
            % 经济 4.77 1599.43  17.85   29.04  7.14   5.46  0
            % 环境 1.02 -0.66    0.69    0.86   3.97   2.44  0.10
            PopObj = 0.*ones(NN,M);
%%
            IntVar = [1:1:D];
            x(:, IntVar) = round(x(:,IntVar)); % Round Integer variables.
%%
c = 0.*ones(NN,const_num);
c = repmat([-3 -2 -1],NN,1);
PopCon = 0.*ones(NN,1);
    for k = 1:NN
        z = x(k,:);  % 第k行
        f1 = 0;     f2 = 0;     f3 = 0;
%%
        for i = 1:length(lower)          
            f1 = f1 + v1(z(i));
            f2 = f2 + v2(z(i));
        end
            PopObj(k,1) =  f1;
            PopObj(k,2) =  f2;
%%          
            DDD = sqrt(D);
        for i = 2+DDD:length(upper)-(2+DDD) % 聚集度
            if (z(i) == z(i+1) || z(i) == z(i-1) || z(i) == z(i-DDD) || z(i) == z(i+DDD) ||...
                z(i) == z(i+DDD-1) || z(i) == z(i+DDD+1) || z(i) == z(i-DDD-1) || z(i) == z(i-DDD+1))
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
    end
            
    end
end