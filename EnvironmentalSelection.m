function Population = EnvironmentalSelection(Population,N,Z,Zmin)
global PopCon
% The environmental selection of NSGA-III
%------------------------------- Copyright --------------------------------
%--------------------------------------------------------------------------
    if isempty(Zmin)
        Zmin = ones(1,size(Z,2));
    end

    %% Non-dominated sorting
    Population_objs = CalObj(Population);
    [FrontNo,MaxFNo] = NDSort(Population_objs,PopCon,N);
    Next = FrontNo < MaxFNo;
    
    %% Select the solutions in the last front
    Last   = find(FrontNo==MaxFNo);
    Choose = LastSelection(Population_objs(Next,:),Population_objs(Last,:),N-sum(Next)+3,Z,Zmin);
    Next(Last(Choose)) = true;
    % Population for next generation
    Population = Population(Next,:);
end

function Choose = LastSelection(PopObj1,PopObj2,K,Z,Zmin)
% Select part of the solutions in the last front
                               % 原来是减号
    PopObj = [PopObj1;PopObj2] + repmat(Zmin,size(PopObj1,1)+size(PopObj2,1),1);
    [N,M]  = size(PopObj);
    N1     = size(PopObj1,1);
    N2     = size(PopObj2,1);
    NZ     = size(Z,1);

    %% Normalization
    % Detect the extreme points
    Extreme = zeros(1,M);
    w       = zeros(M)+eye(M)+1e-6;
    for i = 1 : M
        [~,Extreme(i)] = min(max(PopObj./repmat(w(i,:),N,1),[],2));
    end
    % Calculate the intercepts of the hyperplane constructed by the extreme
    % points and the axes
    Hyperplane = (PopObj(Extreme,:) + 0.01.*eye(M) )\ones(M,1);  % 原来是\   %????????????????????????????
    a = 1./Hyperplane;
    if any(isnan(a))
        a = max(PopObj,[],1)';
    end
    % Normalization
    a = a-Zmin';
    PopObj = PopObj./repmat(a',N,1);
    
    %% Associate each solution with one reference point
    % Calculate the distance of each solution to each reference vector
    Cosine   = 1 - pdist2(PopObj,Z,'cosine');
    Distance = repmat(sqrt(sum(PopObj.^2,2)),1,NZ).*sqrt(1-Cosine.^2);
    Distance = pdist2(PopObj,Z,'euclidean');%上面注释的两个语句也可以哦
    % Associate each solution with its nearest reference point
    [d,pi] = min(Distance',[],1);

    %% Calculate the number of associated solutions except for the last front of each reference point
    rho = hist(pi(1:N1),1:NZ);%除了最后front的每一个参考点，计算关联解的个数
    
    %% Environmental selection
    Choose  = false(1,N2);
    Zchoose = true(1,NZ);
    % Select K solutions one by one
    while sum(Choose) < K-1
        % Select the least crowded reference point
        Temp = find(Zchoose);
        Jmin = find(rho(Temp)==min(rho(Temp)));
        if isempty(Jmin)
            break
        end
        j    = Temp(Jmin(randi(length(Jmin))));
        I    = find(Choose==0 & pi(N1+1:end)==j);
        % Then select one solution associated with this reference point
        if ~isempty(I)
            if rho(j) == 0
                [~,s] = min(d(N1+I));
            else
                s = randi(length(I));
            end
            Choose(I(s)) = true;
            rho(j) = rho(j) + 1;
        else
            Zchoose(j) = false;
        end
    end
end