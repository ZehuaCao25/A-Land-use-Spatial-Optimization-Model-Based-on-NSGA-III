function Offspring = GA(Parent,Parameter)
global encoding lower upper IntVar xl xu etac
%   GA - Genetic operators for real, binary, and permutation based encodings.
%   
%   Off = GA(P) returns the offsprings generated by genetic operators,
%   where P1 is a set of parents. If P is an array of INDIVIDUAL objects,
%   then Off is also an array of INDIVIDUAL objects; while if P is a matrix
%   of decision variables, then Off is also a matrix of decision variables,
%   i.e., the offsprings will not be evaluated. P is split into two subsets
%   P1 and P2 with the same size, and each object/row of P1 and P2 is used
%   to generate TWO offsprings. Different operators are used for real,
%   binary, and permutation based encodings, respectively.
%
%   Off = GA(P,{proC,disC,proM,disM}) specifies the parameters of
%   operators, where proC is the probabilities of doing crossover, disC is
%   the distribution index of simulated binary crossover, proM is the
%   expectation of number of bits doing mutation, and disM is the
%   distribution index of polynomial mutation.
%
%   Example:
%       Off = GA(Parent)
%       Off = GA(Parent.decs,{1,20,1,20})
%
%   See also GAhalf

%------------------------------- Reference --------------------------------

%--------------------------------------------------------------------------

    %% Parameter setting
    if nargin > 1
        [proC,disC,proM,disM] = deal(Parameter{:});
    else  % ��������
        [proC,disC,proM,disM] = deal(0.4, 20,0.8, 20);      % deal(1,20,1,20);
              etac = disC;
    end
    if isa(Parent(1),'INDIVIDUAL')
        calObj = true;
        Parent = Parent.decs;
    else
        calObj = false;
    end
    Parent1 = Parent(1:floor(end/2),:);
    Parent2 = Parent(floor(end/2)+1:floor(end/2)*2,:);
    [N,D]   = size(Parent1);
 
    switch encoding
        case 'binary'
            %% Genetic operators for binary encoding
            % One point crossover
            k = repmat(1:D,N,1) > repmat(randi(D,N,1),1,D);
            k(repmat(rand(N,1)>proC,1,D)) = false;
            Offspring1    = Parent1;
            Offspring2    = Parent2;
            Offspring1(k) = Parent2(k);
            Offspring2(k) = Parent1(k);
            Offspring     = [Offspring1;Offspring2];
            % Bitwise mutation
            Site = rand(2*N,D) < proM/D;
            Offspring(Site) = ~Offspring(Site);
        case 'permutation'
            %% Genetic operators for permutation based encoding
            % Order crossover
            Offspring = [Parent1;Parent2];
            k = randi(D,1,2*N);
            for i = 1 : N
                Offspring(i,k(i)+1:end)   = setdiff(Parent2(i,:),Parent1(i,1:k(i)),'stable');
                Offspring(i+N,k(i)+1:end) = setdiff(Parent1(i,:),Parent2(i,1:k(i)),'stable');
            end
            % Slight mutation
            k = randi(D,1,2*N);
            s = randi(D,1,2*N);
            for i = 1 : 2*N
                if s(i) < k(i)
                    Offspring(i,:) = Offspring(i,[1:s(i)-1,k(i),s(i):k(i)-1,k(i)+1:end]);
                elseif s(i) > k(i)
                    Offspring(i,:) = Offspring(i,[1:k(i)-1,k(i)+1:s(i)-1,k(i),s(i):end]);
                end
            end
        otherwise
            %% Genetic operators for real encoding
            % Simulated binary crossover
            beta = zeros(N,D);
            mu   = rand(N,D);
            beta(mu<=0.5) = (2*mu(mu<=0.5)).^(1/(disC+1));
            beta(mu>0.5)  = (2-2*mu(mu>0.5)).^(-1/(disC+1));
            beta = beta.*(-1).^randi([0,1],N,D);
            beta(rand(N,D)<0.5) = 1;
            beta(repmat(rand(N,1)>proC,1,D)) = 1;
            Offspring = [(Parent1+Parent2)/2+beta.*(Parent1-Parent2)/2
                         (Parent1+Parent2)/2-beta.*(Parent1-Parent2)/2];
            Offspring = round(Offspring); 
%%
%         [NNN] = size(Parent,1);
%         xl1=lower';
%         xu1=upper';
%         rc=randperm(NNN); % ��1��N����������û�
%         for i=1:(NNN/2)
%         parent1=Parent((rc(2*i-1)),:); % ����
%         parent2=Parent((rc(2*i)),:);   % ˫��
%         if (isequal(parent1,parent2))==1 & rand(1)>0.5
%             child1=parent1;
%             child2=parent2;
%         else 
%             for j = 1: D  
%                 if parent1(j)<parent2(j)
%                 beta(j)= 1 + (2.2/( 1+ parent2(j)-parent1(j)))  *  (   min(  (parent1(j)-xl1(j)),(xu1(j)-parent2(j))  )   );
%         else       %�� ������1     % ���������д�λ�ö�����1
%                 beta(j)= 1 + (2.2/( 1+ parent1(j)-parent2(j)))*(min((parent2(j)-xl1(j)),(xu1(j)-parent1(j))));
%                 end   
%            end
%         u=rand(1,D);
%         alpha=2-beta.^-(etac+1);        % dtac: distribution index for crossover
%         betaq=(u<=(1./alpha)).*(u.*alpha).^(1/(etac+1))+(u>(1./alpha)).*(1./(2 - u.*alpha)).^(1/(etac+1));
%         child1=0.5*(((1 + betaq).*parent1) + (1 - betaq).*parent2);
%         child2=0.5*(((1 - betaq).*parent1) + (1 + betaq).*parent2);
%         end
%         child_offspring((rc(2*i-1)),:)=poly_mutation(child1);           % polynomial mutation
%         child_offspring((rc(2*i)),:)=poly_mutation(child2);             % polynomial mutation
%         child_offspring(:, IntVar) = round(child_offspring(:,IntVar));
%         Offspring = child_offspring;
% 

%%
            % Polynomial mutation
            Lower = repmat(lower,2*N,1);
            Upper = repmat(upper,2*N,1);
            Site  = rand(2*N,D) < proM/D;
            mu    = rand(2*N,D);
            temp  = Site & mu<=0.5;
            Offspring       = min(max(Offspring,Lower),Upper);
            Offspring(temp) = Offspring(temp)+(Upper(temp)-Lower(temp)).*((2.*mu(temp)+(1-2.*mu(temp)).*...
                              (1-(Offspring(temp)-Lower(temp))./(Upper(temp)-Lower(temp))).^(disM+1)).^(1/(disM+1))-1);
            temp = Site & mu>0.5; 
            Offspring(temp) = Offspring(temp)+(Upper(temp)-Lower(temp)).*(1-(2.*(1-mu(temp))+2.*(mu(temp)-0.5).*...
                              (1-(Upper(temp)-Offspring(temp))./(Upper(temp)-Lower(temp))).^(disM+1)).^(1/(disM+1)));
            Offspring = round(Offspring);             
                          %%
%               pm = 0.2;   % �������
%               etam = 1;   % cross index??????????????
%             del=min((Offspring-Lower),(Upper-Offspring))./(Upper-Lower);
%             t=rand(1,D);
%             loc_mut=t<pm;        
%             u=rand(1,D);  % V�Ǳ�����
%             
%             
%             delq=(u<=0.5).*((((2*u)+((1-2*u).*((1-del).^(etam+1)))).^(1/(etam+1)))-1)+...
%                 (u>0.5).*(1-((2*(1-u))+(2*(u-0.5).*((1-del).^(etam+1)))).^(1/(etam+1)));
%             c=Offspring+delq.*loc_mut.*(Upper-Lower);
%             Offspring=c;              
            

    end
    if calObj
        Offspring = INDIVIDUAL(Offspring);
    end
    end