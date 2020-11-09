function index = TournamentSelection(K,N,varargin)
%TournamentSelection - Tournament selection.
%
%   P = TournamentSelection(K,N,fitness1,fitness2,...) returns the indices
%   of N solutions by K-tournament selection based on their fitness values.
%   In each selection, the candidate having the MINIMUM fitness1 value will
%   be selected; if more than one candidates have the same minimum value of
%   fitness1, then compare their fitness2 values, and so on.
%
%   Example:
%       P = TournamentSelection(2,100,FrontNo)

%------------------------------- Copyright --------------------------------

%--------------------------------------------------------------------------

    varargin = cellfun(@(S) reshape(S,[],1),varargin,'UniformOutput',false);
    [~,rank] = sortrows([varargin{:}]);
    [~,rank] = sort(rank);
    Parents  = randi(length(varargin{1}),K,N);
    [~,best] = min(rank(Parents),[],1);
    index    = Parents(best+(0:N-1)*K);
end