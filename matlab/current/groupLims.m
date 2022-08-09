function [idxStarts,idxStops,Lengths] = groupLims(X,labelmode)
%Finds the starting indices, stopping indices, and lengths of each group.
%
%   [starts,stops,lengths] = groupLims(X) 
%   [starts,stops,lengths] = groupLims(G,1)
%
%The second syntax assumes that G is a grouping vector resulting from
%either groupConsec() or groupTrue(), meaning that it will ignore groups
%labeled zeros.

 if nargin<2, labelmode=false; end
 
 if isempty(X),[idxStarts,idxStops,Lengths]=deal(0); return; end
 
 iscol=size(X,1)>size(X,2);
 X=X(:).';
 
 if labelmode
    
    G=X; 
     
    dG=diff(G);


    lidxStarts=[G(1)~=0, dG>0];
    idxStarts=find(lidxStarts);    if iscol, idxStarts=idxStarts(:);end
    
    if nargout<2, return; end
    
    lidxStops=[(lidxStarts(2:end) & G(1:end-1)>0)|dG<0,G(end)~=0];
    idxStops=find(lidxStops);  if iscol, idxStops=idxStops(:);end

     
 else

    transitions=diff(X)~=0;

    
    lidxStarts=[1, transitions];
    idxStarts=find(lidxStarts);    if iscol, idxStarts=idxStarts(:);end
    
    if nargout<2, return; end
    
    lidxStops=[transitions, 1];
    idxStops=find(lidxStops);   if iscol, idxStops=idxStops(:);end
     
 end
 
 if nargout>2
  Lengths=idxStops+1-idxStarts;
 end
 
end

