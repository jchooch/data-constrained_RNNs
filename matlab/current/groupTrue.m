function labelmap = groupTrue(xbinary)
%Takes a logical input vector and returns a label map of all the consecutive 
%true elements. False elements arelabeled zero.
% 
%       labelmap = groupTrue(xbinary)       
%                            
%IN:                         
%                            
%    xbinary: an input vector of type logical.               
%               
%OUT:                        
%                            
%    labelmap: same as xbinary but with consecutive trues replaced with numerical
%              group labels.
%
%EXAMPLE:
% 
%     >> xbinary=logical([1 0 0 1 1 1 0 1 1 0]);
%     >> labelmap = group1s(xbinary)
% 
%      labelmap =
% 
%          1     0     0     2     2     2     0     3     3     0

 assert(islogical(xbinary)&&isvector(xbinary),'Input must be a logical vector')
 
 y=reshape(diff([0,xbinary(:).']),size(xbinary));
 y(y<=0)=0;
 
 labelmap=cumsum(y).*xbinary;

end