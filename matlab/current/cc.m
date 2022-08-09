%%% THIS FUNCTION IS FROM THE MTBA LIBRARY DUE TO J. K. Gupta, S. Singh and N. K. Verma
%%% SEE MORE INFORMATION AT: https://www.iitk.ac.in/idea/mtba/

function biClustResult = cc(datamatrix, alpha, delta, number, random)
%Cheng&Church Algorithm
% Cheng, Yizong, and George M. Church. 
% "Biclustering of expression data." 
% Proceedings of the eighth international conference on intelligent systems 
% for molecular biology. Vol. 8. 2000.
%
% Usage
% >> biClustResult =  bimaxBiclust(mat, alpha, delta, numBiclust, rand)
%
% Inputs:
%   datamatrix            - the data matrix (logical)
%   alpha(optional)       - Scaling factor
%   delta(optional)       - Maximum of accepted scores
%   numBiclust(optional)  - number of biclusters to be found
%
% Outputs:
%   biClustResult: A structure consisting of
%       RowxNum     - Logical Matrix which contains 1 in [i,j] if Row i is 
%                     in Bicluster j
%       NumxCol     - Logical Matrix which contains 1 in [i,j] if Col j is 
%                     in Bicluster i
%       ClusterNo   - Number of clusters
%       Clust       - Another structure array containing all clusters with
%                     their respective row and column indices.
%
% Author: Jayesh Kumar Gupta, 2013.
%
% Contact: Jayesh Kumar Gupta http://home.iitk.ac.in/~jayeshkg
%          Indian Institute of Technology, Kanpur, India

%##########################################################################
%       Error Checking and Default Values
%##########################################################################
if nargin < 1
  error('input :  No matrix as input');
end

if nargin < 2
  alpha   = 1.5;
  delta   = 1.0;
  number  = 100;
end

if nargin < 5
  random  = false;
end
%##########################################################################
mat   = datamatrix;
ma    = max(mat(:)); % Useful when random thing done
mi    = min(mat(:));
nrow  = size(mat, 1);
ncol  = size(mat, 2);
x     = false(nrow, number);
y     = false(number, ncol);
STOP  = false;
logr  = true(nrow,1);

for i = 1:number
  if sum(logr)<2
    STOP = true;
    break;
  end
  
  erg = bigcc(mat(logr,:), delta, alpha);
  
  if sum(erg.logr)==0
    STOP = true;
    break;
  else
    x(logr,i) = erg.logr;
    y(i,:) = erg.logc;
    
    if random
      mat(erg.logr, erg.logc) = unifrnd(repmat(mi, ...
        length(find(erg.logr)), length(find(erg.logc))), ...
        repmat(ma, length(find(erg.logr)), length(find(erg.logc)))); % Something Horribly wrong is here
    else
      logr(logr) = ~(logr(logr)&erg.logr);
    end
  end
end

if STOP
  biClustResult.RowxNum   = x(:,1:(i-1));
  biClustResult.NumxCol   = y(1:(i-1),:);
  biClustResult.ClusterNo = i-1;
  for j=1:(i-1)
    rows = find(biClustResult.RowxNum(:,j)>0);
    cols = find(biClustResult.NumxCol(j,:)>0);
    biClustResult.Clust(j) = struct('rows', rows, 'cols', cols); % Completely different thing from in other functions
  end
else
  biClustResult.RowxNum = x;
  biClustResult.NumxCol = y;
  biClustResult.ClusterNo = i;
  for j=1:i
    rows = find(biClustResult.RowxNum(:,j)>0);
    cols =find(biClustResult.NumxCol(j,:)>0);
    biClustResult.Clust(j) = struct('rows', rows, 'cols', cols); % Completely different thing from in other functions
  end
end

end
%##########################################################################

%##########################################################################
% Find biggest bicluster
%##########################################################################
function res = bigcc(mat, delta, alpha)
logr = true(size(mat,1),1);
logc = true(size(mat,2),1);

step1 = cc2(mat, logr, logc, delta, alpha);
step2 = cc1(mat, step1.logr, step1.logc, delta);

if sum(step2.logr)==0
  res = struct('logr',0,'logc',0);
  warning(sprintf('No matrix with score smaller than delta=%d found',...
    delta));
else
  res = cc3(mat, step2.logr, step2.logc);
end
end
%##########################################################################

%##########################################################################
% Algorithm-1 :: CC Single Node Deletion
%##########################################################################
function res = cc1(mat, logr, logc, delta)
while (ccscore(mat(logr,logc))>delta)
  di = rowscore(mat(logr,logc));
  dj = colscore(mat(logr,logc));
  [~,maxdi] = max(di);
  [~,maxdj] = max(dj);
  mdi = false(size(logr(logr))); %}
  mdi(maxdi)=true;         %} Little hackish here
  mdj = false(size(logc(logc))); %}
  mdj(maxdj)=true;        %}
  
  if di(maxdi)>dj(maxdj)
    logr(logr) = ~(logr(logr)&mdi);
  else
    logc(logc) = ~(logc(logc)&mdj);
  end
  
  if (~(sum(logr)>1 && sum(logc)>1))
    break;
  end
end
if (sum(logr)>1 && sum(logc)>1)
  res = struct('logr',logr,'logc',logc);
else
  res = struct('logr',0,'logc',0);
  warning(sprintf('No matrix with score smaller than %d found', delta));
end
end
%##########################################################################

%##########################################################################
% Algorithm-2 :: CC Multiple Node Deletion
%##########################################################################
function res = cc2(mat, logr, logc, delta, alpha)
mdi = 1;
mdj = 1;
h=ccscore(mat((logr),(logc)));
while (h>delta && (sum(mdi)+sum(mdj))>0)
  if sum(logr)>100
    di = rowscore(mat(logr,logc));
    mdi = di>(alpha*h);
    if sum(mdi)<(sum(logr)-1)
      %             logr = logr(logr);
      %             logr(mdi) = false; %Error is here
      logr(logr) = ~(logr(logr)&mdi);
      h=ccscore(mat(logr,logc));
    else
      warning(sprintf('Alpha:%d\t Too Small',alpha));
      mdi = 0;
    end
  else
    mdi=0;
  end
  if sum(logc)>100
    dj = colscore(mat(logr,logc));
    mdj = dj>(alpha*h);
    if sum(mdj)<(sum(logc)-1)
      logc(logc) = ~(logc(logc)&mdj);
    else
      warning(sprintf('Alpha:%d\t Too Small',alpha));
      mdj = 0;
    end
  else
    mdj=0;
  end
  h=ccscore(mat((logr),(logc)));
end
res = struct('logr',logr,'logc',logc);
end
%##########################################################################

%##########################################################################
% Algorithm-3 :: CC Node Addition
%##########################################################################
function res = cc3(mat, logr, logc)
br=1;
ilogr = false(size(logr));
while br>0
  br1 = sum(logc);
  br2 = sum(logr);
  h   = ccscore(mat(logr,logc));
  dj  = addcolscore(mat, logr, logc);
  
  mdj = dj<=h;
  logc(mdj) = true;
  
  h   = ccscore(mat(logr,logc));
  di  = addrowscore(mat,logr,logc);
  idi = inverseaddrowscore(mat,logr,logc);
  
  mdi = di<=h;
  logr(mdi) = true;
  imdi = (idi<=h);
  mat((~(bsxfun(@eq,logr,imdi)))&imdi,:) = -mat((~(bsxfun(@eq,logr,...
    imdi)))&imdi,:);
  logr(imdi) = true;
  
  br=sum(logc)+sum(logr)-br1-br2;
end
res = struct('logr',logr,'logc',logc);
end
%##########################################################################


%##########################################################################
% Helper functions to Calculate scores for Node Deletion
%##########################################################################
function rowScore = rowscore(mat)
% Subtract row mean
a1 = bsxfun(@minus, mat, mean(mat,2));
% Subtract col mean
a2 = bsxfun(@minus, a1, mean(mat));
% Add matrix mean
a3 = bsxfun(@plus, a2, mean2(mat));
% Sum along columns then division by number of columns
rowScore = sum(a3.^2, 2)/size(mat,2);
end

function colScore = colscore(mat)
% Subtract row mean
a1 = bsxfun(@minus, mat, mean(mat,2));
% Subtract col mean
a2 = bsxfun(@minus, a1, mean(mat));
% Add matrix mean
a3 = bsxfun(@plus, a2, mean2(mat));
% Sum along rows then divide by number of rows
colScoreT = sum(a3.^2,1)/size(mat,1);
colScore = colScoreT';
end

% Calculate H(I,J)
function ccScore = ccscore(mat)
% Subtract row mean
a1 = bsxfun(@minus, mat, mean(mat,2));
% Subtract col mean
a2 = bsxfun(@minus, a1, mean(mat));
% Add matrix mean
a3 = bsxfun(@plus, a2, mean2(mat));
ccScore = sumsqr(a3)/(size(mat,1)*size(mat,2));
end

%##########################################################################
% Helper functions to Calculate scores for Node Addition
%##########################################################################
function addrowScore = addrowscore(mat, row, col)
% Subtract row mean
a1 = bsxfun(@minus, mat, mean(mat(:,col),2));
% Subtract col mean
a2 = bsxfun(@minus, a1, mean(mat(row,:)));
% Add matrix mean
a3 = bsxfun(@plus, a2, mean2(mat(row,col)));
% Sum along columns then division by number of columns
addrowScore = sum(a3.^2, 2)/size(mat(row,col),2);
end

function addcolScore = addcolscore(mat, row, col)
% Subtract row mean
a1 = bsxfun(@minus, mat, mean(mat(:,col),2));
% Subtract col mean
a2 = bsxfun(@minus, a1, mean(mat(row,:)));
% Add matrix mean
a3 = bsxfun(@plus, a2, mean2(mat(row,col)));
% Sum along rows then divide by number of rows
addcolScore1 = sum(a3.^2,1)/size(mat(row,col),1);
addcolScore = addcolScore1';
end

function iaddrowscore = inverseaddrowscore(mat, row, col)
% Add row mean
a1 = bsxfun(@plus, -mat, mean(mat(:,col),2));
% Subtract col mean
a2 = bsxfun(@minus, a1, mean(mat(row,:)));
% Add matrix mean
a3 = bsxfun(@plus, a2, mean2(mat(row,col)));
% Sum along columns then division by number of columns
iaddrowscore = sum(a3.^2, 2)/size(mat(row,col),2);
end
