function biClustResult =  bimax(datamatrix, minNoRows,...
                                    minNoCols, numBiclust)
%bimax function for one time running 
%
% Usage
% >> [RowxNum,NumxCol,Number] =  bimaxBiclust(datamatrix, minNoRows, ...
%                                                   minNoCols, numBiclust)
%
% Inputs: 
%   datamatrix            - the data matrix (logical)
%   minNoRows(optional)   - minimal occurences on each row or column
%   minNoCols(optional)   -
%   numBiclust(optional)  - 
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
    minNoRows = 0;
    minNoCols = 0;
    numBiclust = 100;
end

if nargin < 4
    numBiclust = 100;
end
%##########################################################################

noRows = size(datamatrix,1);
noCols = size(datamatrix,2);
[x,y,~,er]= mbimax(datamatrix, noRows, noCols, minNoRows, ...
                    minNoCols, numBiclust);

if er == 1
    warning('Did not calculate all biclusters');
end

number1 = sum(x,1); % Sum of Each Column
number2 = sum(y,2); % Sum of Each Row

number_sum = number1 + number2';
Number = size(find(number_sum > 0),2); % Actual Number of Biclusters

RowxNum = logical(double(x(:,find(number_sum > 0))));
NumxCol = logical(double(y(find(number_sum > 0),:)));


for i=1:Number
 rows = find(RowxNum(:,i)>0);
 cols =find(NumxCol(i,:)>0);
 cluster(i) = struct('rows', rows, 'cols', cols);
end

% Returning result
biClustResult = struct('RowxNum',RowxNum, 'NumxCol',NumxCol, ...
                        'ClusterNo', Number, 'Clust', cluster);
end
