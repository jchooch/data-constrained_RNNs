function p = getScore(dataset,ta,tb,s)
% Function to get score for partial models
% This function computes score for the partial model, which is indicative
% of the significance of the partial model. This function helps in
% selecting best 'l' models.
% 
% The method used for calculating the score is as follows:
% For s=2  : 
%           score = k/n where k is the total number of planted rows found
% otherwise:
%           score is obtained by solving for 'p' in the following equation:
%           For i=1:n : A(i)/((A(i)-B)*p + B) = n
%           The equation can be solved using newton iteration method.
%
% The definitions for A(i) and B can be found in the paper mentioned below:
% Ben-Dor, Amir, Benny Chor, Richard Karp, and Zohar Yakhini. "Discovering
% local structure in gene expression data: the order-preserving submatrix problem." 
% In Proceedings of the sixth annual international conference on
% Computational biology, pp. 49-57. ACM, 2002.
%            
% Inputs:
%        dataset      :      Input matrix [n*m] of co-occurences of n 
%                            instances (rows) and m features (columns).
%        ta           :      Column with lowest rank.
%        tb           :      Column with highest rank.
%        s            :      Number of columns being considered at a time.
%
% Outputs:
%         p           :      Significance score for partial models.
% 
% Author: Sumanik Singh, 2013
%        
% Contact: sumanik@iitk.ac.in, sumaniksingh@gmail.com
%          Department of Electrical Engineering, 
%          Indian Institute of Technology, Kanpur, India
%% Implementation of 'getScore' function
n = size(dataset,1);
m = size(dataset,2);
% Finding the ranked matrix D for input dataset
[so,ind] = sort(dataset,2);
[s1,D] = sort(ind,2);

% Calculate exact value of p if s=2 because the model is already complete.
if s==2
    k=0;
    for row = 1:size(dataset,1)
        if D(row,ta) < D(row,tb)
            k=k+1;
        end
    end
    p = k/n;
% else do newton iteration to estimate p
else
    gapSize = s-2;
    commonDenominator = nchoosek(m,s);
    % Calculating B[i] or B
    B = 1/(nchoosek(m,2)*2);
    % Calculating A[i] and diff[i] = A[i]-B
    for i=1:n
        gi = D(i,tb) - D(i,ta) - 1;
        if gi >= gapSize
            A(i) = nchoosek(gi,gapSize)/commonDenominator;
        else
            A(i)=0;
        end
        diff(i) = A(i) - B;
    end
    
    % Starting Point for Newton Iteration
    p = 0.05;              
    % Do at most 20 newton iteration steps to get p
    for iteration = 1:20
        f = -n;
        dfdp = 0;
        for i = 1:n
            denominator(i) = diff(i)*p + B;
            if denominator(i)==0 && A(i)~=0
                disp('WARNING: Division by 0 in Newton Iteration.');
                f = 0;
                dfdp = 1;
            else if A(i)~=0
                    f = f + (A(i)/denominator(i));
                    dfdp = dfdp - (A(i)*diff(i)/(denominator(i)^2));
                end
            end
        end
        
        if dfdp==0
            p = 0;
        else
            q = p;
            p = p - (f/dfdp);
            if abs(f/dfdp) < 0.0001
                break;
            end
            if p<0 && iteration<20
                p = q/2;
            end
        end
    end
    % special modification of newton algorithm, required to prevent from being caught in negative range 
    if p<0 || p>1
        p=0;
    end
end
end
