function ridx=roulette_distribution(P,N)
% P is a distribution
% N is the number of values to randomize;

% P=[0.1,0.8,0.1];
% x=[10,20,30];
% N=1000;
 
FP=cumsum(P); % cumulative distribution (CDF)

p=rand(N,1); % randomize uniform variable
ridx=find_indx_upper_in_array(p,FP); % transform to indexes into the CDF
