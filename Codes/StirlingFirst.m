function s = StirlingFirst(n,k)
%
% This function computes signed Stirling numbers of the 1st kind, as described by Wikipedia.
%
% Input:
%   n = argument 1
%   k = argument 2
%
% Output:
%   s = s(n,k)
%
%       where s is the Stirling number of the first kind, defined by the recurrence relation:
%
%       s(n,k) = -(n-1) * s(n-1 , k) + s(n-1 , k-1)
%
%       with the initial conditions
%
%       s(0,0) = 1,   s(n,0) = s(0,k) = 0, if n > 0
%
% See the Wikipedia article: https://en.wikipedia.org/wiki/Stirling_numbers_of_the_first_kind

%% Check for the terminating conditions
if n == 0
    if k == 0
        s = 1;
    else %k > 0
        s = 0;
    end
else %n > 0
    if k == 0
       s = 0; 
    else %k > 0
        %this is the general case ==> use the recursive formula
        s = -(n-1)*StirlingFirst(n-1,k) + StirlingFirst(n-1,k-1);
    end
end