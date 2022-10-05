function x=ar1(n,a,c)

% ------------------------------------------------------------------------
% Input:
%        n=time series length
%        a=parameter
%        c=constant
% Output:
%        x=numerical realization of length n for an autoregressive model of
%        order 1 with parameter a and constant c
% Model:
%        x_{t+1}=c+a*x_{t}+e_{t} with e_{t}~N(0,1)
% Example:
%        x=ar1(10^4,0.8,8)
% ------------------------------------------------------------------------

% Initialize
x=zeros(n,1);
x0=rand(1,1);
x(1)=c+a*x0+randn(1,1);

% Simulate
for i=2:n
    x(i)=c+a*x(i-1)+randn(1,1);
end