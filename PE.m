function entropy=PE(ts,wl,lag)

% ------------------------------------------------------------------------
% Input:
%        ts=time series
%        wl=order for the ordinal symbolization of ts
%        lag=lag for the ordinal symbolization of ts
% Output:
%        entropy=permutation entropy
% Example:
%        ts=randn(1000,1);
%        PE_randn=PE(ts,3,1)
% ------------------------------------------------------------------------

m=length(ts)-(wl-1)*lag;
pk=hist_indices(ts,wl,lag)/m;
pk=pk(pk~=0);
entropy=-dot(pk,log2(pk))/(wl-1);

function hist_indcs=hist_indices(ts,wl,lag)
indcs=perm_indices(ts,wl,lag);
hist_indcs=hist(indcs,1:1:factorial(wl));

function indcs=perm_indices(ts,wl,lag)
m=length(ts)-(wl-1)*lag;
indcs=zeros(m,1);
for i=1:wl-1
st=ts(1+(i-1)*lag :m+(i-1)*lag);
for j=i:wl-1
indcs=indcs+(st>ts(1+j*lag:m+j*lag));
end
indcs=indcs*(wl-i);
end
indcs=indcs + 1;