function KL = ComputeKLDivergence(p,q,tol)
% Computes symmetric KL Divergence
% Assumes that p and q have same domain, otherwise makes no sense
% tol = parameter saying when to ignore small values
if length(p) ~= length(q)
    error('Invalid pair for KL divergence');
end
len = length(p);
pDiv = 0; qDiv = 0;
for i = 1:len
    if min(p(i),q(i)) > tol
        pDiv = pDiv+p(i)*log(p(i)/q(i));
        qDiv = qDiv +q(i)*log(q(i)/p(i));
    end
end
KL = (pDiv+qDiv)/2;
end