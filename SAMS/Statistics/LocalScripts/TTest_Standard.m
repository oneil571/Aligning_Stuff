function [f,p]=TTest_Standard(Group1,Group2)
%
% [h p] =hotelT2(v1,v2)
%
% Computes the Hotelling's T-square statistic for given vector measures v1 and v2
% v1,v2    :  vector measures  
%                For instance if v1 is the size of 2562 x 3 x 47, 
%                v1(:,1,12) is the x-coordinate of the 12-th subject.
%                v1 and v2 do not have to have the same number of subjects.
%
%  h : Hotelling's T2 statistic
%  p : The corresponding p-value
%
% (C) Moo K. Chung 2007-2008
%     Department of Biostatistics and Medical Informatics
%     Waisman Laboratory for Brain Imaging
%     University of Wisconsin-Maison
%  
% email://mkchung@wisc.edu

n1 = length(Group1); n2 = length(Group2);
d = size(Group1{1}.V,1);
V_G1 =  zeros(d,Group1{1}.nV,n1);
V_G2 =  zeros(d,Group1{1}.nV,n2);

for k = 1:n1
    V_G1(:,:,k) = Group1{k}.V;
end
for k = 1:n2
    V_G2(:,:,k) = Group2{k}.V;
end
t2 = zeros(1,Group1{1}.nV);

for i=1:Group1{1}.nV
    W_G1=squeeze(V_G1(:,i,:))';
    W_G2=squeeze(V_G2(:,i,:))';
    m1 = mean(W_G1,1)'; m2 = mean(W_G2,1)';
    if n1 < 2
        C1 = zeros(d,d);
    else
        C1 = cov(W_G1);
    end
    if n2 < 2
        C2 = zeros(d,d);
    else
        C2 = cov(W_G2);
    end
    pooledCov = ((n1-1)*C1+(n2-1)*C2)/(n1+n2-2);
    t2(i)=(m2-m1)'*pinv(pooledCov)*(m2-m1);
end

t2 = n1*n2/(n1+n2)*t2;           %T2 statistics, separated only for education
f=(n1+n2-d-1)/(d*(n1+n2-2))*t2;
p=1-fcdf(f, d, n1+n2-d-1 );
