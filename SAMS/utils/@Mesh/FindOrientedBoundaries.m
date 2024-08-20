function [BL,BI] = FindOrientedBoundaries(G)
%Computes the list of boundaries in order
%BL lists the boundaries in order
%BI lists the next boundary in the chain

[BV,BE] = G.FindBoundaries;
nv = size(BV);
BL = [];
BI = [];
bn = 1;
while nnz(BV) > 0
    vind = find(BV,1);
    v=BV(vind);
    BL(bn,:)=0;
    BE(v,2)=0;
    i=1;
    while v ~= 0
        BL(bn,i) = v;
        i=i+1;
        vind=find(BV==v);
        BV(vind)=0;
        BE(BE==v)=0;
        v=BE(v,BE(v,:)~=0);
    end
    bn=bn+1;
end

for i=1:bn-1
    BLi=BL(i,:);
    ind=find(BLi);
    BI{i}=zeros(size(ind));
    BI{i}(BLi(ind)) = ind;
end