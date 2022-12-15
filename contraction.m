function [IC, OC] = contraction(gamma,M,W)
n = length(gamma);
nVar = length(M);
IC =zeros(n,nVar);
OC =zeros(n,nVar);
for i=1:n
    IC(i,:) = insideContraction(gamma(i),M,W);
    OC(i,:) = outsideContraction(gamma(i),M,W);
end
end