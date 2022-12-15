function E = expansion(beta,M,R)
  %n = length(beta);
  nVar = length(M);
  E =zeros(1,nVar);
  %for i=1:n
     %E(i,:) = M + beta(i)*(R - M);
  %end
  E = M + beta.*(R - M);
end