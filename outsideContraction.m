function OC = outsideContraction(gamma,M,R)
  n = length(gamma);
  nVar = length(M);
  OC = zeros(n,nVar);
  for i=1:n
     OC(i,:) = M + gamma(i)*(R - M);
  end
end