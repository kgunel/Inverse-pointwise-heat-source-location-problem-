function IC = insideContraction(gamma,M,R)
  n = length(gamma);
  nVar = length(M);
  IC = zeros(n,nVar);
  for i=1:n
     IC(i,:) = M - gamma(i)*(R-M);
  end
end