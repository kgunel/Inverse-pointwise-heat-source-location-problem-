function S = shrink(delta,X,Y)
  n = length(delta);
  nVar = length(X);
  S = zeros(n,nVar);
  for i=1:n
    S(i,:) = X + delta(i)*(Y - X);
  end
end