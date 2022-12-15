function y = Cost(position, pfunc, bfunc, measuredT, x, y, low, up, nx, ns, epss)
   candidate_source_x = position(1);
   candidate_source_y = position(2);
   % T_tilde is the temperature calculated using candidate heat source
   [T_tilde, x_tilde, y_tilde] = fd2poisson(pfunc, bfunc, candidate_source_x, candidate_source_y, low, up, nx, epss);
   %[T_tilde, ~, ~] = fd2poisson(pfunc, bfunc, candidate_source_x, candidate_source_y, low, up, nx, epss);
   [m, n]= size(T_tilde);
   totalE = 0; %counter = 0;
   for i=2:m-1
       for j=2:n-1
             totalE = totalE + (T_tilde(i,j) - measuredT((i-1)*ns+1,(j-1)*ns+1))^2 ;
       end
   end
y = totalE ;
end