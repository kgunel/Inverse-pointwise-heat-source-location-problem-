function [xmin,fmin,nFeval,nExp,nIC,nOC,nShrink,iter,BestCost]=NMS(objfnc, nVar, VarMin, VarMax, maxIter, maxnFeval, tol, u, x, y, varEps)
%
% Nelder-Mead Simplex Algorithm for
% solving the unconstrained optimization problem:
%         min f(x).
%
% It uses the adaptive parameters introduced in the following paper:
%
%
% Inputs:
%  xinit--initial guess
%  tol--tolerance for termination (Recommended value: 10^-4)
%  maxnFeval--maximum number of function evaluations
%  myfunction--objective function. (In myfunction(x), x is a ROW vector)
%
% Outputs:
%  xmin--approximate optimal solution at termination. It is a row vector.
%  fmin--minimum function value at termination
%  ct--number of function evaluations at termination
%

%% Settings
ns = 25;                        % Sub-mesh size for measured T
nx = floor(size(x,1)/ns)-1;     % Mesh size for measured T
pfunc = 'f';
bfunc = 'g';
measuredT = u;

 % set up parameters
 alpha=1; beta = 2; gamma=0.5; delta=0.5;


dim = nVar;
X = VarMin + (VarMax - VarMin)*rand(dim+1,nVar);
for i=1:dim+1
     FX(i)=feval(objfnc, X(i,:), pfunc, bfunc, measuredT, x, y, VarMin, VarMax, nx, ns, varEps);
end
ct=dim+1;
[FX,I]=sort(FX);
X=X(I,:);

%% Main iteration
BestCost = zeros(1,maxIter);
nExp=0; nIC=0; nOC=0; nShrink=0;
%gcf=figure;
for iter = 1:maxIter
%      if max(max(abs(X(2:dim+1,:)-X(1:dim,:)))) <= tol
%          break;
%      end
%     
%     if ct>maxnFeval
%         break;
%     end
  
    % Centroid of the dim best vertices
    M = mean(X(1:dim,:));  
    
    xref=(1+alpha)*M- alpha*X(dim+1,:);
    Fref=feval(objfnc, xref, pfunc, bfunc, measuredT, x, y, VarMin, VarMax, nx, ns, varEps);
    ct=ct+1;
    if Fref<FX(1)
        % expansion
        xexp=(1+alpha*beta)*M-alpha*beta*X(dim+1,:);
        Fexp=feval(objfnc, xexp, pfunc, bfunc, measuredT, x, y, VarMin, VarMax, nx, ns, varEps);
        ct=ct+1;
        nExp=nExp+1;
        if Fexp < Fref
            X(dim+1,:)=xexp;
            FX(dim+1)=Fexp;
        else
            X(dim+1,:)=xref;
            FX(dim+1)=Fref;
        end
    else
        if Fref<FX(dim)
            % accept reflection point
            X(dim+1,:)=xref;
            FX(dim+1)=Fref;
        else
            if Fref<FX(dim+1)
                % Outside contraction
                xoc=(1+alpha*gamma)*M-alpha*gamma*X(dim+1,:);
                Foc=feval(objfnc, xoc, pfunc, bfunc, measuredT, x, y, VarMin, VarMax, nx, ns, varEps);
                ct=ct+1;
                nOC = nOC+1;
                if Foc<=Fref
                    X(dim+1,:)=xoc;
                    FX(dim+1)=Foc;
                else
                    % shrink
                    for i=2:dim+1
                        X(i,:)=X(1,:)+ delta*(X(i,:)-X(1,:));
                        FX(i)=feval(objfnc, X(i,:), pfunc, bfunc, measuredT, x, y, VarMin, VarMax, nx, ns, varEps);
                        nShrink=nShrink+1;
                    end
                    ct=ct+dim;
                end
            else
                %inside contraction
                xic=(1-gamma)*M+gamma*X(dim+1,:);
                Fic=feval(objfnc, xic, pfunc, bfunc, measuredT, x, y, VarMin, VarMax, nx, ns, varEps);
                ct=ct+1;
                nIC=nIC+1;
                if Fic<FX(dim+1)
                    X(dim+1,:)=xic;
                    FX(dim+1)=Fic;
                else
                    % shrink
                    for i=2:dim+1
                        X(i,:)=X(1,:)+ delta*(X(i,:)-X(1,:));
                        FX(i)=feval(objfnc, X(i,:), pfunc, bfunc, measuredT, x, y, VarMin, VarMax, nx, ns, varEps);
                        nShrink=nShrink+1;
                    end
                    ct=ct+dim;
                end
            end
        end
    end
    [FX,I]=sort(FX);
    X=X(I,:);
    BestCost(iter) = FX(1);
    if mod(iter,100)==0
       fprintf('Iteration = %4d\t :\t Best Cost = %5.4e\n', iter, BestCost(iter));
    end
end
xmin=X(1,:);
fmin=FX(1);
nFeval=ct;