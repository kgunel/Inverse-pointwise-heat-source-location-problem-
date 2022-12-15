function y = stepSearch(fhandle, objfnc, X, Y,  a, b, tol)
fhandle = fcnchk(fhandle);
stepSize = 0.1;
n = 2;

% Golden section search
x1 = b - stepSize*(b - a);
x2 = b + stepSize*(b - a);
while abs(b-a)>= tol
    params = [x1 x2];
    if isequal(fhandle,@contraction)
        [IC, OC] = feval(fhandle, params, X, Y);
        for i=1:n
            fIC(i) = feval(objfnc,IC(i,:));
            fOC(i) = feval(objfnc,OC(i,:));
            if fIC(i)<fOC(i)
                fParams(i,:) = IC(i,:);
            else
                fParams(i,:) = OC(i,:);
            end
        end
    else
        fParams = feval(fhandle, params, X, Y);
    end
    
    for i=1:2
       Costs(i) = feval(objfnc,fParams(i,:));
    end
    
    if Cost(2)>Cost(1) % f(x2)>f(x1)
        b = x2;
        x2 = x1;
        x1 = b - stepSize*(b - a);
    else
        a = x1;
        x1 = x2;
        x2 = b + stepSize*(b - a);
    end
end
if Cost(1) < Cost(2)
    y = x1;
else
    y = x2;
end