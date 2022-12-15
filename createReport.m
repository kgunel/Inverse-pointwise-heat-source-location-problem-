function createReport(Sol,minInd,maxInd,fmin,T,fid)

% Check underflowMeanFlag error
tol = 240;
underflowMeanFlag=0; underflowStdFlag=0;
overflowMeanFlag=0; overflowStdFlag=0;
if mean(fmin) < 10^(-tol)
    underflowMeanFlag = 1;
end
if std(fmin) < 10^(-tol)
    underflowStdFlag = 1;
end

if mean(fmin) > 10^tol
    overflowMeanFlag = 1;
end
if std(fmin) > 10^tol
    overflowStdFlag = 1;
end


if (underflowMeanFlag==0) && (overflowMeanFlag==0)
    MeanCost = num2str(mean(fmin));
else
    if (underflowMeanFlag==1)
        s = num2str(mean(10^tol*fmin));
        ind = strfind(s,'e');
        if ~isempty(ind)
                realPart = [s(1:6) 'e'];
                exponentialPart = num2str(str2double(s(ind+1:end))-tol);
                MeanCost = [realPart exponentialPart];
        else
            MeanCost ='0.000e+000';
        end
    else % Overflow error
        s = num2str(mean(10^(-tol)*fmin));
        ind = strfind(s,'e');
        if ~isempty(ind)
                realPart = [s(1:6) 'e'];
                exponentialPart = num2str(str2double(s(ind+1:end))+tol);
                MeanCost = [realPart exponentialPart];
        else
            MeanCost ='0.000e+000';
        end
    end
end


if (underflowStdFlag==0)  && (overflowStdFlag==0)
    StdCost = num2str(std(fmin));
else
    if (underflowStdFlag==1)
        s = num2str(std(10^tol*fmin));
        ind = strfind(s,'e');
        if ~isempty(ind)
            realPart = [s(1:6) 'e'];
            exponentialPart = num2str(str2double(s(ind+1:end))-tol);
            StdCost = [realPart exponentialPart];
        else
            StdCost ='0.000e+000';
        end
    else % Overflow
        s = num2str(std(10^(-tol)*fmin));
        ind = strfind(s,'e');
        if ~isempty(ind)
            realPart = [s(1:6) 'e'];
            exponentialPart = num2str(str2double(s(ind+1:end))+tol);
            StdCost = [realPart exponentialPart];
        else
            StdCost ='0.000e+000';
        end
    end
end

% Display result in command windows
fprintf('\n\t\t Best Cost : %4.3e\n\t\t Worst Cost : %4.3e\n\t\t Mean of Costs : %s %c %s\n',...
    Sol(minInd).fmin,Sol(maxInd).fmin,MeanCost,char(177),StdCost);
fprintf('\t\t # of Cost function evaluation : %d\n', Sol(minInd).nFeval);
fprintf('\t\t # of Expansion point evaluation : %d\n',Sol(minInd).nExp);
fprintf('\t\t # of Inside contraction point evaluation : %d\n',Sol(minInd).nIC);
fprintf('\t\t # of Outside contraction point evaluation : %d\n',Sol(minInd).nOC);
fprintf('\t\t # of Shrink point evaluation : %d\n',Sol(minInd).nShrink);
fprintf('\t\t Mean of elapsed time : %4.3e %c %4.3e\n',mean(T),char(177),std(T));

% Add results into Results.txt
fprintf(fid,'\n\t\t Best Cost : %4.3e\n\t\t Worst Cost : %4.3e\n\t\t Mean of Costs : %s %c %s\n',...
    Sol(minInd).fmin,Sol(maxInd).fmin,MeanCost,char(177),StdCost);
fprintf(fid,'\t\t # of Cost function evaluation : %d\n', Sol(minInd).nFeval);
fprintf(fid,'\t\t # of Expansion point evaluation : %d\n',Sol(minInd).nExp);
fprintf(fid,'\t\t # of Inside contraction point evaluation : %d\n',Sol(minInd).nIC);
fprintf(fid,'\t\t # of Outside contraction point evaluation : %d\n',Sol(minInd).nOC);
fprintf(fid,'\t\t # of Shrink point evaluation : %d\n',Sol(minInd).nShrink);
fprintf(fid,'\t\t Mean of elapsed time : %4.3e %c %4.3e\n',mean(T),char(177),std(T));
end