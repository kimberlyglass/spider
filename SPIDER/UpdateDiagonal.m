function diagMat=UpdateDiagonal(diagMat, num, alpha, step);

diagMat(1:(num+1):end)=nan;
% diagstd=nanstd(diagMat,1); % this line works only in MATLAB, replaced with the two below for compatibility with Octave
DM=reshape(diagMat(~isnan(diagMat)), num-1, num);
diagstd=std(DM,1);
diagMat(1:(num+1):end)=diagstd*num*exp(2*alpha*step);
