function normMat=NormalizeNetwork(X);

mu0=mean(X(:));
std0=std(X(:));
mu1=mean(X);
std1=std(X,1);
mu2=mean(X,2);
std2=std(X,1,2);

Z1=(X-repmat(mu1, size(X,1), 1))./repmat(std1, size(X,1), 1);
Z2=(X-repmat(mu2, 1, size(X,2)))./repmat(std2, 1, size(X,2));
normMat=Z1/sqrt(2)+Z2/sqrt(2);

% checks and defaults for missing data
Z0=(X-mu0)/std0;
f1=isnan(Z1); f2=isnan(Z2);
normMat(f1)=Z2(f1)/sqrt(2)+Z0(f1)/sqrt(2);
normMat(f2)=Z1(f2)/sqrt(2)+Z0(f2)/sqrt(2);
normMat(f1 & f2)=2*Z0(f1 & f2)/sqrt(2);
