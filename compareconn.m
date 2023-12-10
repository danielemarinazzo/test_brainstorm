function fitFC=compareconn(A,B)
N=length(A);
Isubdiag = find(tril(ones(N),-1));
fitFC=corr(A(Isubdiag),B(Isubdiag),'type','Spearman');