function fitFC=compareconn(A,B)
N=length(A);
Isubdiag = find(tril(ones(N),-1));
fitFC=corr2(A(Isubdiag),B(Isubdiag));