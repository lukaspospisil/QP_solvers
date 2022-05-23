clear all

addpath('solver/')
addpath('solver/mprgp/')
addpath('solver/smalse/')

n = 10;
mE = 5;
mI = 2;

BE = randn(mE,n);
BI = randn(mI,n);

A = randn(n,n);
A = A'*A;
A = A - min(eig(A))*eye(n);
b = randn(n,1);

x0 = randn(n,1);

%quadprog
disp('------------------------------------')
x1 = quadprog(A,-b,BI,zeros(size(BI,1),1),BE,zeros(size(BE,1),1),[],[],x0);
disp(['f = ' num2str(0.5*x1'*A*x1 - b'*x1)])
disp(['eq = ' num2str(norm(BE*x1,2))])
disp(['ineq = ' num2str(norm(max(BI*x1,0),2))])

%Dostal
disp('------------------------------------')
[x2,it2] = solvedual(A,b,BI,BE,x0);
disp(['f    = ' num2str(0.5*x2'*A*x2 - b'*x2)])
disp(['eq   = ' num2str(norm(BE*x2,2))])
disp(['ineq = ' num2str(norm(max(BI*x2,0),2))])
disp(['it   = ' num2str(it2)])

