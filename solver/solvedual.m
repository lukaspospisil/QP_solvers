function [x,it] = solvedual(K,f,BI,BE,x0)
%SOLVEDUAL Summary of this function goes here
%   Detailed explanation goes here

[U,D] = eig(K);
lambda = diag(D);

% Kernel
lambdazero = 1e-6;
R = U(:,abs(lambda) < lambdazero); 

% pseudoinverse
lambdapinv = zeros(size(lambda));
lambdapinv(abs(lambda) >= lambdazero) = 1./lambda(abs(lambda) >= lambdazero);
Kpinv = U*diag(lambdapinv)*U';

% pinv test
if false
    disp('PSEUDOINVERSE TEST')
    disp(['- err1 = ' num2str(norm(Kpinv - Kpinv*K*Kpinv,'fro'))])
    disp(['- err2 = ' num2str(norm(K - K*Kpinv*K,'fro'))])
end

B = [BE;BI];
mE = size(BE,1);
mI = size(BI,1);

% indexes (what is what)
idxE = 1:mE;
idxI = mE+1:mE+mI;

F = B*Kpinv*B';
G = R'*B';
d = B*Kpinv*f;
e = R'*f;

Q = G'*((G*G')\G);
%Q = G'*inv(G*G')*G;
P = eye(size(Q,1)) - Q;

rho = 0; % regularization parameter?
lambda_tilde = G\e;

% dualization
Adual = P*F*P + rho*Q;
bdual = P*(d - F*lambda_tilde);
BEdual = G;
cEdual = zeros(size(G,1),1);
lbdual = -Inf*ones(size(bdual)); 
lbdual(idxI) = -lambda_tilde(idxI);
%[lambda,~,~,output,lagmult] = quadprog(Adual,-bdual,[],[],BEdual,cEdual,lbdual,[],[]);
%alpha = lagmult.eqlin;

myopts = SmalseOptions();
myopts.normA = gersgorin(Adual);
myopts.normBB = gersgorin(BEdual*BEdual');

x0 = ones(size(bdual));
myopts.maxit = 50;
myopts.dispdebug = false;
[lambda, it, alpha] = smalse_m(Adual,bdual,x0,BEdual,lbdual,myopts);

disp(['eq_dual   = ' num2str(norm(BEdual*lambda,2))])
disp(['ineq_dual = ' num2str(norm(min(lambda(idxI),0),2))])

% reconstruction
x = Kpinv*(f - B'*(lambda+lambda_tilde)) - R*alpha;

end

