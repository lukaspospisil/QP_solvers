function [x_sol, it, lambda_sol, it_in_out] = smalse_m(A,b,x0,BE,lb,opts)
% SMALE-M algorithm

rho = opts.rho; % rho > 0
beta = opts.beta; % 0 < beta < 1
M0 = opts.M; % M0 > 0

lambda = zeros(size(BE,1),1);
M = M0;

it = 1; % iteration counter

x = x0;

it_in_out = [];

gp_out = Inf;

BExc = BE*x;

%nmb_hess_mult = nmb_hess_mult_mprgp; % counter of hessian multiplication
%it_mprgp_all(it) = it_mprgp; % counter of MPGP iterations

% set settings of inner solver
opts_inner = MprgpOptions();
opts_inner.smalsem = true;
opts_inner.smalsem_obj.eta = opts.eta;
opts_inner.smalsem_obj.BE = BE;

while (norm(BExc) > opts.myeps || norm(gp_out) > opts.myeps) && it < opts.maxit
    x_old = x;
    lambda_old = lambda;

    normA = opts.normA + rho*opts.normBB;
    opts_inner.smalsem_obj.M = M;

    % use MPRGP to solve inner problem
%    [x, it_in] = mprgp(A + rho*BE'*BE, b - BE'*lambda, lb, x, normA, opts_inner);
    [x, it_in, gp_out, hess_mult, gp_norms, fi_norms, beta_norms] = mprgp(A + rho*BE'*BE, b - BE'*lambda, lb, x, normA, opts_inner);
    
    it_in_out(end+1) = it_in;
    
    BExc = BE*x;
    lambda = lambda + rho*(BExc);
    
    if get_L(A,b,x,lambda,rho,BExc) < get_L(A,b,x_old,lambda_old,rho,BE*x_old) + rho/2*dot(BExc,BExc)
        M = beta*M; 
%        rho = rho/beta;
    end

    if opts.dispdebug
       disp('-------------------------------------'); 
       disp(['SMALSE-M it = ' num2str(it)]); 
       disp([' it_in = ' num2str(it_in)])
       disp([' M     = ' num2str(M) ]); 
       disp([' rho   = ' num2str(rho) ]); 
       disp([' L     = ' num2str(get_L(A,b,x,lambda,rho,BExc))])
       disp([' norm(lambda) = ' num2str(norm(lambda)) ]); 
       disp([' norm(BEx) = ' num2str(norm(BE*x)) ]); 
       disp([' norm(qp) = ' num2str(norm(gp_out)) ]); 
%       disp([' it_mprgp = ' num2str(it_mprgp) ]); 
%       verify_kkt_mprgp( A + rho*B'*B,(b - B'*lambda_old),l,x );
  
    end

    it = it + 1; 
    
end

if it == opts.maxit
   disp('!!! error in SMALSE-M iterations '); 
   disp([' outer: ' num2str(norm(BE*x))]);
   disp([' inner: ' num2str(norm(gp_out))]);
end

%it_mprgpout = sum(it_mprgp_all);
x_sol = x; 
lambda_sol = lambda;  

% say goodbye
if(false)
    disp(['outer it: ' num2str(it)])
    disp(['inner it: ' num2str(it_mprgpout)])
    disp(' ')
    
    if it_mprgpout==0
       keyboard 
    end
end
    
end

function [L] = get_L(A,b,x,lambda,rho,Bxc)
    L = get_f(A,b,x) + lambda'*(Bxc) + rho/2*norm(Bxc)^2;
    
end

function [fx] = get_f(A,b,x)
    fx = 1/2*dot(A*x,x) - dot(b,x);
end
