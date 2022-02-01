function [x, it] = mycg(A,b,x0,myeps,maxit)
%  function [x,it] = mycg(A,b,x0,myeps,maxit)
%
%  Conjugate gradient method
%
%  x  ... solution
%  it ... number of performed iterations
%
%  A  ... system matrix
%  b  ... right-hand side vector
%  x0 ... initial approximation
%  myeps ... precision in stopping criteria |Ax - b| < myeps
%  maxit ... maximum number of iterations
%

x = x0; % set initial approximation
g = A*x-b; % compute gradient
p = g; % first A-orthonormal vector

gg = dot(g,g);

it = 0; % initial number of iterations

while and(norm(g) > myeps, it < maxit)
    Ap = A*p;
    pAp = dot(Ap,p);

    alpha = gg/pAp;
    x = x - alpha*p;
    g = g - alpha*Ap;
    
    gg_old = gg;
    gg = dot(g,g);
    
    beta = gg/gg_old;
    p = g + beta*p;  

    it = it + 1;
end

end
