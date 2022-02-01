function [x it hess_mult gp_norms fi_norms beta_norms] = mprgp(A, b, lb, ub, x0, normA, maxit, myeps)
% implementace algoritmu MPRGP - minimalizacni algoritmus
% s linearnim omezenim
%
% VSTUP:
% A - symetricka positivne definitni matice
% b - vektor pravych stran
% l - vektor linearnich omezeni x(i) >= l(i)
% u - upper bound x(i) <= u(i)
% eps - pozadovana presnost
%
% VYSTUP:
% x - reseni
% it - pocet hlavnich iteraci
%

% my_opt.gamma = 1
% my_opt.normA = norm(A)
% my_opt.my_eps
% my_opt.max_it

hess_mult = 0;

% 0. Initialization
gamma = 1; % magic parameter
alpha_bar = 1.95/normA; % 0 < alpha_bar <= 2*norm(inv(A));
x = get_projection_box(x0,lb,ub); % x \in \Omega
it = 0;

g = A*x - b; hess_mult = hess_mult + 1;
fi = get_fi(x,g,lb,ub);
beta = get_beta(x,g,lb,ub);
gp = fi + beta;

p = fi;

gp_norms(1) = norm(gp);
fi_norms(1) = norm(fi);
beta_norms(1) = norm(beta);

while my_stopping_criterion(x,norm(gp),myeps) && it < maxit
    if dot(beta,beta) <= gamma^2*dot(fi,fi)
        % 1. Proportional x_k. Trial conjugate gradient step
        Ap = A*p; hess_mult = hess_mult + 1;
        
        alpha_cg = (g'*p)/(p'*Ap);
        y = x - alpha_cg*p;
        alpha_f = get_alpha_f(x,p,lb,ub);
        
        if min([alpha_cg, alpha_f]) == inf
            disp('MPRGP: No solution.')
        end
        
        if alpha_cg <= alpha_f
            % 2. Conjugate gradient step
            x = y;
            g = g - alpha_cg*Ap;
            fi = get_fi(x,g,lb,ub);
            beta = fi'*Ap/(p'*Ap);
            p = fi - beta*p;
        else
            % 3. Expansion step
            x = x-alpha_f*p;
            g = g - alpha_f*Ap;
            fi = get_fi(x,g,lb,ub);
            x = get_projection_box(x-alpha_bar*fi,lb,ub);
            g = A*x-b; hess_mult = hess_mult + 1;
            p = get_fi(x,g,lb,ub);
        end
    else
        % 4. Proportioning step
        Abeta = A*beta; hess_mult = hess_mult + 1;
        
        alpha_cg = g'*beta/(beta'*Abeta);
        x = x - alpha_cg*beta;
        g = g - alpha_cg*Abeta;
        p = get_fi(x,g,lb,ub);
    end
    fi = get_fi(x,g,lb,ub);
    beta = get_beta(x,g,lb,ub);
    gp = fi+beta;
    
    gp_norms(it+1) = norm(gp);
    fi_norms(it+1) = norm(fi);
    beta_norms(it+1) = norm(beta);
    
    it = it + 1;
    
end

end

function fi = get_fi(x,g,l,u)
% box constraints
fi = g;

if ~isempty(l)
    fi(x <= l) = 0;
end
if ~isempty(u)
    fi(x >= u) = 0;
end

end

function beta = get_beta(x,g,l,u)
beta = zeros(size(x));

if ~isempty(l)
    i = find(x == l);
    beta(i) = min(g(i),0);
end

if ~isempty(u)
    j = find(x == u);
    beta(j) = max(g(j),0);
end
end

function alpha_f = get_alpha_f(x,p,l,u)

alpha_f1 = inf;
if ~isempty(l)
    i = find(p > 0);
    if min(size(i)) > 0
        alpha_f1 = min((x(i)-l(i))./p(i));
    end
end

alpha_f2 = inf;
if ~isempty(u)
    j = find(p < 0);
    if min(size(j)) > 0
        alpha_f2 = min((x(j)-u(j))./p(j));
    end
end

alpha_f = min(alpha_f1,alpha_f2);
end

function [result] = my_stopping_criterion(x,norm_gp,myeps)
result = (norm_gp > myeps);
end



