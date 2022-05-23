function [x, it, gp, hess_mult, gp_norms, fi_norms, beta_norms] = mprgp(A, b, lb, x0, normA, opts)
% implementace algoritmu MPRGP - minimalizacni algoritmus
% s linearnim omezenim
%
% VSTUP:
% A - symetricka positivne definitni matice
% b - vektor pravych stran
% l - vektor linearnich omezeni x(i) >= l(i)
% eps - pozadovana presnost
%
% VYSTUP:
% x - reseni
% it - pocet hlavnich iteraci
%

hess_mult = 0;

% 0. Initialization
gamma = opts.gamma; % magic parameter
alpha_bar = 1.95/normA; % 0 < alpha_bar <= 2*norm(inv(A));
x = get_projection_lowerbound(x0,lb); % x \in \Omega
it = 0;

g = A*x - b; hess_mult = hess_mult + 1;
fi = get_fi(x,g,lb);
beta = get_beta(x,g,lb);
gp = fi + beta;

p = fi;

gp_norms(1) = norm(gp);
fi_norms(1) = norm(fi);
beta_norms(1) = norm(beta);

while my_stopping_criterion(x,norm(gp),opts) && it < opts.maxit
    if dot(beta,beta) <= gamma^2*dot(fi,fi)
        % 1. Proportional x_k. Trial conjugate gradient step
        Ap = A*p; hess_mult = hess_mult + 1;
        
        alpha_cg = (g'*p)/(p'*Ap);
        y = x - alpha_cg*p;
        alpha_f = get_alpha_f(x,p,lb);
        
        if min([alpha_cg, alpha_f]) == inf
            disp('MPRGP: No solution.')
        end
        
        if alpha_cg <= alpha_f
            % 2. Conjugate gradient step
            x = y;
            g = g - alpha_cg*Ap;
            fi = get_fi(x,g,lb);
            beta = fi'*Ap/(p'*Ap);
            p = fi - beta*p;
        else
            % 3. Expansion step
            x = x-alpha_f*p;
            g = g - alpha_f*Ap;
            fi = get_fi(x,g,lb);
            x = get_projection_lowerbound(x-alpha_bar*fi,lb);
            g = A*x-b; hess_mult = hess_mult + 1;
            p = get_fi(x,g,lb);
        end
    else
        % 4. Proportioning step
        Abeta = A*beta; hess_mult = hess_mult + 1;
        
        alpha_cg = g'*beta/(beta'*Abeta);
        x = x - alpha_cg*beta;
        g = g - alpha_cg*Abeta;
        p = get_fi(x,g,lb);
    end
    fi = get_fi(x,g,lb);
    beta = get_beta(x,g,lb);
    gp = fi+beta;
    
    gp_norms(it+1) = norm(gp);
    fi_norms(it+1) = norm(fi);
    beta_norms(it+1) = norm(beta);
    
    it = it + 1;
    
end

end

function fi = get_fi(x,g,lb)
% box constraints
fi = g;

if ~isempty(lb)
    fi(x <= lb) = 0;
end

end

function beta = get_beta(x,g,lb)
beta = zeros(size(x));

if ~isempty(lb)
    i = find(x == lb);
    beta(i) = min(g(i),0);
end

end

function alpha_f = get_alpha_f(x,p,lb)

alpha_f = inf;
if ~isempty(lb)
    i = find(p > 0);
    if min(size(i)) > 0
        alpha_f = min((x(i)-lb(i))./p(i));
    end
end

end

function [result] = my_stopping_criterion(x,norm_gp,opts)
    if opts.smalsem
        result = (norm_gp > max([opts.myeps, min([norm(opts.smalsem_obj.BE*x),opts.smalsem_obj.eta] )]) );
    else
        result = (norm_gp > opts.myeps);
    end
end



