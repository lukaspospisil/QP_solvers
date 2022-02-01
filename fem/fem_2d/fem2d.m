function [A,b] = fem1d(P,E,GammaD,GammaN,uD,gN)

%       [A,b] = MKP(P,E,GammaD,GammaN,tau,f,uD,gN)
%       Sestavi soustavu metodou konecnych prvku pro okrajovou ulohu
%       -div[tau(x)*grad(u'(x))] = f(x) v Omega
%                           u(x) = uD(x) na GammaD
%                tau(x)*du(x)/dn = gN(x) na GammaN
%       + podm. prechodu
%       P ... souradnice uzlu (matice 2 x n)
%       E ... indexy uzlu trojuhelnikove site (matice 3 x m)
%       GammaD ... indexy Dirichletovych uzlu (vektor nD)
%       GammaN ... dvojice indexu Neumannovych uzlu (matice 2 x mN)
%       tau,f ... po trojuhelnicich konst. (vektory m)
%       uD ... po Dir. useckach afinni (vektor nD)
%       gN ... po Neum. useckach konst. (vektor mN)


n = size(P,2);
A = sparse(n,n);
b = zeros(n,1);

m = size(E,2);
nD = length(GammaD);
mN = size(GammaN,2);

% vypocet objemovych integralu
for k=1:m
    ek = E(:,k);
    Pk = P(:,ek);
    
    tauk = 1;
    %vypocet lokalni matice tuhosti
    Ak = fem2d_alok_gradgrad(Pk,tauk);

    fk = -10; 
    %prava strana, akorat objem jehlanu
    bk = fem2d_blok_id(Pk,fk);
    
    A(ek,ek) = A(ek,ek) + Ak;
    b(ek) = b(ek) + bk;
end

for i=1:length(GammaD)
    A(:,GammaD(i)) = zeros(size(A,1),1);
    A(GammaD(i),:) = zeros(1,size(A,1));
    A(GammaD(i),GammaD(i)) = 1;
end

% symetrize matrix.. it seems that there are some numerical errors
A = (A+A')/2;

end
