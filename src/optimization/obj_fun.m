function [e, g]= obj_fun(X,p,Y,objfun,defo,K_q)
        [n,d] = size(X.center);
        P.center = reshape(p(1:n*d),d,n)';
        P.vector{1} = reshape(p(n*d+1:2*n*d),d,n)';
        P.vector{2} = reshape(p(2*n*d+1:3*n*d),d,n)';
        ept = cost(X,P,Y,objfun,defo,K_q);
        e = ept.cost;
        G = gsgrad(X,P,Y,objfun,defo,K_q);
        g = [reshape(G.center',n*d,1); reshape(G.vector{1}',n*d,1);...
    reshape(G.vector{2}',n*d,1)];
end