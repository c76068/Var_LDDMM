function [f,g] = obj_cost(p)
% Simple wrapper function for cost that can be called by Hanso bfgs.
  global X_1 Y_1 objfun_1 defo_1 
   
  % m: number of diracs
%    if isstruct(X_1) == 0
%        [~,d] = size(X_1);
%        m = length(p)/d;
%    else
%        [~,d] = size(X_1.center);
%        m = length(p)/((length(X_1.vector)+1)*d);     
%    end
    
    [m,d] = size(X_1.center);
   
    P = vec2stru(p,m,d);
    
%     ept = cost(X_1,P,Y_1,objfun_1,defo_1);
%     f = ept.cost;
%    
%    if nargout > 1 % gradient required
%        
%       G = gsgrad(X_1,P,Y_1,objfun_1,defo_1);
%       g = stru2vec(G);     
%    end

    if nargout > 1 % gradient required
        [ept,G] = cost(X_1,P,Y_1,objfun_1,defo_1);
        g = stru2vec(G);  
    else
        ept = cost(X_1,P,Y_1,objfun_1,defo_1);
    end
    
     f = ept.cost;
end
