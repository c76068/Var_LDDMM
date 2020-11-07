% function [f,g] = obj_vfnorm(X,y,defo,objfun)
function [f,g] = obj_vfnorm(y)
  global Z objfun_1 
   
   if isstruct(Z) == 0
       [~,d] = size(Z);
       m = length(y)/d;
   else
       [~,d] = size(Z.center);
       m = length(y)/((length(Z.vector)+1)*d);     
   end
   
    Y = vec2stru(y,m,d);
    f = varifoldnorm(Z,Y,objfun_1);
%    f = Two_vec_vfnorm(Z,Y,objfun_1);
%    f = Two_vec_vfnorm(X,Y,defo,objfun);
   
   if nargout > 1 % gradient required
      G = dvarifoldnorm(Y,Z,objfun_1);
      g = stru2vec(G);
%     G = dTwo_vec_vfnorm(Y,Z,objfun_1);
%     G = dTwo_vec_vfnorm(Y,X,defo,objfun);
      
   end
end

