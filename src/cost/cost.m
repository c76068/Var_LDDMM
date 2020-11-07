function  [ept,Ep,H_p] = cost(X,P,Target,objfun,defo)
% This function  computes the cost.
%
% Inputs :
%   X.center: is a (n x d) matrix containing the points.
%   X.vector: a cell contains first  and second vectors of the 2-vector.
%   mom.center: is a (n x d) matrix containing the spatial momentums.
%   mom.vector: a cell contains momentums about first  and second vectors of the 2-vector.
%   defo: is a structure of deformations.


% Target: similar to X
% K_q: Positive definite matrix used in computing the Hamiltonian energy.
% Cost can be computed without inputuing K_q

%Output:
%ept.X,ept.P: evolution of states and momemtum variables
%ept.ham: Hamiltonian energy
%ept.dat: data attachment term
%ept.cost: Total cost


   [n,d]=size(X.center); 
   
   [ept.X,ept.P]=forward(X,P,defo,1);
   [ept.ham,H_p] = Ham(X,P,defo);
   
   Z = ept.X{end};
   if isfield(X,'weight') ~= 0
      Z.weight = X.weight;
   end
   
    ept.dat = varifoldnorm(Z,Target,objfun);
    ept.cost = ept.ham + objfun.lambda*ept.dat;
    
    if nargout >1
       
        
       dG_1 = dvarifoldnorm(Z,Target,objfun);
       dpfinal.center = zeros(n,d);
       
       if isfield(X,'vector') ~= 0
          for i=1:length(X.vector) 
             dpfinal.vector{i} = zeros(n,d);
          end
       end
       
       [~,dpinit]=backward_tan(dG_1,dpfinal,ept,defo);
       
       Ep = add_XY_h(H_p,dpinit,objfun.lambda);
    end

end


