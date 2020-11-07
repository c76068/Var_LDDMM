function [X_evol,P_evol]=forward(X_init,P_init,defo,tf)
% [X_evol,P_evol]=FORWARD(X_init,P_init,defo,final_time) compute
% Forward integration of the Hamiltonian flow from initial coupled
% configuration of varifold/momentums.
%
% Input :
%  X_init.center : initial coordinates of the points (position) in a (n x d) matrix.
%  X_init.vector : a cell contains first  and second vectors of the initial 2-vector.
%  P_init.center : is a (n x d) matrix containing the initial spatial momentums.
%  P_init.vector : a cell contains initial momentums about first  and second vectors of the 2-vector..
%  defo : structure containing the parameters of deformations (kernel_size_mom,method,nstep,...)
%  final_time : final time (optional, and fixed by default to 1)
%
% Output
%  X_evol.center : a cell list containing evolution path of X.center ( points_evol{i} is a n x d matrix and i ranges from 1 to defo_options.nstep+1)
%  X_evol.vector{i} : a cell containing evolution of X.vector{i}
%  P_evol.center : a cell list containing evolution path of spatial momentums ( nomentums_evol{i} is a n x d matrix and i ranges from 1 to defo_options.nstep+1)
%  P_evol.vector{i} : a cell containing evolution of momentums aboutX.vector{i}

   method = defo.odemethod;
%  method = 'rk4';
   if nargin == 3
    tf=1;
   end

% if strcmp(defo.action,'normalized')
%     weight = sqrt(Inn(X_init.vector{1},X_init.vector{1}));
%     X_init.vector{1} = X_init.vector{1}./weight;
% end
%  

   
   
if strcmp(method,'middle_point')
    [X_evol,P_evol]=mpt(X_init,P_init,defo,tf);
end

 if strcmp(method,'rk4')
     [X_evol,P_evol]=rk4(X_init,P_init,defo,tf);
 end
 
%  if strcmp(defo.action,'normalized')
%     for t=1:length(X_evol)
%        X_evol{t}.vector{1} =  weight.*X_evol{t}.vector{1};    
%     end
%  end

end


function [X_evol,P_evol]=mpt(X_init,P_init,defo,tf)


[~,d] =size(X_init.center);
if nargin == 3
    tf=1;
end

X_evol=cell(1,defo.nb_euler_steps+1);
P_evol=cell(1,defo.nb_euler_steps+1);


dt=tf/defo.nb_euler_steps;

X_evol{1} = X_init;
P_evol{1} = P_init;

if strcmp(defo.action,'normalized')
    weight = sqrt(dot(X_init.vector{1},X_init.vector{1},2));
end

for i=1:defo.nb_euler_steps
    
    
    [nX1,nP1]=fdh(X_evol{i},P_evol{i},defo,dt/2);
    [nX2,nP2]=fdh(nX1,nP1,defo,dt);
    
    X_evol{i+1} = add_XY_h(X_evol{i},add_XY_h(nX2,nX1,[1,-1]),1);
    P_evol{i+1} = add_XY_h(P_evol{i},add_XY_h(nP2,nP1,[1,-1]),1);
    
    if strcmp(defo.action,'normalized')
       X_evol{i+1} = direc_normalize(X_evol{i+1},weight); %normalize the length of direction vector
    end
end
end

function [X_evol,P_evol]=rk4(X_init,P_init,defo,tf)

[~,d] =size(X_init.center);
if nargin == 3
    tf=1;
end

X_evol=cell(1,defo.nb_euler_steps+1);
P_evol=cell(1,defo.nb_euler_steps+1);


dt=tf/defo.nb_euler_steps;

X_evol{1} = X_init;
P_evol{1} = P_init;

if strcmp(defo.action,'normalized')
    weight = sqrt(dot(X_init.vector{1},X_init.vector{1},2));
end

for i=1:defo.nb_euler_steps
    
       
    tempX = X_evol{i}; tempP =P_evol{i};
    X_evol{i+1} = X_evol{i}; P_evol{i+1} = P_evol{i};
    nX=cell(1,4);
    nP=cell(1,4);
    
    for ind=1:4
       
    [nX{ind},nP{ind}]=fdh(tempX,tempP,defo,dt); 
    
      h = [1 1/2 1/2 1];
     
     if ind<4      
      tempX = add_XY_h(X_evol{i},nX{ind},h(ind+1)); tempP = add_XY_h(P_evol{i},nP{ind},h(ind+1));
     end
     
      X_evol{i+1} = add_XY_h(X_evol{i+1},nX{ind},1/(6*h(ind)));
      P_evol{i+1} = add_XY_h(P_evol{i+1},nP{ind},1/(6*h(ind)));
      
       if strcmp(defo.action,'normalized')
          X_evol{i+1} = direc_normalize(X_evol{i+1},weight); %normalize the length of direction vector
       end
    
    end
    
end
end

function [nX,nP]=fdh(X,P,defo,h)
% This function implements an elementary Euler Step
%
% Inputs :
% Input :
%  X.center : points (position) in a (n x d) matrix.
%  X.vector : a cell contains first  and second vectors of the 2-vector.
%  P.center : is a (n x d) matrix containing the spatial momentums.
%  P.vector : a cell contains momentums about first  and second vectors of the 2-vector..
%  defo: is a structure of deformations.
%  h is the time step.

% Outputs
%  nX.center : (n x d) matrix containing the new points.
%  nX.vector : a cell contains first  and second vectors of the new 2-vector.
%  nP.center : (n x d) matrix containing the new the momentums.
%  nP.vector : a cell contains new momentums about first  and second vectors of the 2-vector
 

  [dqHr,dpHr] = dHr_tan(X,P,defo);
  
 if  strcmp(defo.odemethod,'middle_point')
   
     nX = add_XY_h(X,dpHr,h);
     nP = add_XY_h(P,dqHr,-h);
     
 end
 
 if  strcmp(defo.odemethod,'rk4')
     
     nX = add_XY_h(X,dpHr,[0,h]);
     nP = add_XY_h(P,dqHr,[0,-h]);

 end

    
end


function X = direc_normalize(X,r)
 [~,d] = size(X.center);
 m = length(X.vector);
 
 if m==1
    norm_V = sqrt(dot(X.vector{1},X.vector{1},2));
    X.vector{1} = repmat(r,1,d).*X.vector{1}./repmat(norm_V,1,d);
 elseif m==2
     
 end

end
