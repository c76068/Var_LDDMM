function [dxinit,dpinit]=backward_tan(dxfinal,dpfinal,ept,defo)

% [dxinit,dpinit]=BACKWARD(dxfinal,dpfinal,shooting,defo) performs the
% backward integration for the tangential Hamiltonian flow.
%
% Input :
%   dxfinal.center : gradient in (point) wrt the final position
%   dxfinal.vector : gradient in (direction vectors) wrt the final directions
%   dpfinal : gradient in  (momenta) wrt the final momenta
%   containing dpfinal.cneter and dpfinal.vector
%   ept : a  structure containing evolution path of the states and momentum
%
% Ouput
%   dxinit: initials of the covaraiables of the states

%   dpinit: initials of the covaraiables of the momemtum 
%
% See also : ddHrtP_tan, forward_tan, dHr_tan
%

 method = defo.odemethod;
 
if strcmp(method,'middle_point')
%     [x_evol,d_evol,p1_evol,p2_evol]=bmpt(x_init,d_init,p1_init,p2_init,defo,tf);
    [dxinit,dpinit]=bmpt(dxfinal,dpfinal,ept,defo);
end

 if strcmp(method,'rk4')
     [dxinit,dpinit]=brk4(dxfinal,dpfinal,ept,defo);
 end
end

function [dxinit,dpinit]=bmpt(dxfinal,dpfinal,ept,defo)

% middle point method
% Step size
h=-1/defo.nb_euler_steps;

% Initiatialze dxinit (output) to dxfinal before backward integration
dxinit = dxfinal;
dpinit = dpfinal;

% pause()
 for i=defo.nb_euler_steps+1:-1:2
   % backward midpoint method :    
   
   [nPx{1},nPp{1}] = fddh2(ept.X{i},ept.P{i},dxinit,dpinit,defo,h/2);
   [nPx{2},nPp{2}] = fddh2(add_XY_h(ept.X{i-1},ept.X{i},[0.5 0.5]),add_XY_h(ept.P{i-1},ept.P{i},[0.5 0.5]),nPx{1},nPp{1},defo,h);

   dxinit = add_XY_h(dxinit,add_XY_h(nPx{2},nPx{1},-1),1);
   dpinit = add_XY_h(dpinit,add_XY_h(nPp{2},nPp{1},-1),1);
   
 end
 
end


function [dxinit,dpinit]=brk4(dxfinal,dpfinal,ept,defo)
% rk4 scheme
% Step size
h=-1/defo.nb_euler_steps;

% Initiatialze dxinit (output) to dxfinal before backward integration
dxinit= dxfinal;
dpinit = dpfinal;

% pause()
for i=defo.nb_euler_steps+1:-1:2
   %  method :    
    [nPx{1},nPp{1}] = fddh2(ept.X{i},ept.P{i},dxinit,dpinit,defo,h);
    [nPx{2},nPp{2}] = fddh2(add_XY_h(ept.X{i-1},ept.X{i},[0.5 0.5]),add_XY_h(ept.P{i-1},ept.P{i},[0.5 0.5]),add_XY_h(dxinit,nPx{1},1/2),add_XY_h(dpinit,nPp{1},1/2),defo,h);
    [nPx{3},nPp{3}] = fddh2(add_XY_h(ept.X{i-1},ept.X{i},[0.5 0.5]),add_XY_h(ept.P{i-1},ept.P{i},[0.5 0.5]),add_XY_h(dxinit,nPx{2},1/2),add_XY_h(dpinit,nPp{2},1/2),defo,h);
    [nPx{4},nPp{4}] = fddh2(ept.X{i-1},ept.P{i-1},add_XY_h(dxinit,nPx{3},1),add_XY_h(dpinit,nPp{3},1),defo,h);
    
     w = [1/6 1/3 1/3 1/6];
    for ind =1:4
       dxinit = add_XY_h(dxinit,nPx{ind},w(ind));
       dpinit = add_XY_h(dpinit,nPp{ind},w(ind));
    end
    
    %pause()
end
end

function [nPx,nPp]=fddh2(X,P,Px,Pp,defo,h)
% Basic Euler step 
   [dPq, dPp] = ddHrtP_tan(X,P,Px,Pp,defo);
%    [dPx,dPv,dPp1,dPp2] = ddHrtP_tan(x,v,p1,p2,Px,Pv,Pp1,Pp2,defo);
   
   
   if  strcmp(defo.odemethod,'middle_point')
    
     nPx = add_XY_h(Px,dPq,h);
     nPp = add_XY_h(Pp,dPp,h);
          
   end
   
   if strcmp(defo.odemethod,'rk4')
    
     nPx = add_XY_h(X,dPq,[0,h]);
     nPp = add_XY_h(P,dPp,[0,h]);  
    
   end
end