
function [dPq, dPp] = ddHrtP_tan(X,P,Px,Pp,defo)
% [dPx,dPp] = DDhRTP(x,p,Px,Pp,defo) computes the adjoints Hamiltonian 
% system usingthe methods given in defo.method.
%
% Inputs :
% Input :
%  X.center : points (position) in a (n x d) matrix.
%  X.vector : a cell contains first  and second vectors of the initial 2-vector.
%  P.center : is a (n x d) matrix containing the initial spatial momentums.
%  P.vector : a cell contains initial momentums about first  and second vectors of the 2-vector..
%  Px.center : is a (n x d) matrix (adjoint variable of X.center)
%  Px.vector : a cell containing adjoint variables of X.vector
%  Pp.center : is a (n x d) matrix (adjoint variable of P.center) 
%  Pp.vector : a cell containing adjoint variables of P.vector
%defo : structure containing the parameters of deformations (kernel_size_mom,method,nstep,...)

%
% Outputs
%   dPq.center : (n x d) matrix containing the update of Px.center.
%   dPq.vector : a cell containing the update of Px.vector.
%   dPp.center : (n x d) matrix containing the update of Pp.center.
%   dPp.vector : a cell containing the update of Pp.vector.

% finite diff perturbation, use double to increase precision

%hh = 1e-7;
  switch defo.method
     case 'keops'
        hh = 1e-3;
     case 'matlab'
        hh = 1e-7;
  end

defo.precision = 'double';
% diff=@(ph,mh,hh) (ph-mh) / (2*hh);

     
     [dqH_pm, dpH_pm] = dHr_tan(add_XY_h(X,Pp,hh),add_XY_h(P,Px,-hh),defo);
     [dqH_mp, dpH_mp] = dHr_tan(add_XY_h(X,Pp,-hh),add_XY_h(P,Px,hh),defo);
     
     dPq = add_XY_h(dqH_pm,dqH_mp,[1 -1]/(2*hh));
     dPp = add_XY_h(dpH_pm,dpH_mp,[1 -1]/(2*hh));
%    [dxxmhP,dxvmhP,dxp1mhP,dxp2mhP] =  dHr_tan(x-hh*Pp1,v-hh*Pp2,p1+hh*Px,p2+hh*Pv,defo);
%    [dxxphP,dxvphP,dxp1phP,dxp2phP] =  dHr_tan(x+hh*Pp1,v+hh*Pp2,p1-hh*Px,p2-hh*Pv,defo);

%    dPx = diff(dxxphP,dxxmhP,hh);
%    dPv = diff(dxvphP,dxvmhP,hh);
%    dPp1 = diff(dxp1phP,dxp1mhP,hh);
%    dPp2 = diff(dxp2phP,dxp2mhP,hh);


end

