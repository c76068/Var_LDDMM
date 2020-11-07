function [Y,summary] = proj2_M_dirac(X,Y,objfun,options)
% Find the closest linear combination of M diracs to \mu_X, the linear combination of
% diracs that represents X
%Input:
% X: Original surface, include X.center, X.vector
% Y: initial linear combination of M diracs
%Output:
% Y: closest linear combination of M diracs to \mu_X
  global Z objfun_1
  Z = X; objfun_1 = objfun;
  [n,d] = size(Y.center);

  y_0 = stru2vec(Y);
%   fun =@(a)obj_vfnorm(X,a,defo,objfun);
%    options = optimoptions(@fminunc,'Display','iter','Algorithm','quasi-newton');
 
if nargin<4
%  options = optimoptions('fminunc','Algorithm','quasi-newton','SpecifyObjectiveGradient',true,'Display','iter','MaxIterations',10000);
%    options.niter = 200;
   options.record_history =true;
   options.nvec = 10;
   options.prtlevel = 2;
end
%  options = optimoptions(@fminunc,'Display','iter','MaxFunctionEvaluations',1000*length(A_0),'MaxIterations',1000);
%   y = fminunc(fun,y_0,options);
  [y,summary] = perform_bfgs('obj_vfnorm', y_0, options);

  Y = vec2stru(y,n,d);
end