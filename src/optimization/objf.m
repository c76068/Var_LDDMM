function feval = objf(X,y,objfun)

   [~,d] = size(X.center);
   Y.center = y(1:d,:)';
   Y.vector{1} = y(d+1:2*d)';
   Y.vector{2} = y(2*d+1:3*d)';
   feval=Two_vec_vfnorm(X,Y,objfun);
end