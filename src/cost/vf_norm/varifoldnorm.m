function g = varifoldnorm(X,Y,objfun)
    
   switch objfun.method
    case 'keops'
        g = vfnorm_keops(X,Y,objfun);
    case 'matlab'
        g = vfnorm_mat(X,Y,objfun);
   end
end

function g = vfnorm_mat(X,Y,objfun)
      
    if isfield(X,'vector') == 0
       m = 0;     
    else         
       m = length(X.vector);    
    end

     switch m
         case 0
             g = Zero_vec_vfnorm_mat(X,Y,objfun);
         case 1
             g = One_vec_vfnorm_mat(X,Y,objfun);
         case 2
             g = Two_vec_vfnorm_mat(X,Y,objfun);
     end 
 
end

function g = vfnorm_keops(X,Y,objfun)

 PXX = shape_scp_prodspace_keops(X,X,objfun);
 PYY = shape_scp_prodspace_keops(Y,Y,objfun);
 PXY = shape_scp_prodspace_keops(X,Y,objfun);
 g= PXX + PYY - 2* PXY;

end