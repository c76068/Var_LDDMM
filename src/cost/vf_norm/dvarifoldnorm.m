function g=dvarifoldnorm(X,Y,objfun)
   
   switch objfun.method
    case 'keops'
        g = dvfnorm_keops(X,Y,objfun);
    case 'matlab'
        g = dvfnorm_mat(X,Y,objfun);
   end

end

function g = dvfnorm_mat(X,Y,objfun)
    
    if isfield(X,'vector') == 0
       m = 0;     
    else         
       m = length(X.vector);    
    end

     switch m
         case 0
             g = dZero_vec_vfnorm_mat(X,Y,objfun);
         case 1
             g = dOne_vec_vfnorm_mat(X,Y,objfun);
         case 2
             g = dTwo_vec_vfnorm_mat(X,Y,objfun);
     end
 
end

function g = dvfnorm_keops(X,Y,objfun)
 
 GXX = dshape_scp_prodspace_keops(X,X,objfun);
 GXY = dshape_scp_prodspace_keops(X,Y,objfun);
 g = add_XY_h(GXX,GXY,[2,-2]);
 
end