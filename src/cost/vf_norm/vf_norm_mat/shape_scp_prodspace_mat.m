function res=shape_scp_prodspace_mat(X,Y,objfun)

    if isfield(X,'vector') == 0
       m = 0;     
    else         
       m = length(X.vector);    
    end
    
    switch m 
        case 0
            [~,res] = Zero_vec_vfnorm_mat(X,Y,objfun);
        case 1
            [~,res] = One_vec_vfnorm_mat(X,Y,objfun);
        case 2
            [~,res] = Two_vec_vfnorm_mat(X,Y,objfun);
    end
    
end

