function  [ham,ham_grad] = Ham_mat(X,P,defo)

   
   [n,d]=size(X.center);
%    m = length(X.vector);

   if isfield(defo,'Kq') == 0
      K_q = Kqop(X,defo); 
   else   
      K_q = defo.Kq;
   end
      
    p = stru2vec(P);
    ham_grad = K_q*p;
    ham= 0.5*(p'*ham_grad);
    ham_grad = vec2stru(ham_grad,n,d);

end

