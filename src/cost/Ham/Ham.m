function  [ham,ham_grad] = Ham(X,P,defo)


   switch defo.method
    case 'keops'
        HR = @Ham_keops;
    case 'matlab'
        HR = @Ham_mat;
   end
   
   [ham,ham_grad] = HR(X,P,defo);

end