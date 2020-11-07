 function K_q = Kqop(X,defo)
    [n,d]=size(X.center);
    if isfield(X,'vector') == 0
       m = 0;     
    else         
       m = length(X.vector);    
    end
%     m = length(X.vector);
   
    switch m
       case 0
          x = X.center;
          S=zeros(n);
          for l=1:d
             S=S+(repmat(x(:,l),1,n)-repmat(x(:,l)',n,1)).^2;
          end
          I = eye(d);
          Gamma_0 = rho(S,0,defo.kernel_size_mom); 
          K_q = kron(Gamma_0,I);
          
       case 1
   
          x = X.center; v =X.vector{1}; 
          % Calcul de A=exp(-|x_i -x_j|^2/(lam^2))
          S=zeros(n); % S = (|x_i-x_j|^2)_{i,j}
          for l=1:d
             S=S+(repmat(x(:,l),1,n)-repmat(x(:,l)',n,1)).^2;
          end
          Gamma_0 = rho(S,0,defo.kernel_size_mom);
          Gamma_1 = rho(S,1,defo.kernel_size_mom);
          Gamma_2 = rho(S,2,defo.kernel_size_mom);
    
          %Inner product matrices
          Cxv = x*v'; %Inner products of x and v
          Cv = v*v';
    
          A = repmat(diag(Cxv),1,n)'-Cxv;
          H = 2*Gamma_1.*A;
          J =4*Gamma_2.*(A').*A-2*Gamma_1.*Cv;

          I = eye(d);
          K_q = [kron(Gamma_0,I) kron(H,I);kron(H,I)' kron(J,I)];
    
          if strcmp(defo.action,'normalized')
             v_nor = v./sqrt(Inn(v,v)); %normalize v
             PJ = zeros(n*d);
          
             for i=1:n
                PJ((i-1)*d+1:i*d,(i-1)*d+1:i*d) = eye(d)-v_nor(i,:)'*v_nor(i,:);
             end
          
            Coc = zeros(2*n*d);
            Coc(1:n*d,1:n*d) = eye(n*d);
            Coc(n*d+1:2*n*d,n*d+1:2*n*d) = PJ;
            K_q = Coc'*K_q*Coc;
          end
    
       case 2
   
       
       % Calcul de A=exp(-|x_i -x_j|^2/(lam^2))
       S=zeros(n); % S = (|x_i-x_j|^2)_{i,j}
       for l=1:d
           S=S+(repmat(X.center(:,l),1,n)-repmat(X.center(:,l)',n,1)).^2;
       end
       Gamma_0 = rho(S,0,defo.kernel_size_mom);
       Gamma_1 = rho(S,1,defo.kernel_size_mom);
       Gamma_2 = rho(S,2,defo.kernel_size_mom);
    
       %Inner product matrices
       scp_x1 = X.center*X.vector{1}'; % <x_i,u^1_k>
       scp_x2 = X.center*X.vector{2}'; % <x_i,u^2_k>
       scp_11 = X.vector{1}*X.vector{1}'; % <u_i,u_k>
       scp_22 = X.vector{2}*X.vector{2}';
       scp_12 = X.vector{1}*X.vector{2}';
    
       C_1 = repmat(diag(scp_x1),1,n)'-scp_x1; 
       C_2 = repmat(diag(scp_x2),1,n)'-scp_x2;
       D_1 = -C_1'; D_2 = -C_2';
    
       A_1 = 2*Gamma_1.*C_1;  A_2 = 2*Gamma_1.*C_2;
       B_11 = -(4*Gamma_2.*D_1.*C_1 + 2*Gamma_1.*scp_11);
       B_12 = -(4*Gamma_2.*D_1.*C_2 + 2*Gamma_1.*scp_12);
       B_22 = -(4*Gamma_2.*D_2.*C_2 + 2*Gamma_1.*scp_22);
    
       I = eye(d);
       K_q = [kron(Gamma_0,I) kron(A_1,I) kron(A_2,I);kron(A_1,I)' kron(B_11,I) kron(B_12,I);...
             kron(A_2,I)' kron(B_12,I)' kron(B_22,I)];
    

    end
    
 end