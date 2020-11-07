function [dqHr,dpHr] = dHr_tan(X,P,defo)
% This function  computes the reduced Hamiltonian system.
%
% Inputs :
%   X.center: is a (n x d) matrix containing the points.
%   X.vector: a cell contains first  and second vectors of the 2-vector.
%   P.center: is a (n x d) matrix containing the spatial momentums.
%   P.vector: a cell contains momentums about first  and second vectors of the 2-vector.
%   defo: is a structure of deformations.


% Outputs
%   dHr.center : (n x d) matrix containing the gradient of Hr wrt to x at point (X,P).
%   dHr.vector : a cell containing the gradient of Hr wrt to X.vector at point (X,P).
%   dHr.pcenter: (n x d) matrix containing the gradient of Hr wrt to P.center at point (X,P).
%   dHr.pvector :a cell containing the gradient of Hr wrt to P.vector.vector at point (X,P).

% if strcmp(defo.action,'pushforward')
%  [dqHr,dpHr] = dHr_tanpush(X,P,defo);   
% end
switch defo.action
    case 'pushforward'
        [dqHr,dpHr] = dHr_tanpush(X,P,defo);
    case 'normalized'      
        [dqHr,dpHr] = dHr_tannormal(X,P,defo);
end



% if strcmp(defo.action,'normalized')
%  dHr = dHr_tannormal(X,P,defo); 
% end

end

function [dqHr,dpHr] = dHr_tanpush(X,P,defo)
 
 switch defo.method
    case 'keops'
        DHR = @dHr_tan_push_keops;
    case 'matlab'
        DHR = @dHr_tan_push_mat;
end

[dqHr,dpHr] = DHR(X,P,defo);
end

function [dqHr,dpHr] = dHr_tannormal(X,P,defo)
   
   switch defo.method
     case 'keops'
        DHR = @dHr_tan_normal_keops;
     case 'matlab'
        DHR = @dHr_tan_normal_mat;
   end

[dqHr,dpHr] = DHR(X,P,defo);
end

function [dqHr,dpHr] = dHr_tan_push_keops(X,P,defo)
     %here f = dHr
%     [n,d]=size(X.center);
%      m=length(X.vector);
     
    [~,d]=size(X.center);
   
    if isfield(X,'vector') == 0
       m = 0;     
    else         
       m = length(X.vector);    
    end
    
  switch m
      case 0
         
         V = Kernel('Exp(-p*SqNorm2(x-y))*p1',['x=Vx(' num2str(d) ')'],['y=Vy(' num2str(d) ')'],['p1=Vy(' num2str(d) ')'],'p=Pm(1)');
         dV_ad = GradKernel(V,'x',['h1=Vx(' num2str(d) ')']);
         
         x = X.center'; p = P.center'; 
         oos2 = 1/defo.kernel_size_mom^2;
         
         dpHr.center = V(x,x,p,oos2)';
         dqHr.center = dV_ad(x,x,p,oos2,p)';
         
      case 1
    
         V = Kernel('Exp(-p*SqNorm2(x-y))*p1+Grad(Grad(Exp(-p*SqNorm2(x-y))*p2,y,h),h,u)',...
         ['x=Vx(' num2str(d) ')'],['y=Vy(' num2str(d) ')'],['p1=Vy(' num2str(d) ')'],...
         ['p2=Vy(' num2str(d) ')'],['u=Vy(' num2str(d) ')'],'p=Pm(1)',['h=Vy(' num2str(d) ')']);

   
         dV_ad = GradKernel(V,'x',['h1=Vx(' num2str(d) ')']);
         dV = GradKernel(dV_ad,'h1',['h2=Vx(' num2str(d) ')']);
         dV_2ad = GradKernel(dV_ad,'x',['h2=Vx(' num2str(d) ')']);
    
         x = X.center; v =X.vector{1}; p1 = P.center; p2 = P.vector{1};
    
         oos2 = 1/defo.kernel_size_mom^2;
    
         dpHr.center = V(x',x',p1',p2',v',oos2,v')';
         dpHr.vector{1} = dV(x',x',p1',p2',v',oos2,v',p1',v')';
         dqHr.center = dV_ad(x',x',p1',p2',v',oos2,v',p1')'+dV_2ad(x',x',p1',p2',v',oos2,v',p2',v')';
         dqHr.vector{1} = dV_ad(x',x',p1',p2',v',oos2,v',p2')';  
      
      case 2
    
         V = Kernel('Exp(-p*SqNorm2(x-y))*p0+Grad(Grad(Exp(-p*SqNorm2(x-y))*p1,y,h),h,u1)+Grad(Grad(Exp(-p*SqNorm2(x-y))*p2,y,h),h,u2)',...
         ['x=Vx(' num2str(d) ')'],['y=Vy(' num2str(d) ')'],['p0=Vy(' num2str(d) ')'],...
         ['p1=Vy(' num2str(d) ')'], ['p2=Vy(' num2str(d) ')'], ['u1=Vy(' num2str(d) ')'],...
         ['u2=Vy(' num2str(d) ')'],'p=Pm(1)',['h=Vy(' num2str(d) ')']);


         dV_ad = GradKernel(V,'x',['h1=Vx(' num2str(d) ')']);
         dV = GradKernel(dV_ad,'h1',['h2=Vx(' num2str(d) ')']);
         dV_2ad = GradKernel(dV_ad,'x',['h2=Vx(' num2str(d) ')']);
    
         x = X.center';
         y = X.center';
         u1 = X.vector{1}';
         u2 = X.vector{2}';
         p0 = P.center';
         p1 = P.vector{1}';
         p2 = P.vector{2}';

         oos2 = 1/defo.kernel_size_mom^2;
    
    
         dpHr.center = V(x,y,p0,p1,p2,u1,u2,oos2,u1)';
         dpHr.vector{1} = dV(x,y,p0,p1,p2,u1,u2,oos2,u1,p1,u1)';
         dpHr.vector{2} = dV(x,y,p0,p1,p2,u1,u2,oos2,u1,p1,u2)';
         dqHr.center = dV_ad(x,y,p0,p1,p2,u1,u2,oos2,u1,p0)'+...
         dV_2ad(x,y,p0,p1,p2,u1,u2,oos2,u1,p1,u1)'+...
         dV_2ad(x,y,p0,p1,p2,u1,u2,oos2,u1,p2,u2)';
         dqHr.vector{1} = dV_ad(x,y,p0,p1,p2,u1,u2,oos2,u1,p1)';
         dqHr.vector{2} = dV_ad(x,y,p0,p1,p2,u1,u2,oos2,u1,p2)';
    
  end
  
end

function [dqHr,dpHr] = dHr_tan_push_mat(X,P,defo)
    %here f = dHr
    [n,d]=size(X.center);
   
    if isfield(X,'vector') == 0
       m = 0;     
    else         
       m = length(X.vector);    
    end
 
switch m
    case 0
        
       x = X.center; p = P.center;
        % Calcul de A=exp(-|x_i -x_j|^2/(lam^2))
       S=zeros(n);
       for l=1:d
          S=S+(repmat(x(:,l),1,n)-repmat(x(:,l)',n,1)).^2;
       end
       A = rho(S,0,defo.kernel_size_mom);
       B = rho(S,1,defo.kernel_size_mom);
    
    
       dpHr.center = A*p;
    
    
       dp=zeros(n,d);
       for r=1:d
           % Computation of B=2*|x_i -x_j|*exp(-|x_i -x_j|^2/(lam^2))/(lam^2)
           Br=2*(repmat(x(:,r),1,n)-repmat(x(:,r)',n,1)).*B;
           dp(:,r)=dp(:,r)+ sum(p .* (Br*p) ,2);
       end
       
       dqHr.center = dp;
       
    case 1
       x = X.center; v =X.vector{1}; p1 = P.center; p2 = P.vector{1};
        % Calcul de A=exp(-|x_i -x_j|^2/(lam^2))
       S=zeros(n); % S = (|x_i-x_j|^2)_{i,j}
       for l=1:d
           S=S+(repmat(x(:,l),1,n)-repmat(x(:,l)',n,1)).^2;
       end
       Gamma_0 = rho(S,0,defo.kernel_size_mom);
       Gamma_1 = rho(S,1,defo.kernel_size_mom);
       Gamma_2 = rho(S,2,defo.kernel_size_mom);
       Gamma_3 = rho(S,3,defo.kernel_size_mom);
    
    
       %Inner product matrices
       Cxv = x*v'; %Inner products of x and v
       Cv = v*v';
       Cp1 = p1*p1';
       Cp2 = p2*p2';
       Cp12 = p1*p2';
    
       A = repmat(diag(Cxv),1,n)'-Cxv;
       B = Cxv'- repmat(diag(Cxv),1,n);
       H = (2*Gamma_1.*Cp1 + 4*Gamma_2.*A.*Cp12)-(4*Gamma_2.*B.*(Cp12')+ (8*Gamma_3.*A.*B+4*Gamma_2.*Cv).*Cp2);
       J = -(2*Gamma_1.*(Cp12')+4*Gamma_2.*A.*(Cp2));
       D = repmat(sum(J,2),1,d).*v; 
    
       C = zeros(n,d);
       E = zeros(n,d);
       for r=1:d
   
         X_r = repmat(x(:,r)',n,1)-repmat(x(:,r),1,n);
         C(:,r) = sum(H.*X_r,2);
         E(:,r) = sum(J.*X_r,2);
         
       end
    
       dpHr.center = Gamma_0*p1 +2*(Gamma_1.*A)*p2;
       dpHr.vector{1} = -2*(Gamma_1.*B)*p1 - (4*Gamma_2.*A.*B+2*Gamma_1.*Cv)*p2;
       dqHr.center = -( C +D + (2*Gamma_1.*Cp12 - 4*Gamma_2.*B.*(Cp2))*v );
       dqHr.vector{1} = E - (2*Gamma_1.*Cp2)*v; 
    
    case 2
       % Calcul de A=exp(-|x_i -x_j|^2/(lam^2))
       S=zeros(n); % S = (|x_i-x_j|^2)_{i,j}
       for l=1:d
           S=S+(repmat(X.center(:,l),1,n)-repmat(X.center(:,l)',n,1)).^2;
       end
       Gamma_0 = rho(S,0,defo.kernel_size_mom);
       Gamma_1 = rho(S,1,defo.kernel_size_mom);
       Gamma_2 = rho(S,2,defo.kernel_size_mom);
       Gamma_3 = rho(S,3,defo.kernel_size_mom);
    
    
       %Inner product matrices
       scp_x1 = X.center*X.vector{1}'; % <x_i,u^1_k>
       scp_x2 = X.center*X.vector{2}'; % <x_i,u^2_k>
       scp_11 = X.vector{1}*X.vector{1}'; % <u_i,u_k>
       scp_22 = X.vector{2}*X.vector{2}';
       scp_12 = X.vector{1}*X.vector{2}';
    
       scp_mxx = P.center*P.center'; % <p^x_i,p^x_k>
       scp_mx1 = P.center*P.vector{1}';
       scp_mx2 = P.center*P.vector{2}';
    
       scp_m11 = P.vector{1}*P.vector{1}'; %<p^u1_i,p^u2_k>
       scp_m22 = P.vector{2}*P.vector{2}';
       scp_m12 = P.vector{1}*P.vector{2}';
    
    
       A_1 = repmat(diag(scp_x1)',n,1)-scp_x1; %<x_k-x_i,u^1_k>
       A_2 = repmat(diag(scp_x2)',n,1)-scp_x2;
       B_1 =-A_1'; B_2 = -A_2';
    
       M_1 = A_1.*scp_m11 + A_2.*scp_m12;
       M_2 = A_1.*(scp_m12') + A_2.*scp_m22;
       N_1 = -M_1'; N_2 = -M_2';
    
       H = 2*(Gamma_1.*scp_mxx) + 4*Gamma_2.*((scp_mx1.*A_1 + scp_mx2.*A_2)-...
           (B_1.*(scp_mx1') + scp_11.*scp_m11 + scp_12.*scp_m12 +...
            B_2.*(scp_mx2') + (scp_12.*scp_m12)' +scp_22.*scp_m22 )) -...
            8*Gamma_3.*(B_1.*M_1 + B_2.*M_2);
       J_1 = 2*Gamma_1.*(scp_mx1') + 4*Gamma_2.*M_1;
       J_2 = 2*Gamma_1.*(scp_mx2') + 4*Gamma_2.*M_2;
       D_1 = repmat(sum(J_1,2),1,d).*X.vector{1};
       D_2 = repmat(sum(J_2,2),1,d).*X.vector{2};
    
       C = zeros(n,d); E_1 = zeros(n,d); E_2 = zeros(n,d);
    
       for r=1:d       
            X_r = repmat(X.center(:,r)',n,1)-repmat(X.center(:,r),1,n);
            C(:,r) = sum(H.*X_r,2);
            E_1(:,r) = sum(J_1.*X_r,2); E_2(:,r) = sum(J_2.*X_r,2);
       end
    

    
       dpHr.center = Gamma_0*P.center + 2*(Gamma_1.*A_1)*P.vector{1} + 2*(Gamma_1.*A_2)*P.vector{2};
       dpHr.vector{1} = -2*(Gamma_1.*B_1)*P.center- (4*(Gamma_2.*B_1.*A_1) + 2*(Gamma_1.*scp_11))*P.vector{1} -...
                     (4*(Gamma_2.*B_1.*A_2) +2*(Gamma_1.*scp_12))*P.vector{2};
       dpHr.vector{2} = -2*(Gamma_1.*B_2)*P.center- (4*(Gamma_2.*B_2.*A_1) + 2*(Gamma_1.*(scp_12')))*P.vector{1} -...
                     (4*(Gamma_2.*B_2.*A_2) +2*(Gamma_1.*scp_22))*P.vector{2};  
       dqHr.center = - (C-D_1-D_2 + (2*Gamma_1.*scp_mx1 - 4*(Gamma_2.*N_1))*X.vector{1}...
                              + (2*Gamma_1.*scp_mx2 - 4*(Gamma_2.*N_2))*X.vector{2});
       dqHr.vector{1} = -(E_1 + 2*(Gamma_1.*scp_m11)*X.vector{1} + 2*(Gamma_1.*scp_m12)*X.vector{2});
       dqHr.vector{2} = -(E_2 + 2*(Gamma_1.*(scp_m12'))*X.vector{1} + 2*(Gamma_1.*scp_m22)*X.vector{2});
end

end

function [dqHr,dpHr] = dHr_tan_normal_keops(X,P,defo)
   
   x = X.center; v =X.vector{1}; p1 = P.center; p2 = P.vector{1};
   [n,d]=size(x);
    V = Kernel('Exp(-p*SqNorm2(x-y))*p1+Grad(Grad(Exp(-p*SqNorm2(x-y))*p2,y,h),h,u)',...
   ['x=Vx(' num2str(d) ')'],['y=Vy(' num2str(d) ')'],['p1=Vy(' num2str(d) ')'],...
   ['p2=Vy(' num2str(d) ')'],['u=Vy(' num2str(d) ')'],'p=Pm(1)',['h=Vy(' num2str(d) ')']);

   
    dV_ad = GradKernel(V,'x',['h1=Vx(' num2str(d) ')']);
    dV = GradKernel(dV_ad,'h1',['h2=Vx(' num2str(d) ')']);
    dV_2ad = GradKernel(dV_ad,'x',['h2=Vx(' num2str(d) ')']);
    
    oos2 = 1/defo.kernel_size_mom^2;
    NV = sqrt(Inn(v,v)); %norms of d_i
    pjp2 = PJ(p2,v);  %project p2
    Dxv_d = dV(x',x',p1',pjp2',v',oos2,v',p1',v')';
    
    dpHr.center = V(x',x',p1',pjp2',v',oos2,v')';
    dpHr.vector{1} = PJ(Dxv_d,v);
    dqHr.center = dV_ad(x',x',p1',pjp2',v',oos2,v',p1')'+dV_2ad(x',x',p1',pjp2',v',oos2,v',pjp2',v')';
%     dqHr.vector{1} = dV_ad(x',x',p1',pjp2',v',oos2,v',pjp2')' - Inn(v,p2).*Dxv_d - Inn(v,Dxv_d).*p2;
    dqHr.vector{1} = dV_ad(x',x',p1',pjp2',v',oos2,v',pjp2')' - (Inn(v,p2).*PJ(Dxv_d,v) + Inn(v,Dxv_d).*pjp2)./(NV.^2);
%    PJ(S,T)
end

function [dqHr,dpHr] = dHr_tan_normal_mat(X,P,defo)


%v = v./NV; %normalize v
x = X.center; v =X.vector{1}; p1 = P.center; p2 = P.vector{1};

NV = sqrt(Inn(v,v)); %norms of d_i
Pd2 = PJ(p2,v); %project p2

[n,d]=size(x);
  


    % Calcul de A=exp(-|x_i -x_j|^2/(lam^2))
    S=zeros(n); % S = (|x_i-x_j|^2)_{i,j}
    for l=1:d
        S=S+(repmat(x(:,l),1,n)-repmat(x(:,l)',n,1)).^2;
    end
    Gamma_0 = rho(S,0,defo.kernel_size_mom);
    Gamma_1 = rho(S,1,defo.kernel_size_mom);
    Gamma_2 = rho(S,2,defo.kernel_size_mom);
    Gamma_3 = rho(S,3,defo.kernel_size_mom);
 
    
    %Inner product matrices
    Cxv = x*v'; %Inner products of x and v <x_i,v_k>
    Cv = v*v';    % <v_i,v_k>
    Cp1 = p1*p1'; % 
  
    Cpd12 = p1*Pd2';
    Cpd2 = Pd2*Pd2';
    
    A = repmat(diag(Cxv),1,n)'-Cxv;
    B = Cxv'- repmat(diag(Cxv),1,n);
 
    dpHr.center  = Gamma_0*p1 +2*(Gamma_1.*A)*Pd2;
    Dvd = -2*(Gamma_1.*B)*p1 - (4*Gamma_2.*A.*B+2*Gamma_1.*Cv)*Pd2;
    dpHr.vector{1}   = PJ(Dvd,v);
    

    H = (2*Gamma_1.*Cp1 + 4*Gamma_2.*A.*Cpd12)-(4*Gamma_2.*B.*(Cpd12')+ (8*Gamma_3.*A.*B+4*Gamma_2.*Cv).*Cpd2);
    J = -(2*Gamma_1.*(Cpd12')+4*Gamma_2.*A.*(Cpd2));
    D = repmat(sum(J,2),1,d).*v;

    C = zeros(n,d);
    E = zeros(n,d);

    for r=1:d
         X_r = repmat(x(:,r)',n,1)-repmat(x(:,r),1,n);
         C(:,r) = sum(H.*X_r,2);
         E(:,r) = sum(J.*X_r,2);       
    end
    
   dqHr.center = -(C +D + (2*Gamma_1.*Cpd12 - 4*Gamma_2.*B.*(Cpd2))*v);
%     dqHr.vector{1} =E-(2*Gamma_1.*Cpd2)*v- Inn(Dvd,v).*p2 - Inn(p2,v).*Dvd;
%      dqHr.vector{1} =E-(2*Gamma_1.*Cpd2)*v- (Inn(Dvd,v).*PJ(p2,v) + Inn(p2,v).*PJ(Dvd,v))./(NV.^2);
   dqHr.vector{1} =E-(2*Gamma_1.*Cpd2)*v- (Inn(Dvd,v).*Pd2 + Inn(p2,v).*PJ(Dvd,v))./(NV.^2);
end

% function IN = Inn(S,T)
% [~,d] = size(S);
% IN = repmat(sum(S.*T,2),1,d);
% end

% function PJ = PJ(S,T)
% PJ=S-Inn(S,T).*T;
% end