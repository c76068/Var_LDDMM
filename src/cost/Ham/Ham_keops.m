function  [ham,ham_grad] = Ham_keops(X,P,defo)
   
 switch defo.action
    case 'normalized'
        HAM = @Ham_normal_keops;
    case 'pushforward'
        HAM = @Ham_push_keops;
 end
 [ham,ham_grad] = HAM(X,P,defo);
 
end


function  [ham,ham_grad] = Ham_push_keops(X,P,defo)

   
%    [n,d]=size(X.center);
%    m = length(X.vector);
   
    [~,d]=size(X.center);
   
    if isfield(X,'vector') == 0
       m = 0;     
    else         
       m = length(X.vector);    
    end
   
   switch m
       case 0
           V = Kernel('Exp(-p*SqNorm2(x-y))*p1',['x=Vx(' num2str(d) ')'],['y=Vy(' num2str(d) ')'],['p1=Vy(' num2str(d) ')'],'p=Pm(1)');
%            dV_ad = GradKernel(V,'x',['h1=Vx(' num2str(d) ')']);
           
           x = X.center'; p = P.center'; 
           oos2 = 1/defo.kernel_size_mom^2;
         
           Vx = V(x,x,p,oos2)';
           ham = 0.5*(sum(sum(P.center.*Vx)));
           ham_grad.center = Vx;
           
       case 1
             
           V = Kernel('Exp(-p*SqNorm2(x-y))*p1+Grad(Grad(Exp(-p*SqNorm2(x-y))*p2,y,h),h,u)',...
           ['x=Vx(' num2str(d) ')'],['y=Vy(' num2str(d) ')'],['p1=Vy(' num2str(d) ')'],...
           ['p2=Vy(' num2str(d) ')'],['u=Vy(' num2str(d) ')'],'p=Pm(1)',['h=Vy(' num2str(d) ')']);

   
           dV_ad = GradKernel(V,'x',['h1=Vx(' num2str(d) ')']);
           dV = GradKernel(dV_ad,'h1',['h2=Vx(' num2str(d) ')']);

           x = X.center';
           y = X.center';
           u = X.vector{1}';
    
           p1 = P.center';
           p2 = P.vector{1}';
           pjp2 = PJ(P.vector{1},X.vector{1})';
           oos2 = 1/defo.kernel_size_mom^2;
           dum = u;
    
           switch defo.action     
              case 'pushforward'
                  
                 Vx = V(x,y,p1,p2,u,oos2,u)';            
                 dVx_u = dV(x,y,p1,p2,u,oos2,u,p1,u)';
                 
              case 'normalized'
%             dobj = V(x,y,p1,PJ(P.vector{1},X.vector{1})',u,oos2,dum)'; 
                 Vx = V(x,y,p1,pjp2,u,oos2,dum)';
                 dVx_u = dV(x,y,p1,pjp2,u,oos2,dum,dum,u)';
                 dVx_u = PJ(dVx_u,X.vector{1});
           end
           
           
    
           ham = 0.5*(sum(sum(P.center.*Vx)) + sum(sum(P.vector{1}.*dVx_u)));
    
           ham_grad.center = Vx;
           ham_grad.vector{1} = dVx_u;
   
       case 2
   
           V = Kernel('Exp(-p*SqNorm2(x-y))*p0+Grad(Grad(Exp(-p*SqNorm2(x-y))*p1,y,h),h,u1)+Grad(Grad(Exp(-p*SqNorm2(x-y))*p2,y,h),h,u2)',...
           ['x=Vx(' num2str(d) ')'],['y=Vy(' num2str(d) ')'],['p0=Vy(' num2str(d) ')'],...
           ['p1=Vy(' num2str(d) ')'], ['p2=Vy(' num2str(d) ')'], ['u1=Vy(' num2str(d) ')'],...
           ['u2=Vy(' num2str(d) ')'],'p=Pm(1)',['h=Vy(' num2str(d) ')']);


           dV_ad = GradKernel(V,'x',['h1=Vx(' num2str(d) ')']);
           dV = GradKernel(dV_ad,'h1',['h2=Vx(' num2str(d) ')']);
          %dV_2ad = GradKernel(dV_ad,'x',['h2=Vx(' num2str(d) ')']);
    
           x = X.center';
           y = X.center';
           u1 = X.vector{1}';
           u2 = X.vector{2}';
           p0 = P.center';
           p1 = P.vector{1}';
           p2 = P.vector{2}';

           oos2 = 1/defo.kernel_size_mom^2;
    
    
           Vx = V(x,y,p0,p1,p2,u1,u2,oos2,u1)';
           dVx_ad_p1 = dV_ad(x,y,p0,p1,p2,u1,u2,oos2,u1,p1)';
           dVx_ad_p2 = dV_ad(x,y,p0,p1,p2,u1,u2,oos2,u1,p2)';
    
           ham= 0.5*(sum(sum(P.center.*Vx)) + sum(sum(X.vector{1}.*dVx_ad_p1)) + sum(sum(X.vector{2}.*dVx_ad_p2)) );
    
           if nargout >1
              ham_grad.center = V(x,y,p0,p1,p2,u1,u2,oos2,u1)';
              ham_grad.vector{1} = dV(x,y,p0,p1,p2,u1,u2,oos2,u1,p1,u1)';
              ham_grad.vector{2} = dV(x,y,p0,p1,p2,u1,u2,oos2,u1,p1,u2)';
           end  
    
   end
   
end

function  [ham,ham_grad] = Ham_normal_keops(X,P,defo)

   
   [n,d]=size(X.center);
   
    V = Kernel('Exp(-p*SqNorm2(x-y))*p1+Grad(Grad(Exp(-p*SqNorm2(x-y))*p2,y,h),h,u)',...
   ['x=Vx(' num2str(d) ')'],['y=Vy(' num2str(d) ')'],['p1=Vy(' num2str(d) ')'],...
   ['p2=Vy(' num2str(d) ')'],['u=Vy(' num2str(d) ')'],'p=Pm(1)',['h=Vy(' num2str(d) ')']);

   
    dV_ad = GradKernel(V,'x',['h1=Vx(' num2str(d) ')']);
    dV = GradKernel(dV_ad,'h1',['h2=Vx(' num2str(d) ')']);

    x = X.center';
    y = X.center';
    u = X.vector{1}';
    
    p1 = P.center';
    p2 = P.vector{1}';
    
    oos2 = 1/defo.kernel_size_mom^2;
    pjp2 = PJ(p2',u')';
    
    Vx = V(x,y,p1,pjp2,u,oos2,u)';
    PJdVx_u = PJ(dV(x,x,p1,pjp2,u,oos2,u,p1,u)',u');
    
    %dp1Hr = V(x',x',p1',pjp2',v',oos2,v')';
    %dp2Hr = PJ(Dxv_d,v);
    
    ham= 0.5*(sum(sum(P.center.*Vx)) + sum(sum(P.vector{1}.*PJdVx_u)));
    
    ham_grad.center = Vx;
    ham_grad.vector{1} = PJdVx_u;
   
end

% function IN = Inn(S,T)
% [~,d] = size(S);
% IN = repmat(sum(S.*T,2),1,d);
% end
% 
% function PJ = PJ(S,T)
% PJ=S-Inn(S,T).*T;
% end