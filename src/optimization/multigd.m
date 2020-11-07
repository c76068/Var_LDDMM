% function [shooting,E] = multigd(matching,opt,defo,objfun)
function[shooting,hist,E]=multigd(Template,Target,mom,opt,defo,objfun)

temp=zeros(1,3);
[temp(1),~] = size(opt.maxiteration);
[temp(2),~] = size(objfun.kernel_size_geom);
[temp(3),~] = size(objfun.kernel_size_grass);
% temp
% pause()
if ~(temp(1)==temp(2) && temp(2)==temp(3))
    fprintf('error')
     pause()
     return
end

% p1 = matching.mom1;
% p2 = matching.mom2;
gopt =opt;
gobjfun = objfun;
% gmatching =matching;
E=[];

for i=1:temp(1)
 gopt.maxiteration = opt.maxiteration(i,:);
 gobjfun.kernel_size_geom = objfun.kernel_size_geom(i,:);
 gobjfun.kernel_size_grass = objfun.kernel_size_grass(i,:);
 
 [shooting,hist,Ei]=gradientdescent(Template,Target,mom,gopt,defo,gobjfun);
 mom = shooting.P;
%  mom.p1 = shooting.mom1{1};
%  mom.p2 = shooting.mom2{1};
 E=[E 0 Ei];
end

end

function [shooting,hist,E] = gradientdescent(Template,Target,mom,opt,defo,objfun)
%Input:
%defo : structure containing the parameters of deformations (kernel_size_mom,method,nstep,...)

%Output:
%shooting.X: optimal trajectory of Template
%shooting.P: momentum evolution with optimal initial momentum
%shooting.del: gradient descend step 
%shooting.cost
%hist.P_ini
%hist.grad

% del = opt.stepsize;
% [n,d]=size(Template.center);
% [Ny,~] = size(Target.center);

if strcmp(defo.method,'matlab')
 defo.Kq = Kqop(Template,defo);
end

Eold = -5;

shooting = cost(Template,mom,Target,objfun,defo);
E(1) = shooting.cost;
shooting.del = opt.stepsize;
% P_old = shooting.P{1};

iter=1;
hist.P_ini{iter} = mom;

for iter = 1:opt.maxiteration
% while abs(E(iter)-Eold)>1e-6 && iter < opt.maxiteration
%       while  iter < opt.maxiteration
        Eold = E(iter);
       
        %compute gradient        
        
  
        shooting.grad = gsgrad(Template,hist.P_ini{iter},Target,objfun,defo); 
        hist.grad{iter} =shooting.grad;
        
        
        shooting = adjstep(shooting,Target,hist.P_ini{iter},defo,objfun);
        
        %If it decrease the energy, then increase the stepsize by 1.2
        if shooting.cost-Eold<0
            shooting.del = 1.2*shooting.del;
        end
        
        %If the energy does not decrease, divide the stepsize by 2
        sto =0;
        while shooting.cost-Eold>0 && sto<41
            sto =sto+1;
            shooting.del=shooting.del/2;
            shooting = adjstep(shooting,Target,hist.P_ini{iter},defo,objfun);
        end
        
        iter=iter+1;
 
        E(iter) = shooting.cost;
        shooting.iter = iter;
%  shooting.stepsize =del;
        hist.P_ini{iter} = shooting.P{1};
        

optcon_g = norm([shooting.grad.center(:);shooting.grad.vector{1}(:);shooting.grad.vector{2}(:)],inf); 
        
fprintf('iter=%d\t',iter-1)
fprintf('stepsize=%d\t',shooting.del)
fprintf('Energy=%12.6f\t',E(iter))
fprintf('Energy decreased by =%12.6f\n',E(iter-1)-E(iter))
% fprintf('norm of grad_p.center = %12.6f\t',norm(shooting.grad.center,'fro'))
% fprintf('norm of grad_p.vector = %12.6f\n',norm(shooting.grad.vector{1},'fro')+norm(shooting.grad.vector{2},'fro'))
fprintf('sup norm of gradient = %12.6f\t',optcon_g)
fprintf('\r\n')

 
 if abs(E(iter)-Eold)<1e-6 && optcon_g<1e-6
    break; 
 end
 
 if shooting.del<1e-16
    break; 
 end
end
  
  
end
 

function [adjoutput] = adjstep(shooting,Target,P_old,defo,objfun)
%Doing one step of gradient descend and adjust the stepsize if there are 
%NAN in forward equations
%
%Input:
% ept.X ept.P: evolution of X and P 
% ept.grad: gradient of cost function
% K_q : symmetric positive definite matrix used to compute cost
% del: stepsize
%
%Output:
%adjoutput.X, adjoutput.P, adjoutput.cost, adjoutput.del

% P_old = shooting.P{1};
Template = shooting.X{1};

P_new = add_XY_h(P_old,shooting.grad,-shooting.del);

ept =cost(Template,P_new,Target,objfun,defo); 

Tes = ept.X{defo.nb_euler_steps+1}.center+ept.X{defo.nb_euler_steps+1}.vector{1}+...
                ept.P{defo.nb_euler_steps+1}.center+ept.P{defo.nb_euler_steps+1}.vector{1}+...
                ept.P{defo.nb_euler_steps+1}.vector{2};
          
            tinf =0;
       while sum(sum(isnan(Tes)+isinf(Tes)))>0
            shooting.del = shooting.del/2;
            
            P_new = add_XY_h(P_old,shooting.grad,-shooting.del);
            ept =cost(Template,P_new,Target,objfun,defo);

            Tes = ept.X{defo.nb_euler_steps+1}.center+ept.X{defo.nb_euler_steps+1}.vector{1}+...
                ept.P{defo.nb_euler_steps+1}.center+ept.P{defo.nb_euler_steps+1}.vector{1}+...
                ept.P{defo.nb_euler_steps+1}.vector{2};
            tinf =tinf +1;
       end
 adjoutput = ept;
 adjoutput.del = shooting.del;
 adjoutput.grad = shooting.grad;

end

% function lsoutput=linesearch(shooting,Template,Target,K_q,del,defo,objfun)
%  %This is the line search function
%  
%  
% [adjoutput,del] = adjstep(shooting,Template,Target,K_q,del,defo,objfun);
% 
% %  mom_old.p1 = adjoutput.mom1{1};
% %  mom_old.p2 = adjoutput.mom2{1};
% 
% mom_old.p1 = shooting.mom1{1};
% mom_old.p2 = shooting.mom2{1};
%  Energy = adjoutput.cost;
%  
%  eta = zeros(2,1);
%  
% 
%  lsoutput.points = adjoutput.x;
%  lsoutput.direction = adjoutput.v;
%  lsoutput.mom1 = adjoutput.mom1;
%  lsoutput.mom2 = adjoutput.mom2;
%  
%  for i=1:4
%   
% %  switch i  
% %      case 1
% %       eta(1) = del(1);
% %       eta(2) = del(2)/2;
% %      case 2
% %       eta(1) = del(1);
% %       eta(2) = 2*del(2);
% %      case 3
% %       eta(1) = del(1)/2;
% %       eta(2) = del(2);   
% %      case 4
% %       eta(1) = 2*del(1);
% %       eta(2) = del(2);    
% %  end
% 
% switch i  
%      case 1
%       eta(1) = del(1)/2;
%       eta(2) = del(2)/2;
%      case 2
%       eta(1) = 1.25*del(1);
%       eta(2) = 1.25*del(2);
%      case 3
%       eta(1) = 1.5*del(1);
%       eta(2) = 1.5*del(2);   
%      case 4
%       eta(1) = 2*del(1);
%       eta(2) = 2*del(2);    
%  end
%  mom_new.p1 = mom_old.p1 - eta(1)*(shooting.grad{1});
%  mom_new.p2 = mom_old.p2 - eta(2)*(shooting.grad{2});
%  
%  ept =cost(Template,mom_new,Target,objfun,defo,K_q);
% 
%  if ept.cost < Energy 
%      Energy = ept.cost;
%      del = eta;
%      lsoutput.points = ept.x;
%      lsoutput.direction = ept.v;
%      lsoutput.mom1 = ept.mom1;
%      lsoutput.mom2 = ept.mom2;
%  end
%  
%  end       
%  
% 
%         lsoutput.energy = Energy;
%         lsoutput.stepsize = del;       
% end
% 
