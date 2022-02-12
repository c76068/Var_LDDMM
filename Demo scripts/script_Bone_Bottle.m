%% Quantization-----------------------------------
clear all
close all

addpath(genpath('../src/'))
 
%Load curves data----
template = import_shape_vtk('../Data/bone_bottle/bottle.vtk');
template.x=template.x(:,1:2);

target = import_shape_vtk('../Data/bone_bottle/bone.vtk'); 
target.x=target.x(:,1:2);


%Plot bottle and bone curves------
fig = figure(1);
hold on 
plot([template.x(:,1) ; template.x(1,1)] ,[template.x(:,2) ; template.x(1,2)],'b','LineWidth',2)
plot([target.x(:,1) ; target.x(1,1)],[target.x(:,2) ; target.x(1,2)],'r','LineWidth',2)
axis([20 140 -10 170])
axis off


%Transform the curves to discrete varifolds---------------------------------------
[X.center,X.vector] = mesh2var(template.x,template.G);
[Y.center,Y.vector] = mesh2var(target.x,target.G); 



%Set options for data attachment term-------------------------
objfun.kernel_geom = 'gaussian';% spatial kernel type, possible values:'gaussian' or 'cauchy'
objfun.kernel_size_geom = [10];% spatial kernel size
objfun.kernel_grass = 'gaussian_oriented';%Grassmanian/orientation kernel type, possible values: 'linear', 'gaussian_oriented', 'gaussian_unoriented', 'binet'
objfun.kernel_size_grass = [2];%Grassmanian/orientation kernel size
objfun.method = 'matlab';% compute kernel operation in data attachment term using matlab or keops, possible values: 'keops' or 'matlab'

%Set options for L-BFGS-----------------------------
options.record_history =true;
options.maxit = 1000; % max number of iterations(default 1000) (applies to each starting vector)
options.nvec = 10; %0 for full BFGS matrix update, otherwise specifies number of vectors to save and use in the limited memory updates (default: 0)
options.prtlevel = 2; %one of 0 (no printing), 1 (minimal), 2 (verbose) (default: 1)

X_sp = cell(1,3);
num = [25 40 150];
for i=1:3
    N = num(i);
    %uniformly subsample curve and transform it to discrete varifold as
    %initial varifold for optimization----------------
    ind = round(linspace(1,size(X.center,1),N));
    [X_rs.center, X_rs.vector] = mesh2var(template.x(ind,:), [1:N;2:N 1]');
   
    %Projection optimization--------------------
    X_sp{i} = proj2_M_dirac(X,X_rs,objfun,options);
end



%plot compressed varifolds-------------------
for i = 1:3 
    N = num(i);
    fig = figure(2);
    clf
    hold on
    quiver(X_sp{i}.center(:,1),X_sp{i}.center(:,2),X_sp{i}.vector{1}(:,1),X_sp{i}.vector{1}(:,2),0,'Linewidth',2)
    axis([20 140 -10 170])
    axis off

    print(fig,['Bottle_quan',num2str(N)],'-dpng')
end

clear p i ind

save('Bone_Bottle_sparse_unif_geogauss_10_grassgaussori_2.mat','X_sp','Y', 'template', 'target')

%% Registration---------------------------------------


clear all
close all

addpath(genpath('../src/'))

load('Bone_Bottle_sparse_unif_geogauss_10_grassgaussori_2.mat')

defo.kernel_size_mom = 10;
defo.nb_euler_steps =30;
defo.method='matlab';
defo.odemethod = 'rk4';

objfun.kernel_geom = 'gaussian';
objfun.kernel_size_geom = [15;10];%Two stage of registrations with kernel sizes 15 and 10
objfun.kernel_grass = 'gaussian_oriented';
objfun.kernel_size_grass = [2;2];%Two stage of registrations with kernel sizes 2 and 2
objfun.method = 'matlab';
objfun.lambda=1;

options.record_history =true;
options.maxit = 500;
options.nvec = 10;
options.prtlevel = 2;

tic
parfor i=1:3
   
   P = struct();
   P.center = zeros(size(X_sp{i}.center));
   P.vector{1} = zeros(size(X_sp{i}.center));
      
   [P_sp{i},summary_sp{i}] = multiBFGS(X_sp{i},Y,P,defo,objfun,options);
  
end
toc
total_time = toc;


save('Bone_Bottle_temp_sparse_matching_ker_10_itr500_lam1.mat')

%% Display registration results----------------------------

clear all
close all

addpath(genpath('../src/'))
load('Bone_Bottle_temp_sparse_matching_ker_10_itr500_lam1.mat')
%load('Bone_Bottle_sparse_unif_geogauss_10_grassgaussori_2.mat')

objfun_1 = objfun;
objfun.kernel_size_geom = objfun.kernel_size_geom(2,:);
objfun.kernel_size_grass = objfun.kernel_size_grass(2,:);

%plot registration results:------------

[Xgr,Ygr]=ndgrid(40:2:120,-10:2:170);
nxgrid=size(Xgr,1);
nygrid=size(Xgr,2);
gr=[Xgr(:) Ygr(:)];

M = [25 40 150];
for i= 1:3

    ind = M(i);
    Z = X_sp{i}; P_op = P_sp{i};
    obj_evol = shoot_and_flow(Z,P_op,defo,template.x);
    gr_evol = shoot_and_flow(Z,P_op,defo,gr);


    fig = figure(1);
    for t=defo.nb_euler_steps+1

        clf
        hold on 
        plot([obj_evol{t}(:,1) ; obj_evol{t}(1,1)] ,[obj_evol{t}(:,2) ; obj_evol{t}(1,2)],'b','LineWidth',2)
        plot([target.x(:,1) ; target.x(1,1)],[target.x(:,2) ; target.x(1,2)],'r','LineWidth',2)

        axis([20 140 -10 170])
        axis off

        xt=reshape(gr_evol{t}(:,1),size(Xgr));
        yt=reshape(gr_evol{t}(:,2),size(Ygr));

        for x=1:nxgrid
           h=plot(xt(x,:),yt(x,:));
           set(h,'Color',.6*[1,1,1]);
        end

        for y=1:nygrid
           h=plot(xt(:,y),yt(:,y));
           set(h,'Color',.6*[1,1,1]);
        end

        plot([template.x(:,1) ; template.x(1,1)] ,[template.x(:,2) ; template.x(1,2)],'b--')
        pause(0.03)
    end
    print(fig,['Bottle_quan_matching',num2str(M(i))],'-dpng')

end

