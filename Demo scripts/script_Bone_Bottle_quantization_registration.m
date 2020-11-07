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
X = curve2var(template); Y = curve2var(target);

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

%Project varifold X to space of M diracs, M = 3,...,size(X.center,1)
norm_X = sqrt(shape_scp_prodspace_keops(X,X,objfun)); % varifold norm of X 

Bottle = cell(1,size(X.center,1));
X_rs = cell(1,size(X.center,1));
X_sp = cell(1,size(X.center,1));


parfor i=3:size(X.center,1)
    %uniformly subsample curve and transform it to discrete varifold as
    %initial varifold for optimization----------------
    ind = round(linspace(1,size(X.center,1),i));
    Bottle{i}.x = template.x(ind,:);  Bottle{i}.G = [1:i;2:i 1]';
    X_rs{i} = curve2var(Bottle{i});
    dist_XaX = sqrt(varifoldnorm(X_rs{i},X,objfun)); % distance between uniformly subsampled X and X 
    unif_rel_err_X(i) = dist_XaX/norm_X; % relative error 
    
    %Projection optimization--------------------
    X_sp{i} = proj2_M_dirac(X,X_rs{i},objfun,options);
    dist_XaX = sqrt(varifoldnorm(X_sp{i},X,objfun)); % distance between sparse X and X 
    rel_err_X(i) = dist_XaX/norm_X; % relative error 

end


%Project varifold Y to space of M diracs, M = 3,...,size(Y.center,1)
norm_Y = sqrt(shape_scp_prodspace_keops(Y,Y,objfun)); % square norm of Y 

Bone = cell(1,size(Y.center,1));
Y_rs = cell(1,size(Y.center,1));
Y_sp = cell(1,size(Y.center,1));

parfor i=3:size(Y.center,1)
    ind = round(linspace(1,size(Y.center,1),i));
    Bone{i}.x = target.x(ind,:);  Bone{i}.G = [1:i;2:i 1]';
    Y_rs{i} = curve2var(Bone{i});
    dist_YaY = sqrt(varifoldnorm(Y_rs{i},Y,objfun)); 
    unif_rel_err_Y(i) = dist_YaY/norm_Y; 

    Y_sp{i} = proj2_M_dirac(Y,Y_rs{i},objfun,options);
    dist_YaY = sqrt(varifoldnorm(Y_sp{i},Y,objfun)); 
    rel_err_Y(i) = dist_YaY/norm_Y;
end


clear p i ind

save('Bone_Bottle_sparse_unif_geogauss_10_grassgaussori_2.mat')


%plot compressed varifolds-------------------
for N = [25 40 150]  
    fig = figure(2);
    clf
    hold on
    quiver(X_sp{N}.center(:,1),X_sp{N}.center(:,2),X_sp{N}.vector{1}(:,1),X_sp{N}.vector{1}(:,2),0,'Linewidth',2)
    axis([20 140 -10 170])
    axis off

    print(fig,['Bottle_quan',num2str(N)],'-dpng')
end

%plot relative quantization error------------ 
fig = figure(3);
hold on
t= 10;
plot(t:length(rel_err_X),rel_err_X(t:end),'b','LineWidth',2)
plot(t:length(rel_err_X),unif_rel_err_X(t:end),'g','LineWidth',2)
legend('Quantization','Uniform subsample')
xlabel('N')
ylabel('Relative quantization error')
set(gca,'FontSize',13,'FontWeight','bold','Fontname','Ubuntu');
print(fig,'Bone_Bottle_quan_err','-dpng')

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

M = [5:5:50 60:10:200 220:20:368 368];

tic
parfor i=1:length(M) 
   ind = M(i); 
   
   P = struct();
   P.center = zeros(size(X_sp{ind}.center));
   P.vector{1} = zeros(size(X_sp{ind}.center));
   
   
   [P_rs{i},summary_rs{i}] = multiBFGS(X_rs{ind},Y,P,defo,objfun,options);    
   [P_sp{i},summary_sp{i}] = multiBFGS(X_sp{ind},Y,P,defo,objfun,options);
  
end
toc
total_time = toc;


save('Bone_Bottle_temp_sparse_matching_ker_10_itr500_lam1.mat')

%% Display registration results----------------------------

clear all
close all

addpath(genpath('../src/'))
load('Bone_Bottle_temp_sparse_matching_ker_10_itr500_lam1.mat')

objfun_1 = objfun;
objfun.kernel_size_geom = objfun.kernel_size_geom(2,:);
objfun.kernel_size_grass = objfun.kernel_size_grass(2,:);

clear i ind rel_err_X rel_err_Y total_time_X total_time_Y total_time dist_XY norm_X norm_Y

for i=1:length(M)
    %quantization
    Z = X_sp{M(i)}; P_op = P_sp{i};
    Xt = shoot_and_flow(Z,P_op,defo,template.x);
    cur_deformed.x = Xt{end}; cur_deformed.G = template.G;
    U = curve2var(cur_deformed); %transform the discrete curve into 1-D varifold
    E_sp(i) = Ham(Z,P_op,defo) + objfun.lambda*varifoldnorm(U,Y,objfun);
    
    %%uniformly subsample
    Z = X_rs{M(i)}; P_op = P_rs{i};
    Xt = shoot_and_flow(Z,P_op,defo,template.x);
    cur_deformed.x = Xt{end}; cur_deformed.G = template.G;
    U = curve2var(cur_deformed); %transform the discrete curve into 1-D varifold
    E_rs(i) = Ham(Z,P_op,defo) + objfun.lambda*varifoldnorm(U,Y,objfun);
    
    
end

%plot difference of energies----------
fig=figure(1);
hold on
plot(M,E_sp-E_sp(end),'b','Linewidth',2)
plot(M,E_rs-E_rs(end),'g','Linewidth',2)
legend('Quantization','Uniform subsample')
set(gca,'FontSize',13,'FontWeight','bold','Fontname','Ubuntu');
xlabel('N')
ylabel('E(v^N) - E(v^*)')
axis([-1 400 0 6000])

%plot registration results:------------

[Xgr,Ygr]=ndgrid(40:2:120,-10:2:170);
nxgrid=size(Xgr,1);
nygrid=size(Xgr,2);
gr=[Xgr(:) Ygr(:)];


for i= [5 8 20]

    ind = M(i);
    Z = X_sp{ind}; P_op = P_sp{i};
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

