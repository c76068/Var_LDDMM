%% Registration
clear all
close all


addpath(genpath('../src/'))

%Load the amygdala surfaces--------------------
template = import_shape_ply('../Data/amygdala/0354387_1_6_amyg_R_up2NN.ply');
target = import_shape_ply('../Data/amygdala/1449400_1_4_amyg_R_up2NN.ply');


figure(1)
hold on

trisurf(template.G,template.x(:,1),template.x(:,2),template.x(:,3),'Facecolor','r','FaceAlpha',0.7)
trisurf(target.G,target.x(:,1),target.x(:,2),target.x(:,3),'Facecolor','b','FaceAlpha',0.7)
view(3)

%Match the center of two surfaces-----------------------------------
template.x = template.x + (mean(target.x) - mean(template.x));

figure(2)
hold on

trisurf(template.G,template.x(:,1),template.x(:,2),template.x(:,3),'Facecolor','r','FaceAlpha',0.5)
trisurf(target.G,target.x(:,1),target.x(:,2),target.x(:,3),'Facecolor','b','FaceAlpha',0.5)
view(3)

%Transform the surfaces to discrete varifolds---------------------------------------
[X.center,X.vector] = mesh2var(template.x,template.G);
[Y.center,Y.vector] = mesh2var(target.x,target.G); 
 
%Set options for deformation and numerical ODE----------------------------------------
defo.kernel_size_mom = 4.75; %deformation kernel size
defo.nb_euler_steps =15; % number of steps in the (for||back)ward integration
defo.odemethod = 'rk4'; %numerical scheme for ODE
defo.method='keops';% compute kernel operation in deformation using matlab or keops, possible values: 'keops' or 'matlab'

%Set options for data attachment term-------------------------
objfun.kernel_geom = 'gaussian';% spatial kernel type, possible values:'gaussian' or 'cauchy'
objfun.kernel_size_geom = [3];% spatial kernel size
objfun.kernel_grass = 'gaussian_oriented'; %Grassmanian/orientation kernel type, possible values: 'linear', 'gaussian_oriented', 'gaussian_unoriented', 'binet'
objfun.kernel_size_grass = [1.1]; %Grassmanian/orientation kernel size
objfun.method='keops'; % compute kernel operation in data attachment term using matlab or keops, possible values: 'keops' or 'matlab'
objfun.lambda=10;

%Set options for L-BFGS-----------------------------
options.record_history = true;
options.maxit = 500; % max number of iterations(default 1000) (applies to each starting vector)
options.nvec = 10; %0 for full BFGS matrix update, otherwise specifies number of vectors to save and use in the limited memory updates (default: 0)
options.prtlevel = 2; %one of 0 (no printing), 1 (minimal), 2 (verbose) (default: 1)


%registration and optimization----------------------
[P_op,summary]= registration(X,Y,defo,objfun,options);

% shoot the discrete varifold X and flow the vertices of the source surface------------- 
[obj_evol, X_evol, P_evol] = shoot_and_flow(X,P_op,defo,template.x);

save('example1_keops_float.mat')

%% Display registration result: evolution of deformed surface
clear all
close all

load('example1_keops_float.mat')

fig=figure(1);
for t = [1 4 8 12 16]

clf

trisurf(target.G,target.x(:,1),target.x(:,2),target.x(:,3),'Facecolor','r','FaceAlpha',0.5)
view(70,2)
 
hold on

trisurf(template.G,obj_evol{t}(:,1),obj_evol{t}(:,2),obj_evol{t}(:,3),'Facecolor','b','FaceAlpha',0.3)

xlabel('x')
ylabel('y')
zlabel('z')

axis([85 110 123 142 98 118])
axis off

print(fig,['example1_30',num2str(t)],'-dpng')
pause(0.3)

end


%% Display registration result: evolution of deformed discrete varifold
clear all
close all

load('example1_keops_float.mat')


X.normal = cross(X.vector{1},X.vector{2});
Y.normal = cross(Y.vector{1},Y.vector{2});


rng(5)
ind = round(linspace(1,size(Y.center,1),400));

Z.center = Y.center(ind,:);
Z.vector{1} = Y.vector{1}(ind,:);
Z.vector{2} = Y.vector{2}(ind,:);
Z.normal = Y.normal(ind,:);


rng(5)
ind = round(linspace(1,size(X.center,1),400));


fig=figure(2);
for t = [1 4 8 12 16]

    
    clf
    hold on
    view(70,2)
    for i=1:size(Z.center,1)
        A = 0.9*orth([Z.vector{1}(i,:)',Z.vector{2}(i,:)']);
        B = [Z.center(i,:)+A(:,1)';Z.center(i,:)+A(:,2)';Z.center(i,:)-A(:,1)';Z.center(i,:)-A(:,2)'];
        G = [1 2 3;3 4 1];
        trisurf(G,B(:,1),B(:,2),B(:,3),'Facecolor','r','FaceAlpha',0.3,'edgecolor','none')
        quiver3(Z.center(:,1),Z.center(:,2),Z.center(:,3),Z.normal(:,1),Z.normal(:,2),Z.normal(:,3),'color','r','LineWidth',1.6) 
    end

    U.center = X_evol{t}.center(ind,:);
    U.vector{1} = X_evol{t}.vector{1}(ind,:);
    U.vector{2} = X_evol{t}.vector{2}(ind,:);
    U.normal = cross(U.vector{1},U.vector{2});
    
    for i=1:size(U.center,1)
        A = 0.9*orth([U.vector{1}(i,:)',U.vector{2}(i,:)']);
        B = [U.center(i,:)+A(:,1)';U.center(i,:)+A(:,2)';U.center(i,:)-A(:,1)';U.center(i,:)-A(:,2)'];
        G = [1 2 3;3 4 1];
        trisurf(G,B(:,1),B(:,2),B(:,3),'Facecolor','b','FaceAlpha',0.3,'edgecolor','none')
        quiver3(U.center(:,1),U.center(:,2),U.center(:,3),U.normal(:,1),U.normal(:,2),U.normal(:,3),'color','b','LineWidth',1.6)  
    end
    
xlabel('x')
ylabel('y')
zlabel('z')
axis([85 110 123 142 98 118])
axis off
print(fig,['example1_40',num2str(t)],'-dpng')
pause(0.3)
end

