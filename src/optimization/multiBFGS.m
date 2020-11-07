function[P,summary]=multiBFGS(X,Y,P_ini,defo,objfun,options)
    if isfield(defo,'action') == 0
        defo.action='pushforward';
    end

    global X_1 Y_1 objfun_1 defo_1
    X_1 = X; Y_1 = Y; 
    objfun_1 = objfun; 
    defo_1 = defo;
 
    [n,d] = size(P_ini.center);
    [temp,~] = size(objfun.kernel_size_geom);
 
    if nargin<6
        options.record_history =true;
        options.nvec = 10;
        options.prtlevel = 2;
    end
  
    if strcmp(defo_1.method,'matlab')
        defo_1.Kq = Kqop(X_1,defo_1);
    end
  
    for i=1:temp
        %  gopt.maxiteration = opt.maxiteration(i,:);
        objfun_1.kernel_size_geom = objfun.kernel_size_geom(i,:);
        objfun_1.kernel_size_grass = objfun.kernel_size_grass(i,:);

        p_0 = stru2vec(P_ini);
        [p,summary{i}] = perform_bfgs('obj_cost', p_0, options);

        P_ini = vec2stru(p,n,d);

    end

    P = P_ini;

end

