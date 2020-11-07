function [P_op,summary]= registration(X,Y,defo,objfun,options)
    if isfield(defo,'action') == 0
        defo.action='pushforward';
    end
    
    P.center = zeros(size(X.center));
    if isfield(X,'vector') ~= 0
        for i=1:length(X.vector)
            P.vector{i} = zeros(size(X.center));
        end
    end
    
    [P_op,summary] = multiBFGS(X,Y,P,defo,objfun,options);

end