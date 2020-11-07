function P = vec2stru(p,n,d)
%convert a 3nd x 1 array into a structure


        m = length(p)/(n*d) -1;
        
        if m == 0
              P.center = reshape(p(1:n*d),d,n)'; 
        else 
              P.center = reshape(p(1:n*d),d,n)';
              
              for i = 1:m
                P.vector{i} = reshape(p(i*n*d+1:(i+1)*n*d),d,n)';
              end

        end
end