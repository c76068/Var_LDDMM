function p = stru2vec(P)
%convert into a structure into a 3nd x 1 array 

%  if isstruct(P) == 1
%     [n,d] = size(P.center);
%     p = reshape(P.center',n*d,1);
%     for i = 1:length(P.vector)
%       p = [p; reshape(P.vector{i}',n*d,1)];
%     end
% 
%  else
%     [n,d] = size(P);
%     p = reshape(P',n*d,1);
%  end
    
   [n,d] = size(P.center);
   if isfield(P,'vector') == 0
       p = reshape(P.center',n*d,1);    
   else         
      p = reshape(P.center',n*d,1);
      for i = 1:length(P.vector)
         p = [p; reshape(P.vector{i}',n*d,1)];
      end  
   end

end