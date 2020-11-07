function Z = add_XY_h(X,Y,h)
% compute: 
if length(h)<2
    a=1; b=h;
else
    a = h(1); b=h(2);
end

Z.center = a*X.center+ b*Y.center;

if isfield(X,'vector') ~= 0
   for i=1:length(X.vector)
      Z.vector{i} = a*X.vector{i} + b*Y.vector{i};
   end
end

end
