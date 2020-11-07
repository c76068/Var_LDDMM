function IN = Inn(S,T)
[~,d] = size(S);
IN = repmat(sum(S.*T,2),1,d);
end
