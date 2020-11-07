function PJ = PJ(S,T)
%project S to orthogonal complement of T

T = T./sqrt(Inn(T,T)); %normalize T
PJ=S-Inn(S,T).*T;
end