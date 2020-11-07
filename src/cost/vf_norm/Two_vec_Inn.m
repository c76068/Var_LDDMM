function G = Two_vec_Inn(X,Y)
%Compute inner product of sets of two vectors
%Input:
%   X.vector : a cell contains first  and second vectors of the 2-vector in the first shape.
%   Y.vector : a cell contains first  and second vectors of the 2-vector in the second shape.
%   
% Output
%   G : a real matrix.

[F_X,d] = size(X.vector{1});
[F_Y,~] = size(Y.vector{1});

U = zeros(2*F_X,d); V = zeros(2*F_Y,d);
U(1:2:end,:) =X.vector{1}; U(2:2:end,:) =X.vector{2};
V(1:2:end,:) =Y.vector{1}; V(2:2:end,:) =Y.vector{2};

% C = mat2cell(U*(V'), repmat(2,1,F_X), repmat(2,1,F_Y));
% G = cellfun(@det,C);

A = U*(V');
A_11 = A(1:2:end,1:2:end);
A_12 = A(1:2:end,2:2:end);
A_21 = A(2:2:end,1:2:end);
A_22 = A(2:2:end,2:2:end);
G=A_11.*A_22 -A_12.*A_21;
end