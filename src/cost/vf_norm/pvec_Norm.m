function G = pvec_Norm(XiX)
%Compute the norm of p-vectors

[F_X,~] = size(XiX{1});
m=length(XiX);

A=zeros(m,m,F_X);
for i=1:m
    for j=1:m
        A(i,j,:)=sum(XiX{i}.*XiX{j},2);
    end
end

if m==1
    G=squeeze(A);
    
elseif m==2
    G=squeeze(A(1,1,:).*A(2,2,:)-A(1,2,:).*A(2,1,:));
    
else
    G=zeros(F_X,1);
    for k=1:F_X
        G(k)=det(A(:,:,k));
    end
    
end
G=sqrt(G);
end