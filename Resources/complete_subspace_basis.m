function [B] = complete_subspace_basis(b1)
%COMPLETE_SUBSPACE_BASIS Summary of this function goes here
%   Detailed explanation goes here
N=length(b1);
B=zeros(N,N-2);
jj=1;
for ii=1:N
    if b1(ii)==0
        B(ii,jj)=1;
        b1(ii)=1;
        jj=jj+1;
    end
end
% if sum(b1)==N
%     disp("Basis completed correctly");
% end

end

