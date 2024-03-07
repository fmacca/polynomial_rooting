function [discr,cnd] = discriminant_cond(a)
%DISCRIMINANT_COND The function returns the conditional number of the discriminant
%matrix
N=length(a)-1;
D=[convmtx(a',N-1);convmtx(polyder(a),N)]; %Discriminant matrix
cnd=cond(D); %Conditional number of discriminant matrix
discr=det(D); %Discriminant

end

