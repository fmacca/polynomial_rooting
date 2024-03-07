function [a_transf,t] = poly2D_affinetransf(a)
%POLY2D_AFFINETRANSF This function takes as input a polynomial of degree 2
%it returns as output the polynomial where the term of degree 2 has been
%"eliminated" via an affine transformation z = z' + t
%In addition to that, it returns the shift term t
t=-a(2)/2;
a_transf=[1; 0; -a(2)^2/4+a(3)];
end

