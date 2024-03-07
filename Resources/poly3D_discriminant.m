function [Delta] = poly3D_discriminant(a)
%POLY2D_DISCRIMNINANT This function takes as input a monic polynomial of
%degree 3
%it returns as output the polynomial discriminant
Delta = a(2)^2*a(3)^2-4*a(3)^3-4*a(2)^3*a(4)-27*a(4)+18*a(2)*a(3)*a(4);
end

