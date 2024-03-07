function [Delta] = poly2D_discriminant(a)
%POLY2D_DISCRIMNINANT This function takes as input a monic polynomial of degree 2
%it returns as output the polynomial discriminant in the reduced form
Delta = a(2)^2/4-a(3);
end

