function min_dist = min_dist_roots(r)
%min_dist_roots This function takes as input a
%vector of roots and returns the minimum of the
%pairwise distance between them
N=length(r);
C=nchoosek(1:N,2); %all combinations of N numbers, taken 2 at a time
min_dist=abs(r(C(1,1))-r(C(1,2))); %initialize distance with the first combination
for k=2:size(C,1)
    if abs(r(C(k,1))-r(C(k,2))) < min_dist
        min_dist=abs(r(C(k,1))-r(C(k,2)));%upade min_dist when needed
end

end

