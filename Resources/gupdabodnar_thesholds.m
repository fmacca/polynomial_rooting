function [threshold] = gupdabodnar_thesholds(p,alpha)
% Thresholds of acceptance of the Gupda-Bodnar test from their paper

if alpha==0.05
    switch p
        case 2
            threshold=4.22;
        case 4
            threshold=10.25;
        case 6
            threshold=15.01;
        case 8
            threshold=19.34;
        case 10
            threshold=23.15;
        case 12
            threshold=26.84;
        case 14
            threshold=30.60;
        otherwise
            threshold=Inf;
            disp('Threshold value out of table!')
    end
else
    threshold=Inf;
    disp('alpha out of table!')
end


end