function gini = gini_dist(Dist,Grid)
% GINI_DIST computes gini coefficients from discrete distribution function
GridPts = length(Grid);
for i=1:GridPts
    if i==1
        S(i) = Dist(i)*Grid(i);
    else
        S(i) = S(i-1) + Dist(i)*Grid(i);
    end
end

Numerator = Dist(1)*(0+S(1));
for i=2:GridPts
    Numerator = Numerator + Dist(i)*(S(i-1)+S(i));
end
gini = 1-Numerator/S(GridPts);
end