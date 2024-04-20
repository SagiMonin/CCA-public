function plotParabolid(parab,loc_y,loc_x,fig_num)
% a*x^2 + b*x + c + d*y^2 + e * y + f*x*y
[Epx, Epy,deter] = getXYPraboloid(parab);
cord_x = linspace(-3,3,100) + Epx(loc_y,loc_x);
cord_y = linspace(-3,3,100) + Epy(loc_y,loc_x);

[X,Y] = meshgrid(cord_x,cord_y);
z = parab.a(loc_y,loc_x)* X.^2 + parab.b(loc_y,loc_x)* X + parab.c(loc_y,loc_x) +...
    parab.d(loc_y,loc_x)* Y.^2 + parab.e(loc_y,loc_x)* Y + parab.f(loc_y,loc_x) .*X.*Y;

if nargin>3
    figure(fig_num); contour(X,Y,z);
else 
    figure(133); contour(X,Y,z);
end

mat_parab = [parab.a(loc_y,loc_x), parab.f(loc_y,loc_x)/2; ...
            parab.f(loc_y,loc_x)/2, parab.d(loc_y,loc_x)];
        
eig_val = eig(mat_parab);
display(['ratio ', num2str(eig_val(1)/eig_val(2))]);




end

