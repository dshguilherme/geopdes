function s_xx = test_plate_with_hole_g_xy(x, y)

% Inner radius and applied stress
T = 10;
R = 1;
[theta, r] =  cart2pol(x,y);
theta = theta +pi;


s_xx = T*( 1 -((R./r).^2).*(1.5*cos(2*theta) ...
    +cos(4*theta)) +1.5*((R./r).^4).*cos(4*theta));


end