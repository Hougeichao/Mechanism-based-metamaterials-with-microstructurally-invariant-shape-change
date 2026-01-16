function[R] = Rot3D(theta,r)

ux = r(1);
uy = r(2);
uz = r(3);
C = cos(theta);
S = sin(theta);
t = 1 - cos(theta);

R = [t*ux^2 + C,t*ux*uy - S*uz, t*ux*uz + S*uy;
    t*ux*uy + S*uz, t*uy^2 + C, t*uy*uz - S*ux;
    t*ux*uz - S*uy, t*uy*uz + S*ux, t*uz^2 + C];



% function[R] = Rot(a,b,c)
% 
% R = [cos(a)*cos(b),cos(a)*sin(b)*sin(c)-sin(a)*cos(c),cos(a)*sin(b)*cos(c)+sin(a)*sin(c)
%     ; sin(a)*cos(b),sin(a)*sin(b)*sin(c)+cos(a)*cos(c),sin(a)*sin(b)*cos(c)-cos(a)*sin(c)
%     ;-sin(b),cos(b)*sin(c),cos(b)*cos(c)];
% 
% end

