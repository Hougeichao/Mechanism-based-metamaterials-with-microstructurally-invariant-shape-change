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
end
