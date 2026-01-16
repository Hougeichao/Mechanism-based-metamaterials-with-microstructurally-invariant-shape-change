tic; ccc;
e1 = [1;0;0];
e2 = [0;1;0];
e3 = [0;0;1];

% Geometry of the pattern
s1 = e1;
s2 = Rot3D(60/180*pi,e3)*e1;
s3 = Rot3D(120/180*pi,e3)*e1;
s4 = -e1;
s5 = -s2;
s6 = -s3;
s7 = Rot3D(-105/180*pi,e3)*e1;
s8 = Rot3D(15/180*pi,e3)*e1;
s9 = Rot3D(120/180*pi,e3)*s8;
s10 = Rot3D(75/180*pi,e3)*e1;
s11 = -s8;
s12 = Rot3D(120/180*pi,e3)*s11;

Ident = Rot3D(0,e1);
omega_0 = 0.1*pi;
rho_1 = 2*atan((sin(-22.5/180*pi))/(sin(82.5/180*pi))*tan(omega_0/2));
R1_0 = Ident;
R2_0 = Rot3D(omega_0,-s3);
R3_0 = Rot3D(omega_0,-s3)*Rot3D(rho_1,s8);
R4_0 = Rot3D(rho_1,s4);
R5_0 = R3_0;


%%
n_w  = 37;
n_eq_1 = 52;
n_eq_2 = 4;
R = cell(5,1);R{1} = R1_0;R{2} = R2_0;R{3} = R3_0;R{4} = R4_0;R{5} = R5_0;
s_vectors = cell(12,1);s_vectors{1}=s1;s_vectors{2}=s2;s_vectors{3}=s3;s_vectors{4}=s4;s_vectors{5}=s5;s_vectors{6}=s6;s_vectors{7}=s7;s_vectors{8}=s8;s_vectors{9}=s9;s_vectors{10}=s10;s_vectors{11}=s11;s_vectors{12}=s12;

M_1 = zeros(3*n_eq_1, 3*n_w);
M_2 = zeros(3*n_eq_2, 3*n_w);

eqns_1 = [
    1 2  1 2  3;
    5 6  1 2  3;
    9 10 1 2  3;
    13 14  1 2  3;
    18 19  1 2  3;
    23 24  1 2  3;
    27 28  1 2  3;
    32 33  1 2  3;

    2 3  2 3  8;
    6 7  2 3  8;
    10 11 2 3 8;
    14 15 2 3 8;
    19 20 2 3 8;
    24 25 2 3 8;
    28 29 2 3 8;
    33 34 2 3 8;
    
    3 4 3 4 7;
    7 8 3 4 7;
    11 12 3 4 7;
    15 16 3 4 7;
    20 21 3 4 7;
    25 26 3 4 7;
    29 30 3 4 7;
    34 35 3 4 7;

    1 4 1 4 4;
    5 8 1 4 4;
    9 12 1 4 4;
    13 16 1 4 4;
    18 21 1 4 4;
    23 26 1 4 4;
    27 30 1 4 4;
    32 35 1 4 4;

    2 13 2 1 6;
    6 18 2 1 6;
    10 23 2 1 6;
    14 27 2 1 6;
    19 32 2 1 6;
    24 37 2 1 6;

    4 5 4 1 1;
    8 9 4 1 1;
    16 18 4 1 1;
    21 23 4 1 1;
    30 32 4 1 1;
    35 37 4 1 1;

    16 17 4 5 10;
    21 22 4 5 10;
    30 31 4 5 10;
    35 36 4 5 10;
    6 17 2 5 11;
    10 22 2 5 11;
    19 31 2 5 11;
    24 36 2 5 11;


];


for eq = 1:size(eqns_1,1)
    i  = eqns_1(eq,1);
    j  = eqns_1(eq,2);
    Ri = R{eqns_1(eq,3)};
    Rj = R{eqns_1(eq,4)};
    s_temp  = s_vectors{eqns_1(eq,5)};

    M_1 = add_equation(M_1, eq, i, j, Ri, Rj, s_temp);
end


eqns_2 = [
    3 5 13 17 3 1 1 5 9 2 5 12;
    7 9 18 22 3 1 1 5 9 2 5 12;
    15 18 27 31 3 1 1 5 9 2 5 12;
    20 23 32 36 3 1 1 5 9 2 5 12;];

for eq = 1:size(eqns_2,1)
    i  = eqns_2(eq,1);
    j  = eqns_2(eq,2);
    p  = eqns_2(eq,3);
    q  = eqns_2(eq,4);
    Ri = R{eqns_2(eq,5)};
    Rj = R{eqns_2(eq,6)};
    Rp = R{eqns_2(eq,7)};
    Rq = R{eqns_2(eq,8)};
    si  = s_vectors{eqns_2(eq,9)};
    sj  = s_vectors{eqns_2(eq,10)};
    sp  = s_vectors{eqns_2(eq,11)};
    sq  = s_vectors{eqns_2(eq,12)};

    M_2 = add_equation_2(M_2, eq, i, j, p, q, Ri, Rj, Rp, Rq, si, sj, sp, sq);
end

MM = [M_1;M_2];
nM = null(MM);

%% Gram Schmidt process

null1 = repmat([1;0;0], n_w, 1);
null1 = null1/norm(null1);
null2 = repmat([0;1;0], n_w, 1);
null2 = null2/norm(null2);
null3 = repmat([0;0;1], n_w, 1);
null3 = null3/norm(null3);

V=[null1,null2,null3];

u1 = nM(:,1);

w = u1;
for i = 1:size(V, 2)
    vi = V(:, i);
    w = w - (dot(w, vi) / dot(vi, vi)) * vi;
end
v_left = w / norm(w)

%% The left null space vector corresponds to the knwon bulk mecahnism, which justifies the assertion that the pattern has a single DOF

toc


