function [E_add,E_adj] = elastic_energy_case1(xi, geo_1, D, p)

a_1 = geo_1(1,:)';a_2 = geo_1(2,:)';a_3 = geo_1(3,:)';
b_1 = geo_1(4,:)';b_2 = geo_1(5,:)';b_3 = geo_1(6,:)';b_4 = geo_1(7,:)';
c_1 = geo_1(8,:)';c_2 = geo_1(9,:)';c_3 = geo_1(10,:)';c_4 = geo_1(11,:)';
d_1 = geo_1(12,:)';d_2 = geo_1(13,:)';d_3 = geo_1(14,:)';d_4 = geo_1(15,:)';d_5 = geo_1(16,:)';

ac_1 = geo_1(17,:)';ac_2 = geo_1(18,:)';ac_3 = geo_1(19,:)';
bc_1 = geo_1(20,:)';bc_2 = geo_1(21,:)';bc_3 = geo_1(22,:)';bc_4 = geo_1(23,:)';
cc_1 = geo_1(24,:)';cc_2 = geo_1(25,:)';cc_3 = geo_1(26,:)';cc_4 = geo_1(27,:)';
dc_1 = geo_1(28,:)';dc_2 = geo_1(29,:)';dc_3 = geo_1(30,:)';dc_4 = geo_1(31,:)';dc_5 = geo_1(32,:)';

a_1p = geo_1(33,:)';a_2p = geo_1(34,:)';a_3p = geo_1(35,:)';a_4p = geo_1(36,:)';
b_1p = geo_1(37,:)';b_2p = geo_1(38,:)';b_3p = geo_1(39,:)';
c_1p = geo_1(40,:)';c_2p = geo_1(41,:)';c_3p = geo_1(42,:)';c_4p = geo_1(43,:)';c_5p = geo_1(44,:)';
d_1p = geo_1(45,:)';d_2p = geo_1(46,:)';d_3p = geo_1(47,:)';d_4p = geo_1(48,:)';

Rm = Rot(-xi);
Rp = Rot(xi);

% E_add = 1/16*( (norm(Rp*(a_4p-b_2p)+Rm*b_1)-norm(a_4p-b_2p+b_1))^2+(norm(Rp*(a_1p-d_1p)+Rm*c_1)-norm(a_1p-d_1p+c_1))^2+(norm(Rp*(a_2p-c_4p)+Rm*d_1)-norm(a_2p-c_4p+d_1))^2+(norm(Rp*(a_3p-d_3p)+Rm*a_1)-norm(a_3p-d_3p+a_1))^2+     ...
%                (norm(Rp*(b_1p-a_1p)+Rm*c_4)-norm(b_1p-a_1p+c_4))^2+(norm(Rp*(b_2p-c_5p)+Rm*b_2)-norm(b_2p-c_5p+b_2))^2+(norm(Rp*(b_3p-c_2p)+Rm*d_3)-norm(b_3p-c_2p+d_3))^2+   ...
%                (norm(Rp*(c_1p-b_1p)+Rm*c_3)-norm(c_1p-b_1p+c_3))^2+(norm(Rp*(c_2p-d_4p)+Rm*d_4)-norm(c_2p-d_4p+d_4))^2+(norm(Rp*(c_3p-a_3p)+Rm*a_3)-norm(c_3p-a_3p+a_3))^2+(norm(Rp*(c_4p-b_3p)+Rm*d_2)-norm(c_4p-b_3p+d_2))^2+(norm(Rp*(c_5p-d_2p)+Rm*b_3)-norm(c_5p-d_2p+b_3))^2+ ...
%                (norm(Rp*(d_4p-a_2p)+Rm*d_5)-norm(d_4p-a_2p+d_5))^2+(norm(Rp*(d_1p-c_1p)+Rm*c_2)-norm(d_1p-c_1p+c_2))^2+(norm(Rp*(d_2p-a_4p)+Rm*b_4)-norm(d_2p-a_4p+b_4))^2+(norm(Rp*(d_3p-c_3p)+Rm*a_2)-norm(d_3p-c_3p+a_2))^2 ...
%     );
% 
% E_adj = 1/16*( (norm(Rm*ac_2-Rp*c_3p)-norm(ac_2-c_3p))^2+(norm(Rm*ac_3-Rp*d_3p)-norm(ac_3-d_3p))^2+(norm(Rm*ac_1-Rp*a_3p)-norm(ac_1-a_3p))^2+ ...
%                (norm(Rm*bc_1-Rp*a_4p)-norm(bc_1-a_4p))^2+(norm(Rm*bc_2-Rp*d_2p)-norm(bc_2-d_2p))^2+(norm(Rm*bc_3-Rp*c_5p)-norm(bc_3-c_5p))^2+(norm(Rm*bc_4-Rp*b_2p)-norm(bc_4-b_2p))^2+ ...
%                (norm(Rm*cc_1-Rp*a_1p)-norm(cc_1-a_1p))^2+(norm(Rm*cc_2-Rp*b_1p)-norm(cc_2-b_1p))^2+(norm(Rm*cc_3-Rp*c_1p)-norm(cc_3-c_1p))^2+(norm(Rm*cc_4-Rp*d_1p)-norm(cc_4-d_1p))^2+ ...
%                (norm(Rm*dc_1-Rp*a_2p)-norm(dc_1-a_2p))^2+(norm(Rm*dc_2-Rp*d_4p)-norm(dc_2-d_4p))^2+(norm(Rm*dc_3-Rp*c_2p)-norm(dc_3-c_2p))^2+(norm(Rm*dc_4-Rp*b_3p)-norm(dc_4-b_3p))^2+(norm(Rm*dc_5-Rp*c_4p)-norm(dc_5-c_4p))^2);

%% asymmetric spring
E_add =(spring(b_1,a_4p-b_2p,Rm,Rp)+spring(c_1,a_1p-d_1p,Rm,Rp)+spring(d_1,a_2p-c_4p,Rm,Rp)+spring(a_1,a_3p-d_3p,Rm,Rp)+ ...
              spring(c_4,b_1p-a_1p,Rm,Rp)+spring(b_2,b_2p-c_5p,Rm,Rp)+spring(d_3,b_3p-c_2p,Rm,Rp)+ ...
              spring(c_3,c_1p-b_1p,Rm,Rp)+spring(d_4,c_2p-d_4p,Rm,Rp)+spring(a_3,c_3p-a_3p,Rm,Rp)+spring(d_2,c_4p-b_3p,Rm,Rp)+spring(b_3,c_5p-d_2p,Rm,Rp)+ ...
              spring(d_5,d_4p-a_2p,Rm,Rp)+spring(c_2,d_1p-c_1p,Rm,Rp)+spring(b_4,d_2p-a_4p,Rm,Rp)+spring(a_2,d_3p-c_3p,Rm,Rp));

E_adj = (spring(ac_2,-c_3p,Rm,Rp)+spring(ac_3,-d_3p,Rm,Rp)+spring(ac_1,-a_3p,Rm,Rp)+spring(bc_1,-a_4p,Rm,Rp)+spring(bc_2,-d_2p,Rm,Rp)+spring(bc_3,-c_5p,Rm,Rp)+spring(bc_4,-b_2p,Rm,Rp)+ ...
              spring(cc_1,-a_1p,Rm,Rp)+spring(cc_2,-b_1p,Rm,Rp)+spring(cc_3,-c_1p,Rm,Rp)+spring(cc_4,-d_1p,Rm,Rp)+spring(dc_1,-a_2p,Rm,Rp)+spring(dc_2,-d_4p,Rm,Rp)+spring(dc_3,-c_2p,Rm,Rp)+spring(dc_4,-b_3p,Rm,Rp)+spring(dc_5,-c_4p,Rm,Rp));


%%
if p == 1
    Ident = Rot(0);
    P_a = [ [0;0], a_1, a_1+a_2 ];
    P_b = [ [0;0], b_1, b_1+b_2, b_1+b_2+b_3 ];
    P_c = [ [0;0], c_1, c_1+c_2, c_1+c_2+c_3 ];
    P_d = [ [0;0], d_1, d_1+d_2, d_1+d_2+d_3, d_1+d_2+d_3+d_4 ];

    T_a = (D - Ident) * (-ac_1);
    T_b = (D - Ident) * (-bc_1);
    T_c = (D - Ident) * (-cc_1);
    T_d = (D - Ident) * (-dc_1);
    P_a_def = P_a + T_a;
    P_b_def = P_b + T_b;
    P_c_def = P_c + T_c;
    P_d_def = P_d + T_d;


    figure();hold on;axis equal;axis off;
    patch(P_a(1,:), P_a(2,:), [200/255 200/255 200/255], 'FaceAlpha', 0.78, 'LineWidth',1.2);
    patch(P_b(1,:), P_b(2,:), [200/255 200/255 200/255], 'FaceAlpha', 0.78, 'LineWidth',1.2);
    patch(P_c(1,:), P_c(2,:), [200/255 200/255 200/255], 'FaceAlpha', 0.78, 'LineWidth',1.2);
    patch(P_d(1,:), P_d(2,:), [200/255 200/255 200/255], 'FaceAlpha', 0.78, 'LineWidth',1.2);

    figure();hold on;axis equal;axis off;
    patch(P_a_def(1,:), P_a_def(2,:), [200/255 200/255 200/255], 'FaceAlpha', 0.78, 'LineWidth',1.2);
    patch(P_b_def(1,:), P_b_def(2,:), [200/255 200/255 200/255], 'FaceAlpha', 0.78, 'LineWidth',1.2);
    patch(P_c_def(1,:), P_c_def(2,:), [200/255 200/255 200/255], 'FaceAlpha', 0.78, 'LineWidth',1.2);
    patch(P_d_def(1,:), P_d_def(2,:), [200/255 200/255 200/255], 'FaceAlpha', 0.78, 'LineWidth',1.2);


    xlabel('x'); ylabel('y');
end



end
