function [E_add,E_adj] = elastic_energy_case2(xi, geo_2)

a_1 = geo_2(1,:)';a_2 = geo_2(2,:)';a_3 = geo_2(3,:)';
b_1 = geo_2(4,:)';b_2 = geo_2(5,:)';b_3 = geo_2(6,:)';
c_1 = geo_2(7,:)';c_2 = geo_2(8,:)';c_3 = geo_2(9,:)';
d_1 = geo_2(10,:)';d_2 = geo_2(11,:)';d_3 = geo_2(12,:)';

ac_1 = geo_2(13,:)';ac_2 = geo_2(14,:)';ac_3 = geo_2(15,:)';
bc_1 = geo_2(16,:)';bc_2 = geo_2(17,:)';bc_3 = geo_2(18,:)';
cc_1 = geo_2(19,:)';cc_2 = geo_2(20,:)';cc_3 = geo_2(21,:)';
dc_1 = geo_2(22,:)';dc_2 = geo_2(23,:)';dc_3 = geo_2(24,:)';

a_1p = geo_2(25,:)';a_2p = geo_2(26,:)';a_3p = geo_2(27,:)';a_4p = geo_2(28,:)';
b_1p = geo_2(29,:)';b_2p = geo_2(30,:)';b_3p = geo_2(31,:)';b_4p = geo_2(32,:)';
b_5p = geo_2(33,:)';b_6p = geo_2(34,:)';b_7p = geo_2(35,:)';b_8p = geo_2(36,:)';

Rm = Rot(-xi);
Rp = Rot(xi);

% E_add = 1/12*( (norm(Rp*(a_1p-b_2p)+Rm*a_1)-norm(a_1p-b_2p+a_1))^2+(norm(Rp*(a_2p-b_4p)+Rm*b_1)-norm(a_2p-b_4p+b_1))^2+(norm(Rp*(a_3p-b_6p)+Rm*c_1)-norm(a_3p-b_6p+c_1))^2+(norm(Rp*(a_4p-b_8p)+Rm*d_1)-norm(a_4p-b_8p+d_1))^2+     ...
%                (norm(Rp*(b_1p-a_2p)+Rm*b_3)-norm(b_1p-a_2p+b_3))^2+(norm(Rp*(b_2p-b_7p)+Rm*a_2)-norm(b_2p-b_7p+a_2))^2+(norm(Rp*(b_3p-a_3p)+Rm*c_3)-norm(b_3p-a_3p+c_3))^2+(norm(Rp*(b_4p-b_1p)+Rm*b_2)-norm(b_4p-b_1p+b_2))^2+   ...
%                (norm(Rp*(b_5p-a_4p)+Rm*d_3)-norm(b_5p-a_4p+d_3))^2+(norm(Rp*(b_6p-b_3p)+Rm*c_2)-norm(b_6p-b_3p+c_2))^2+(norm(Rp*(b_7p-a_1p)+Rm*a_3)-norm(b_7p-a_1p+a_3))^2+(norm(Rp*(b_8p-b_5p)+Rm*d_2)-norm(b_8p-b_5p+d_2))^2);
% 
% E_adj = 1/12*( (norm(Rm*ac_1-Rp*a_1p)-norm(ac_1-a_1p))^2+(norm(Rm*bc_1-Rp*a_2p)-norm(bc_1-a_2p))^2+(norm(Rm*cc_1-Rp*a_3p)-norm(cc_1-a_3p))^2+(norm(Rm*dc_1-Rp*a_4p)-norm(dc_1-a_4p))^2+ ...
%                (norm(Rm*bc_2-Rp*b_1p)-norm(bc_2-b_1p))^2+(norm(Rm*ac_3-Rp*b_2p)-norm(ac_3-b_2p))^2+(norm(Rm*cc_2-Rp*b_3p)-norm(cc_2-b_3p))^2+(norm(Rm*bc_3-Rp*b_4p)-norm(bc_3-b_4p))^2+ ...
%                (norm(Rm*dc_2-Rp*b_5p)-norm(dc_2-b_5p))^2+(norm(Rm*cc_3-Rp*b_6p)-norm(cc_3-b_6p))^2+(norm(Rm*ac_2-Rp*b_7p)-norm(ac_2-b_7p))^2+(norm(Rm*dc_3-Rp*b_8p)-norm(dc_3-b_8p))^2);

%% asymmetric spring
E_add = (spring(a_1,a_1p-b_2p,Rm,Rp)+spring(b_1,a_2p-b_4p,Rm,Rp)+spring(c_1,a_3p-b_6p,Rm,Rp)+spring(d_1,a_4p-b_8p,Rm,Rp)+ ...
              spring(b_3,b_1p-a_2p,Rm,Rp)+spring(a_2,b_2p-b_7p,Rm,Rp)+spring(c_3,b_3p-a_3p,Rm,Rp)+spring(b_2,b_4p-b_1p,Rm,Rp)+spring(d_3,b_5p-a_4p,Rm,Rp)+spring(c_2,b_6p-b_3p,Rm,Rp)+spring(a_3,b_7p-a_1p,Rm,Rp)+spring(d_2,b_8p-b_5p,Rm,Rp));



E_adj = (spring(ac_1,-a_1p,Rm,Rp)+spring(bc_1,-a_2p,Rm,Rp)+spring(cc_1,-a_3p,Rm,Rp)+spring(dc_1,-a_4p,Rm,Rp)+spring(bc_2,-b_1p,Rm,Rp)+spring(ac_3,-b_2p,Rm,Rp)+spring(cc_2,-b_3p,Rm,Rp)+spring(bc_3,-b_4p,Rm,Rp)+ ...
              spring(dc_2,-b_5p,Rm,Rp)+spring(cc_3,-b_6p,Rm,Rp)+spring(ac_2,-b_7p,Rm,Rp)+spring(dc_3,-b_8p,Rm,Rp));

end
