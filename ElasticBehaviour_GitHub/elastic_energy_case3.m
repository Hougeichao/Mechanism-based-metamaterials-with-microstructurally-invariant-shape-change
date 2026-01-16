function [E_add,E_adj] = elastic_energy_case3(xi, geo_3)

a_1 = geo_3(1,:)';a_2 = geo_3(2,:)';a_3 = geo_3(3,:)';a_4 = geo_3(4,:)';
b_1 = geo_3(5,:)';b_2 = geo_3(6,:)';b_3 = geo_3(7,:)';b_4 = geo_3(8,:)';
c_1 = geo_3(9,:)';c_2 = geo_3(10,:)';c_3 = geo_3(11,:)';c_4 = geo_3(12,:)';
d_1 = geo_3(13,:)';d_2 = geo_3(14,:)';d_3 = geo_3(15,:)';d_4 = geo_3(16,:)';
e_1 = geo_3(17,:)';e_2 = geo_3(18,:)';e_3 = geo_3(19,:)';e_4 = geo_3(20,:)';

ac_1 = geo_3(21,:)';ac_2 = geo_3(22,:)';ac_3 = geo_3(23,:)';ac_4 = geo_3(24,:)';
bc_1 = geo_3(25,:)';bc_2 = geo_3(26,:)';bc_3 = geo_3(27,:)';bc_4 = geo_3(28,:)';
cc_1 = geo_3(29,:)';cc_2 = geo_3(30,:)';cc_3 = geo_3(31,:)';cc_4 = geo_3(32,:)';
dc_1 = geo_3(33,:)';dc_2 = geo_3(34,:)';dc_3 = geo_3(35,:)';dc_4 = geo_3(36,:)';
ec_1 = geo_3(37,:)';ec_2 = geo_3(38,:)';ec_3 = geo_3(39,:)';ec_4 = geo_3(40,:)';

a_1p = geo_3(41,:)';a_2p = geo_3(42,:)';a_3p = geo_3(43,:)';
b_1p = geo_3(44,:)';b_2p = geo_3(45,:)';b_3p = geo_3(46,:)';
c_1p = geo_3(47,:)';c_2p = geo_3(48,:)';c_3p = geo_3(49,:)';
d_1p = geo_3(50,:)';d_2p = geo_3(51,:)';d_3p = geo_3(52,:)';
e_1p = geo_3(53,:)';e_2p = geo_3(54,:)';e_3p = geo_3(55,:)';e_4p = geo_3(56,:)';e_5p = geo_3(57,:)';e_6p = geo_3(58,:)';e_7p = geo_3(59,:)';e_8p = geo_3(60,:)';

Rm = Rot(-xi);
Rp = Rot(xi);


% E_add = 1/20*(spring(e_2,a_1p-b_1p,Rm,Rp)+spring(d_2,a_2p-e_8p,Rm,Rp)+spring(a_1,a_3p-d_2p,Rm,Rp)+spring(e_3,b_1p-c_1p,Rm,Rp)+spring(c_2,b_2p-e_6p,Rm,Rp)+spring(d_1,b_3p-a_2p,Rm,Rp)+ ...
%               spring(e_4,c_1p-d_1p,Rm,Rp)+spring(b_2,c_2p-e_4p,Rm,Rp)+spring(c_1,c_3p-b_2p,Rm,Rp)+spring(e_1,d_1p-a_1p,Rm,Rp)+spring(a_2,d_2p-e_2p,Rm,Rp)+spring(b_1,d_3p-c_2p,Rm,Rp)+ ...
%               spring(b_4,e_1p-d_3p,Rm,Rp)+spring(a_3,e_2p-e_7p,Rm,Rp)+spring(c_4,e_3p-c_3p,Rm,Rp)+spring(b_3,e_4p-e_1p,Rm,Rp)+spring(d_4,e_5p-b_3p,Rm,Rp)+spring(c_3,e_6p-e_3p,Rm,Rp)+spring(a_4,e_7p-a_3p,Rm,Rp)+spring(d_3,e_8p-e_5p,Rm,Rp));
% 
% 
% 
% E_adj = 1/20*(spring(ac_1,-d_2p,Rm,Rp)+spring(ac_2,-a_3p,Rm,Rp)+spring(ac_3,-e_7p,Rm,Rp)+spring(ac_4,-e_2p,Rm,Rp)+spring(bc_1,-c_2p,Rm,Rp)+spring(bc_2,-d_3p,Rm,Rp)+spring(bc_3,-e_1p,Rm,Rp)+spring(bc_4,-e_4p,Rm,Rp)+ ...
%               spring(cc_1,-b_2p,Rm,Rp)+spring(cc_2,-c_3p,Rm,Rp)+spring(cc_3,-e_3p,Rm,Rp)+spring(cc_4,-e_6p,Rm,Rp)+spring(dc_1,-a_2p,Rm,Rp)+spring(dc_2,-b_3p,Rm,Rp)+spring(dc_3,-e_5p,Rm,Rp)+spring(dc_4,-e_8p,Rm,Rp)+ ...
%               spring(ec_1,-c_1p,Rm,Rp)+spring(ec_2,-b_1p,Rm,Rp)+spring(ec_3,-a_1p,Rm,Rp)+spring(ec_4,-d_1p,Rm,Rp));

E_add =(spring(e_2,a_1p-b_1p,Rm,Rp)+spring(d_2,a_2p-e_8p,Rm,Rp)+spring(a_1,a_3p-d_2p,Rm,Rp)+spring(e_3,b_1p-c_1p,Rm,Rp)+spring(c_2,b_2p-e_6p,Rm,Rp)+spring(d_1,b_3p-a_2p,Rm,Rp)+ ...
              spring(e_4,c_1p-d_1p,Rm,Rp)+spring(b_2,c_2p-e_4p,Rm,Rp)+spring(c_1,c_3p-b_2p,Rm,Rp)+spring(e_1,d_1p-a_1p,Rm,Rp)+spring(a_2,d_2p-e_2p,Rm,Rp)+spring(b_1,d_3p-c_2p,Rm,Rp)+ ...
              spring(b_4,e_1p-d_3p,Rm,Rp)+spring(a_3,e_2p-e_7p,Rm,Rp)+spring(c_4,e_3p-c_3p,Rm,Rp)+spring(b_3,e_4p-e_1p,Rm,Rp)+spring(d_4,e_5p-b_3p,Rm,Rp)+spring(c_3,e_6p-e_3p,Rm,Rp)+spring(a_4,e_7p-a_3p,Rm,Rp)+spring(d_3,e_8p-e_5p,Rm,Rp));



E_adj = (spring(ac_1,-d_2p,Rm,Rp)+spring(ac_2,-a_3p,Rm,Rp)+spring(ac_3,-e_7p,Rm,Rp)+spring(ac_4,-e_2p,Rm,Rp)+spring(bc_1,-c_2p,Rm,Rp)+spring(bc_2,-d_3p,Rm,Rp)+spring(bc_3,-e_1p,Rm,Rp)+spring(bc_4,-e_4p,Rm,Rp)+ ...
              spring(cc_1,-b_2p,Rm,Rp)+spring(cc_2,-c_3p,Rm,Rp)+spring(cc_3,-e_3p,Rm,Rp)+spring(cc_4,-e_6p,Rm,Rp)+spring(dc_1,-a_2p,Rm,Rp)+spring(dc_2,-b_3p,Rm,Rp)+spring(dc_3,-e_5p,Rm,Rp)+spring(dc_4,-e_8p,Rm,Rp)+ ...
              spring(ec_1,-c_1p,Rm,Rp)+spring(ec_2,-b_1p,Rm,Rp)+spring(ec_3,-a_1p,Rm,Rp)+spring(ec_4,-d_1p,Rm,Rp));

end
