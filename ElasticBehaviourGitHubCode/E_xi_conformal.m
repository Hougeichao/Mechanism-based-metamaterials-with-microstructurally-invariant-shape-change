%% Case-1 Auxetic design energy: pt & added & adjacent
% ccc;
load U_1.mat
load geo_1.mat
U = U_1;
D = [2 0;0. 2];
% d1 = 0.1;d2 = 0.2; D = 1/(1+4*d1*d2)*[2 4*d2;-4*d1 2];
xi_vals_1 = linspace(-pi/6, pi/6, 200);
E_pt_1 = zeros(size(xi_vals_1));
E_add_1 = zeros(size(xi_vals_1));
E_adj_1 = zeros(size(xi_vals_1));
num_spr = length(U);
for k = 1:length(xi_vals_1)
    E_pt_1(k) = 1/num_spr*elastic_energy_pt(xi_vals_1(k), U, D);
end

for k = 1:length(xi_vals_1)
    [E_add_temp,E_adj_temp] = elastic_energy_case1(xi_vals_1(k), geo_1, D, 0);
    E_add_1(k) = E_add_temp;
    E_adj_1(k) = E_adj_temp;

end

figure(); hold on;
plot(xi_vals_1, E_pt_1, 'LineWidth', 2, 'Color',[0 0 0 ]);
plot(xi_vals_1, E_add_1, 'LineWidth', 2, 'Color',[1 0 0]);
plot(xi_vals_1, E_adj_1, 'LineWidth', 2, 'Color',[0 0 1]);


xlabel('\xi', 'Interpreter', 'tex');
ylabel('E(\xi)', 'Interpreter', 'tex');
title('Elastic Energy E(\xi)');
grid on;
% xlim([-1 1]);
% ylim([-0.1 20]);
legend('$E_{\text{pt}}$','$E_{\text{add}}$','$E_{\text{adj}}$','Interpreter', 'latex');
set(gca, 'FontSize', 24);

% [~,~] = elastic_energy_case1(0, geo_1, D, 1);

%% Case-2 Energy: pt & added & adjacent
% ccc;
load U_2.mat
load geo_2.mat
U = U_2;
D = [2 0;0. 2];
xi_vals_2 = linspace(-pi/3, pi/3, 200);
E_pt_2 = zeros(size(xi_vals_2));
E_add_2 = zeros(size(xi_vals_2));
E_adj_2 = zeros(size(xi_vals_2));
num_spr = length(U);
for k = 1:length(xi_vals_2)
    E_pt_2(k) = 1/num_spr*elastic_energy_pt(xi_vals_2(k), U, D);
end

for k = 1:length(xi_vals_2)
    [E_add_temp,E_adj_temp] = elastic_energy_case2(xi_vals_2(k), geo_2);
    E_add_2(k) = E_add_temp;
    E_adj_2(k) = E_adj_temp;

end

figure(); hold on;
plot(xi_vals_2, E_pt_2, 'LineWidth', 2, 'Color',[0 0 0 ]);
plot(xi_vals_2, E_add_2, 'LineWidth', 2, 'Color',[1 0 0]);
plot(xi_vals_2, E_adj_2, 'LineWidth', 2, 'Color',[0 0 1]);


xlabel('\xi', 'Interpreter', 'tex');
ylabel('E(\xi)', 'Interpreter', 'tex');
title('Elastic Energy E(\xi)');
grid on;
% xlim([-1 1]);
% ylim([-0.1 20]);
legend('$E_{\text{pt}}$','$E_{\text{add}}$','$E_{\text{adj}}$','Interpreter', 'latex');
set(gca, 'FontSize', 24);



%% Case-3 Energy: added & adjacent
% ccc;
load U_3.mat
load geo_3.mat
U = U_3;
D = [2 0;0. 2];
xi_vals_3 = linspace(-pi/3, pi/3, 200);
E_pt_3 = zeros(size(xi_vals_3));
E_add_3 = zeros(size(xi_vals_3));
E_adj_3 = zeros(size(xi_vals_3));
num_spr = length(U);
for k = 1:length(xi_vals_3)
    E_pt_3(k) = 1/num_spr*elastic_energy_pt(xi_vals_3(k), U, D);
end

for k = 1:length(xi_vals_3)
    [E_add_temp,E_adj_temp] = elastic_energy_case3(xi_vals_3(k), geo_3);
    E_add_3(k) = E_add_temp;
    E_adj_3(k) = E_adj_temp;

end

figure(); hold on;
plot(xi_vals_3, E_pt_3, 'LineWidth', 2, 'Color',[0 0 0 ]);
plot(xi_vals_3, E_add_3, 'LineWidth', 2, 'Color',[1 0 0]);
plot(xi_vals_3, E_adj_3, 'LineWidth', 2, 'Color',[0 0 1]);


xlabel('\xi', 'Interpreter', 'tex');
ylabel('E(\xi)', 'Interpreter', 'tex');
title('Elastic Energy E(\xi)');
grid on;
% xlim([-1 1]);
% ylim([-0.1 20]);
legend('$E_{\text{pt}}$','$E_{\text{add}}$','$E_{\text{adj}}$','Interpreter', 'latex');
set(gca, 'FontSize', 24);




%%
% E_adj_1 = E_adj_1 / max(E_adj_1);
% % E_adj_2 = E_adj_2 / max(E_adj_2);
% E_adj_3 = E_adj_3 / max(E_adj_3);



figure(); hold on;box on
plot(xi_vals_1, E_adj_1, 'k',  'LineWidth', 2);
plot(xi_vals_2, E_adj_2, 'r', 'LineWidth', 2);
plot(xi_vals_3, E_adj_3, 'b', 'LineWidth', 2);
xlabel('$\xi$', 'Interpreter', 'latex', 'FontSize', 16);
% ylabel('$E_{\mathrm{adj}}$', 'Interpreter', 'latex', 'FontSize', 16);
ylabel('$E_{\mathrm{adj}} / E_{\max}$','Interpreter','latex','FontSize',16);

title('Asymmetric (ReLU) Spring Energy for Adjacent Panels', ...
      'Interpreter', 'latex', 'FontSize', 16);

legend({'$E_{\mathrm{adj}}^{(1)}$', ...
        '$E_{\mathrm{adj}}^{(2)}$', ...
        '$E_{\mathrm{adj}}^{(3)}$'}, ...
        'Interpreter', 'latex', ...
        'FontSize', 14, ...
        'Location', 'northwest');

set(gca, 'FontSize', 14);
grid on;
ax = gca;
ax.FontSize = 16;

% Force scientific notation
ax.YAxis.Exponent = -5;   % or -6 depending on range
ylim([0 1.2e-5]);         % choose tight limits

%%
figure; hold on;

% --- Left axis: small energy ---
yyaxis left

plot(xi_vals_3, E_adj_3, 'Color',[0 176 80]/255, 'LineWidth', 5);
ylabel('$E_{3}$', 'Interpreter', 'latex');
% ylim([0 4.2e-5])
ylim([0 8.4e-4])
set(gca,'YColor',[0 176 80]/255)

yyaxis right

plot(xi_vals_2, E_adj_2, ...
     'k', 'LineWidth', 5)
plot(xi_vals_1, E_adj_1,'Color',[255 140 0]/255, 'LineWidth', 5);
ylabel('$E_{1}$', 'Interpreter', 'latex');
% ylim([0 4.2e-7])
ylim([0 6.8e-6])
set(gca,'YColor',[255 140 0]/255, 'FontSize', 42);
xlabel('$\xi$', 'Interpreter', 'latex');
xlim([-pi/4 pi/4])
% grid on

legend({'$E_{1}$', ...
        '$E_{2}$', ...
        '$E_{3}$'}, ...
        'Interpreter', 'latex', ...
        'FontSize', 24, ...
        'Location', 'northwest');
pbaspect([5.5 6 1])
