function E = elastic_energy_pt(xi, U, D)

Aeff = (Rot(-xi)-Rot(xi))*inv(D)+Rot(xi);
Du = D * U;
AeffDu = Aeff * Du;

Du_norm = vecnorm(Du);
AeffDu_norm = vecnorm(AeffDu);

E = sum((AeffDu_norm - Du_norm).^2);


%% asymmetric spring
diff = AeffDu_norm - Du_norm;     % 1×N
diff(diff < 0) = 0;               % clamp negative values to zero

E = sum(diff.^2);


end
