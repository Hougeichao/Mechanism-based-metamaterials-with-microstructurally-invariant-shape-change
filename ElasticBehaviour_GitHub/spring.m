function E = spring(m,p,Rm,Rp)

    L_def = norm(Rm*m + Rp*p);
    L_nat = norm(m + p);
    
    if L_def > L_nat
        E = (L_def - L_nat)^2;
    else
        E = 0;
        % E = (L_def - L_nat)^2;
    end
    
end