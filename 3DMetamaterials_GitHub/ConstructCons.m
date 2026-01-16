function M = ConstructCons(M, eq, m, n, p, q, smn, spq, R1_0, R2_0)
    rows = (eq-1)*3 + (1:3);

    col_m = (m-1)*3 + (1:3);
    col_n = (n-1)*3 + (1:3);
    col_p = (p-1)*3 + (1:3);
    col_q = (q-1)*3 + (1:3);

    M(rows, col_m) =  R1_0*skewsym(smn);
    M(rows, col_n) = -R1_0*skewsym(smn);
    M(rows, col_p) =  R2_0*skewsym(spq);
    M(rows, col_q) =  -R2_0*skewsym(spq);
end
