function M = add_equation_2(M, eq, i, j, p, q, Ri, Rj, Rp, Rq, si, sj, sp, sq)
    rows = (eq-1)*3 + (1:3);

    col_i = (i-1)*3 + (1:3);
    col_j = (j-1)*3 + (1:3);
    col_p = (p-1)*3 + (1:3);
    col_q = (q-1)*3 + (1:3);

    M(rows, col_i) =  Ri *skewsym( si);
    M(rows, col_j) =  Rj *skewsym( sj);
    M(rows, col_p) =  Rp *skewsym( sp);
    M(rows, col_q) =  Rq * skewsym(sq);
end
