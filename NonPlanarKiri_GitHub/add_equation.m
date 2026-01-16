function M = add_equation(M, eq, i, j, Ri, Rj, s)
    rows = (eq-1)*3 + (1:3);

    col_i = (i-1)*3 + (1:3);
    col_j = (j-1)*3 + (1:3);

    M(rows, col_i) =  Ri * skewsym(s);
    M(rows, col_j) = -Rj *skewsym( s);
end
