function x = solveUsingLU(A, b, itr, verbose)
    M = size(A, 1);
    [L, U] = factorizeLU(A, M);
    Q = L + U - eye(M);
    myfprintf(verbose, "LU Factors for this iteration %i are:", itr);
    myfprintf(verbose, "L: ");
    mydisplay(verbose, L);
    myfprintf(verbose, "U: ");
    mydisplay(verbose, U);
    myfprintf(verbose, "Q: ");
    mydisplay(verbose, Q);
    y = solveUsingForwardSubstitution(L, b, M);
    myfprintf(verbose, "y in Uy = b for iteration %i is:", itr);
    mydisplay(verbose, y);
    x = solveUsingBackwardSubstitution(U, y, M);
    myfprintf(verbose, "x in Lx = y for iteration %i is:", itr);
    mydisplay(verbose, x);
end