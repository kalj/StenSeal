function c = notStupidCoeffs(p, vars)
    c = sym(zeros(1,length(vars)));

    [someCoeffs, someVars] = coeffs(p, vars);

    [~, I] = ismember(someVars, vars);
    c(I) = someCoeffs;
end