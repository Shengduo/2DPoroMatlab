function [BC, fval] = findBCKeepKappaCmass(nuu, nuu0, nu0, c0, B0, G)
    % fun = @(BC)calDiffCmassKappa(BC, nuu, nuu0, c0, nu0, B0, G);
    % diff = fun([B0, c0]);
    % disp("diff = " + num2str(diff));
    x0 = [B0, c0];
    options = optimset('TolX', 1e-14, 'TolFun', 1e-20);
    [BC, fval] = fminsearch(@(BC)calDiffCmassKappa(BC, nuu, nuu0, c0, nu0, B0, G), x0, options);
    tgt = calCmassKappa(nuu0, nu0, c0, B0, G);
    this = calCmassKappa(nuu, nu0, BC(2), BC(1), G);
    disp("Target cmass, kappa: " + num2str(tgt(1)) + ", " + num2str(tgt(2)));
    disp("New cmass, kappa: " + num2str(this(1)) + ", " + num2str(this(2)));
end

function res = calCmassKappa(nuu0, nu0, c0, B0, G)
    alpB0 = 3 * (nuu0 - nu0) / (B0 * (1 + nuu0) * (1 - 2 * nu0));
    kappa0 = c0 / (2 * G * (1 + nu0) * B0 / (3 * alpB0 * (1 - alpB0 * B0) * (1 - 2 * nu0)));
    M0 = (2 * G * (nuu0 - nu0)) / (alpB0 ^ 2 * (1 - 2 * nuu0) * (1 - 2 * nu0));
    K0 = M0 * alpB0 * (1 - alpB0) / B0;
    Ku0 = M0 * alpB0 / B0;
    cmass0 = c0 * (K0 + 4. / 3. * G) / (Ku0 + 4. / 3. * G);
    res = [cmass0, kappa0];
end

function diff = calDiffCmassKappa(BC, nuu, nuu0, c0, nu0, B0, G)
    B = BC(1); 
    c = BC(2);
    res0 = calCmassKappa(nuu0, nu0, c0, B0, G);
    res1 = calCmassKappa(nuu, nu0, c, B, G);
    diff = norm(res1 ./ res0 - 1., 2) ^ 2;
end