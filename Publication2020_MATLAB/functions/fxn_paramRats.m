%%% Translates Sobol parameters to actual parameters

function [alpha, beta, eta, theta] = fxn_paramRats(beta1, thetabeta, etaalpha, DELTA);

beta = beta1;

theta = beta*thetabeta;

alpha = beta * thetabeta * DELTA / (thetabeta - etaalpha);

eta = etaalpha * alpha;

end