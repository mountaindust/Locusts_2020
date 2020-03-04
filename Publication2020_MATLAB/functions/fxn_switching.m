% Locust switching rates, functions of R

function [ksm, kms] = fxn_switching(R, alpha, eta, gamma, beta, theta, delta)
    % pass in:
        % R = row vector, resources for current time step only.
        % greek letters = scalars, switching parameters
    % return: ksm,kms = row vectors, switching rates to/from moving/stationary
    
    ksm = eta-(eta-alpha)*exp(-gamma.*R);
    kms = theta-(theta-beta)*exp(-delta.*R);
end
