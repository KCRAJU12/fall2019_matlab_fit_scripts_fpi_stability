
% simulate some data
X = linspace(0,100,200);
Y = 20./((X-30).^2+20)+0.08*randn(size(X));

% rough guess of initial parameters
a3 = ((max(X)-min(X))/10)^2;
a2 = (max(X)+min(X))/2;
a1 = max(Y)*a3;
a0 = [a1,a2,a3];

% define lorentz inline, instead of in a separate file
lorentz = @(param, x) param(1) ./ ((x-param(2)).^2 + param(3));

% define objective function, this captures X and Y
fit_error = @(param) sum((Y - lorentz(param, X)).^2);

% do the fit
a_fit = fminsearch(fit_error, a0);

% quick plot
x_grid = linspace(min(X), max(X), 1000); % fine grid for interpolation
plot(X, Y, '.', x_grid, lorentz(a_fit, x_grid), 'r')
legend('Measurement', 'Fit')
title(sprintf('a1_fit = %g, a2_fit = %g, a3_fit = %g', ...
    a_fit(1), a_fit(2), a_fit(3)), 'interpreter', 'none')