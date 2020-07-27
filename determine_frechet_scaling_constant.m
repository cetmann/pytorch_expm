function ell_m = determine_frechet_scaling_constant(m, precision)
    % Determines the cutoff values for matrix 1-norms, for which the
    % [m/m]-Pade approximation of the Frechet derivative of the matrix
    % exponential is lower than the given floating point precision.
    % The method is described in
    % Al-Mohy, Awad H., and Nicholas J. Higham. "Computing the Fréchet 
    % derivative of the matrix exponential, with an application to 
    % condition number estimation." SIAM Journal on Matrix Analysis and 
    % Applications 30.4 (2009): 1639-1657.
    % In the above publication, only the values for double precision (and 
    % only up to three significant values) are provided, cf. Table 6.1.
    
    if nargin < 2
        precision = 'single'
    end
    
    
    order = 150;
    syms x;
    f = 0;h=0;
    for k=0:m
       k_ = sym(k);
       coeff = nchoosek(m,k_) * factorial(2*m-k)/factorial(2*m);
       f = f + coeff*x.^k_;
       h = h + coeff*(-x).^k_;
    end
    g = log(exp(-x)*f/h);
    s = taylor(g,'ExpansionPoint',0, 'Order', order);
    
    %%
    c_ = fliplr(coeffs(s, 'All'));
    c = abs(c_);
    g_tilde = 0;
    for k=(2*m+1):(order-1)
       k_ = sym(k);
       coeff = c(k_+1);
       g_tilde = g_tilde + coeff*x.^k_;
    end

    dg_tilde = diff(g_tilde);

    digits(250);
    
    if precision == 'single'
        y = root(dg_tilde-2.^(-24));
    elseif precision == 'double'
        y = root(dg_tilde-2.^(-53));
    end

    y_= vpa(y);
    ell_m = double(y_(1));
end