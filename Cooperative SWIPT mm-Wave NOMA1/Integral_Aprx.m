function out = Integral_Aprx(mu,chi,xi)
%
A11 = exp(-mu*chi)/chi;
A22 = xi*igamma(0,mu*chi);
A33 = 0;
for uu=2:1000
    B11 = ((-1)^uu)*(xi^uu)/(factorial(uu));
    B211 = exp(-(mu*chi));
    B222 = 0;
    for vv = 1:(uu-1)
        temp = (factorial(vv-1))*((-chi)^(uu-vv-1))/...
            ((factorial(uu-1))*(mu^vv));
        B222 = B222 + temp;
    end
    B233 = ((-chi)^(uu-1))/(factorial(uu-1))*(ei(-mu*chi));
    B22 = B211*B222-B233;
    A33 = A33 + B11*B22;
    
end
out = A11-A22+A33;