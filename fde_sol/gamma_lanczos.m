function y = gamma_lanczos(x)
    g = 7;
    coeffs = [
    0.99999999999980993,
    676.5203681218851,
   -1259.1392167224028,
    771.32342877765313,
   -176.61502916214059,
    12.507343278686905,
   -0.13857109526572012,
    9.9843695780195716e-6,
    1.5056327351493116e-7
    ];

    if x<0 && mod(x,1) == 0
        y = inf;
        return
    end
    
    if x < 0.5
            %reflection formula
           y = pi / ((sin(pi*x)) * gamma_lanczos(1-x));
           return
    end

    x = x-1;
    a = coeffs(1);
    for i = 2:length(coeffs)
        a = a + coeffs(i) / (x + i - 1);
    end

    t = x + g + 0.5;
    y = sqrt(2 * pi) * t^(x + 0.5) * exp(-t) * a;
end