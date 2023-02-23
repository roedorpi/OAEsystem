function H = CombFilter(b,F) 
% H frequency amplitude of combfilter
% b(1) amplitude of generator
% b(2) amplitude of reflection
% F DP frequency
% b(3) Delay time (1/(fmax - fmin))


H = abs((b(1)+b(2)*cos(2*pi*F*b(3)))...
    -1i*(b(2)*sin(2*pi*F*b(3))));

