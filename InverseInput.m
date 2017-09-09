function [Omega1, Omega2, Omega3, Omega4] = InverseInput(U1, U2, U3, U4)
Kt = 1.22641e-04; d = 6.48447e-06;    
Omega2_sqr = 0.25*(((U1-2*U2)/Kt)-(U4/d));
Omega1_sqr = Omega2_sqr + 0.5*(((U2-U3)/Kt) + (U4/d));
Omega3_sqr = Omega1_sqr + U3/Kt;
Omega4_sqr = Omega2_sqr + U2/Kt;
%%%
Omega1 = sqrt(Omega1_sqr);
Omega2 = sqrt(Omega2_sqr);
Omega3 = sqrt(Omega3_sqr);
Omega4 = sqrt(Omega4_sqr);