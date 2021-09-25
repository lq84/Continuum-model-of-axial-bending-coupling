syms t L n Lelm h E
Iz=h*t^3/12;                                                               % Moment of inertia (m^4)
A=t*h;                                                                     % Cross section area (m^2)
a0=Lelm/2;
K_elmEuler=[A*E/2/a0,             0,             0, -A*E/2/a0,             0,             0;
          0,  3*E*Iz/2/a0^3,  3*E*Iz/2/a0^2,        0, -3*E*Iz/2/a0^3,  3*E*Iz/2/a0^2;
          0,  3*E*Iz/2/a0^2,      2*E*Iz/a0,        0, -3*E*Iz/2/a0^2,        E*Iz/a0;
   -A*E/2/a0,             0,             0,  A*E/2/a0,             0,             0;
          0, -3*E*Iz/2/a0^3, -3*E*Iz/2/a0^2,        0,  3*E*Iz/2/a0^3, -3*E*Iz/2/a0^2;
          0,  3*E*Iz/2/a0^2,        E*Iz/a0,        0, -3*E*Iz/2/a0^2,     2*E*Iz/a0]           % Euler beam

s=4;c=1/2;
s1=s*(1+c);s2=2*s1;
k_elm1=E*t/Lelm*[1,0,0,-1,0,0;
    0,0,0,0,0,0;
    0,0,0,0,0,0;
    -1,0,0,1,0,0;
    0,0,0,0,0,0;
    0,0,0,0,0,0];
k_elm2=E*t^3/Lelm/12*[0,0,0,0,0,0;
    0,s2/Lelm^2,s1/Lelm,0,-s2/Lelm^2,s1/Lelm;
    0,s1/Lelm,s,0,-s1/Lelm,s*c;
    0,0,0,0,0,0;
    0,-s2/Lelm^2,-s1/Lelm,0,s2/Lelm^2,-s1/Lelm;
    0,s1/Lelm,s*c,0,-s1/Lelm,s];
K_elm=k_elm1+k_elm2