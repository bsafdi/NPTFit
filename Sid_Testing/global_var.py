# Define units, with GeV as base unit
GeV = 10**6;
eV = 10**-9*GeV;
KeV = 10**-6*GeV;
MeV = 10**-3*GeV;
TeV = 10**3*GeV;

Sec = (1/(6.582119*10**-16))/eV; 
Kmps = 3.3356*10**-6;
Centimeter = 5.0677*10**13/GeV;
Meter = 100*Centimeter;
Km = 10**5*Centimeter;
Kilogram = 5.6085*10**35*eV;
Day = 86400*Sec;
Year = 365*Day;
KgDay = Kilogram*Day;
amu = 1.66053892*10**-27*Kilogram;
Mpc = 3.086*10**24*Centimeter;

# Particle and astrophysics parameters
sigma_v=3*10**-26*(Centimeter**3*Sec**-1)
m_chi=100*GeV
omega_m=0.23
omega_lambda=0.723
rho_c=9.9*10**-27*(Kilogram*Meter**-3);
rho_m=omega_m*rho_c
delta_c=1.69
h=0.7
M_s=1.99*10**30*(Kilogram)
M_star=2*10**12*(h**-1*M_s)

# 2MASS redshift parameters
m=1.90
beta=1.75
z0=0.07
dA = 288.4*(Mpc)

# Halo occupation distribution parameters
Mmin = 10**11.68*M_s
sigmaMmin = .15