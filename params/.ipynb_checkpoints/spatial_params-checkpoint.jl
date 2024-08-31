#  domain
L = 26.6666666666666; # [dimensionless], L*D_a = 40000 km 
Nx = 64; # number of grid points
dx = L/Nx; 

# frequencies
κ=2*π*im*fftshift(fftfreq(Nx,Nx))/L;