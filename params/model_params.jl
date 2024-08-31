################ MODEL PARAMETERS ######################

Hbar = 0.22; #ok (HH)
Qbar = 0.9; #ok (QQ)
Γ =  0.6; #ok (GG)
γ = 1.0/π^(1.0/4.0); # comes from projection on first Hermite modes meridionally
Abar = 0.1331; # comes form projection of Abar on first Hermite mode, Stheta/Hbar/γ
Δa = 0.001;
Sref = 0.0220; # external forcing

# dimensionalization constants
Dₐ = 1500; #[km]
Tₐ = 8; #[h]
Θₐ = 15; #[Kelvin]
Uₐ = 50; #[m/s]  

# time varying profiles
SET_TIME_VARIATION = true

# stochastic switch
dostocha = 1;