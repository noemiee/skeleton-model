using Pandas: read_pickle
using DelimitedFiles
################ MODEL PARAMETERS (CONSTANTS) ######################

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


############# LOAD FORCING DATA ############################

# longitudes & dates
longitudes = read_pickle("./data/longitudes.pckl");
dates = read_pickle("./data/dates.pckl");

# extract starting and ending date for forcing profiles
d1_y, d1_m = dates[1].year, dates[1].month;
d2_y, d2_m = dates[end].year, dates[end].month;

# profiles (monthly values)
Sq_ = read_pickle("./data/Sq.pckl"); 
HbarA = read_pickle("./data/HbarA0.pckl"); 
Aₛ_ = HbarA./Hbar;
#Stheta calculation (not yet computed):  -(1-Qbar)/Qbar*Hbar.*Aₛ_ .+ 1/Qbar.*Sq_; 


################ SIMULATION PARAMETERS ######################
MULTIPLE = 3; # forcing profiles will be smoothed using a 'MULTIPLE'-months running mean 
center = Int(floor(MULTIPLE/2.0)); 

REPEAT = 5; # thermalization: repeat the simulation x times to make sure a statistical equilibrium is reached (thermalization data is later deleted)


dates_tmp = repeat(dates, REPEAT); # create repeated array of dates for reference
dates_tmp = dates_tmp[1+center:end-center]; #


MONTHS = length(dates_tmp)-1; # number of months to simulate 
DAYS = 30*(MONTHS); # number of days to simulate (we only have monthly forcing data, so we can choose for simplicity that each month has 30 days)

# initialize the domain
L = 26.6666666666666; # [dimensionless], L*D_a = 40000 km
Nx = 64; # number of grid points
dx = L/Nx; 
Δt = 0.5*dx; # [dimentionless], dt*T_a = 1.7 hours 
T = DAYS*24/Tₐ; # number of hours [dimensionless] 
Nt = floor(Int, T/Δt);
dostocha = 1; # stochastic switch


# spatial frequencies (for fourier transforms)
κ=2*π*im*fftshift(fftfreq(Nx,Nx))/L;


################ ENSO DATA PROCESSING ######################

function ONI(ENSO_val)
    """Compute ONI values by averaging over 3 months."""
    ONI_val = [mean(ENSO_val[i:i+2]) for i in 1:length(ENSO_val)-2]
    return ONI_val
end

function def_ENSO_phases(ENSO_val)
    """Assign ENSO phases based on ONI values."""
    ENSO_phases = [val >= 0.5 ? "EN" : val <= -0.5 ? "LN" : "N" for val in ENSO_val]
    return ENSO_phases
end

function replace_sequences(arr)
    """Input: array whose elements indicate the phase of ENSO (monthly). 
    Output: replaces sequences shorter than 5 'EN'or 5 'LN' with 'N'."""

    sequences_to_keep = []

    # find all sequences of at least 5 LN 
    current_sequence = []

    for (i,element) in enumerate(arr)
        if element == "LN"
            push!(current_sequence, i)
        else
            if length(current_sequence) >= 5
                append!(sequences_to_keep, current_sequence)
            end
            current_sequence = []
        end
    end
    if length(current_sequence) >= 5
        append!(sequences_to_keep, current_sequence)
    end
    
    # find all sequences of at least 5 EN
    current_sequence = []
    for (i,element) in enumerate(arr)
        if element == "EN"
            push!(current_sequence, i)
        else
            if length(current_sequence) >= 5
                append!(sequences_to_keep, current_sequence)
            end
            current_sequence = []
        end
    end
    if length(current_sequence) >= 5
        append!(sequences_to_keep, current_sequence)
    end
    
    for i in 1:length(arr)
        if !(i in sequences_to_keep)
            arr[i]="N";
        end
    end

    return arr
end


#load Nino3.4 index, compute ONI and identify El Niño and La Niña events
ENSO = reshape(readdlm("./data/Nino34.txt"),:,5);
ENSO_ym = ENSO[2:end,1:2]; # year and month
ENSO_val = ENSO[2:end,end]; # anomaly value


# compute ONI values
ONI_val = ONI(ENSO_val);
ENSO_val = ONI_val;
ENSO_ym = ENSO_ym[2:end-1,:];

# create a reference array assagning ENSO phase based on ONI
ENSO_phases = def_ENSO_phases(ENSO_val);

# count only 5 consecutive EN / LN as actual El Niño / La Niña events
ENSO_phases = replace_sequences(ENSO_phases);


################ FORCING PROFILES INTERPOLATED TO THE TIME STEP OF SIMULATION ######################

Sq     = zeros((Nt,Nx)) # at each time step, store a profile of Sq 
Stheta = zeros((Nt,Nx))
Aₛ     = zeros((Nt,Nx))

        
Aₛ_tmp = repeat(Aₛ_, REPEAT)';
Sq_tmp = repeat(Sq_, REPEAT)';
    
# smoothing of the forcing profiles via a MULTIPLE-months running mean 
mean_M_Aₛ_tmp = [mean(Aₛ_tmp[:,i:i+MULTIPLE-1], dims = 2)[:] for i in 1:size(Aₛ_tmp)[2]-MULTIPLE+1];
mean_M_Sq_tmp = [mean(Sq_tmp[:,i:i+MULTIPLE-1], dims = 2)[:] for i in 1:size(Sq_tmp)[2]-MULTIPLE+1];

Aₛ_tmp = mapreduce(permutedims, vcat, mean_M_Aₛ_tmp)';
Sq_tmp = mapreduce(permutedims, vcat, mean_M_Sq_tmp)';
   
    
#####  interpolate from monthly values to the time step of the simulation (1.7 hourly values)
# coarse grid (correspond to Sq_tmp)
space_ = 1:64;
months_coarse = 0:MONTHS; #index
    
# fine grid
hours_ = range(Δt, step = Δt, length = Nt);
months_fine = hours_*Tₐ/24.0/30; # nth simulated hours correpond to mth simulated months

# Interpolate Sq
itp = LinearInterpolation((space_, months_coarse), Sq_tmp);
Sq_interp = [itp(y,x) for y in space_, x in months_fine];
Sq[:,:]=Sq_interp';
    
# Interpolate As
itp = LinearInterpolation((space_, months_coarse), Aₛ_tmp);
Aₛ_interp = [itp(y,x) for y in space_, x in months_fine];
Aₛ[:,:]=Aₛ_interp';
    
# shift mean of Sq
for t in 1:Nt
    Sq[t,:] = Sq[t,:] .- mean(Sq[t,:]) .+ mean(Aₛ[t,:])*Hbar;
end
    
# compute Stheta
Stheta = -(1-Qbar)/Qbar*Hbar.*Aₛ .+ 1/Qbar.*Sq;


# compute temporal mean of profiles
mean_Aₛ = mean(Aₛ, dims = 1);
mean_Sq = mean(Sq, dims = 1);
mean_Stheta = mean(Stheta, dims = 1);