using Pandas: read_pickle
using DelimitedFiles

SET_TIME_VARIATION = true

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


############# FORCING PARAMETERS ############################

# load longitudes
longitudes = read_pickle("./data/longitudes.pckl");

# load dates
dates = read_pickle("./data/dates.pckl");

# starting and ending date for the forcing profiles
d1_y = dates[1].year;
d1_m = dates[1].month;
d2_y = dates[end].year;
d2_m = dates[end].month;


MULTIPLE = 3;
center = Int(floor(MULTIPLE/2.0));
REPEAT = 5; # (the first repeated periods are for thermalization and are later deleted)

# array of dates for reference
dates_tmp = repeat(dates, REPEAT); 
dates_tmp = dates_tmp[1+center:end-center];

################ SIMULATION PARAMETERS ######################

MONTHS = length(dates_tmp)-1; # number of months to simulate 
DAYS = 30*(MONTHS); # number of days to simulate (we only have monthly forcing data, so we can choose for simplicity that each month has 30 days)

# initialize the domain
L = 26.6666666666666; # [dimensionless], L*D_a = 40000 km
Nx = 64; # number of grid points
dx = L/Nx; 
Δt = 0.5*dx; # [dimentionless], dt*T_a = 1.7 hours #
T = DAYS*24/Tₐ; # number of hours made dimensionless 
Nt = floor(Int, T/Δt);
# stochastic switch
dostocha = 1;


# frequencies
κ=2*π*im*fftshift(fftfreq(Nx,Nx))/L;


################ SOME FUNCTIONS   ######################



function replace_sequences(arr)
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


function get_date_index(data_ym, d1_y, d1_m, d2_y, d2_m)
    
    idx1 = -100000;
    for i in eachindex(data_ym[:,1])
        if data_ym[i, 1] == d1_y && data_ym[i, 2] == d1_m
            idx1 = i;
        end
    end
    
    idx2 = -100000;
    for i in eachindex(data_ym[:,1])
        if data_ym[i, 1] == d2_y && data_ym[i, 2] == d2_m
            idx2 = i;
        end
    end
    
    return idx1, idx2
end

function def_ENSO_phases(ENSO_val)
    # El Niño -> anomaly above 0.5 / La Niña -> anomaly below 0.5 
    ENSO_phases = Array{String}(undef, size(ENSO_val));
    for (i, val) in enumerate(ENSO_val)
        if (val >=0.5)
            ENSO_phases[i]="EN"
        elseif (val <= -0.5) 
            ENSO_phases[i]="LN"
        else
            ENSO_phases[i] = "N"
        end
    end
    return ENSO_phases
end

function ONI(ENSO_val)
    ONI_val = zeros(length(ENSO_val)-2);
    for i in 1:length(ENSO_val)-2
        ONI_val[i] = mean(ENSO_val[i:i+2])
    end
    return ONI_val
end


################ FORCING PROFILES AND ENSO DATA ######################


################################# load Nino3.4 index

ENSO = reshape(readdlm("./data/Nino34.txt"),:,5);
ENSO_ym = ENSO[2:end,1:2]; # year and month
ENSO_val = ENSO[2:end,end]; # anomaly value

# select the time slice that we are interested in
#idx1, idx2 = get_date_index(ENSO_ym, d1_y, d1_m, d2_y, d2_m);

#ENSO_ym = ENSO_ym[idx1-1:idx2+1,:];
#ENSO_val = ENSO_val[idx1-1:idx2+1];


# compute ONI values
ONI_val = ONI(ENSO_val);
ENSO_val = ONI_val;
ENSO_ym = ENSO_ym[2:end-1,:];

# create a reference array indicating if the system was in El Niño or La Niña or Neutral state
ENSO_phases = def_ENSO_phases(ENSO_val);

# count only 4 consecutive EN / LN as actual El Niño / La Niña events
ENSO_phases = replace_sequences(ENSO_phases);

########################################

# load profiles for simulations
Sq_ = read_pickle("./data/Sq.pckl"); # monthly values


HbarA = read_pickle("./data/HbarA0.pckl"); # monthly values
Aₛ_ = HbarA./Hbar;

#Stheta_ = -(1-Qbar)/Qbar*Hbar.*Aₛ_ .+ 1/Qbar.*Sq_; 

Sq     = zeros((Nt,Nx)) # store a profile of Sq at each time step
Stheta = zeros((Nt,Nx))
Aₛ     = zeros((Nt,Nx))

ENSO_phases_monthly = Array{String}(undef, MONTHS+1);

################################

if SET_TIME_VARIATION == true
    #Sq_tmp     = zeros((Nx, MONTHS+1)); # monthly values of Sq (+1 necessary for interpolation, see below)
    #Aₛ_tmp     = zeros((Nx, MONTHS+1));
        
    Aₛ_tmp = repeat(Aₛ_, REPEAT)';
    Sq_tmp = repeat(Sq_, REPEAT)';
    
    # MULTIPLE-months running mean for slower change in forcing profiles (essentially smoothing forcing profiles)
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
    
    
    #########################################
    
    # select the corresponding ENSO reference 
    #ENSO_phases = ENSO_phases[1+Int(floor(MULTIPLE/2)):size(ENSO_phases)[1]-MULTIPLE+1+Int(floor(MULTIPLE/2))]; 
    
    #ENSO_phases_monthly[1:num_vals*REPEAT] = repeat(ENSO_phases, REPEAT);
    
    
else
    REPEAT = Nt;
    Sq[:,:]     = repeat(mean_Sq, REPEAT);
    Stheta[:,:] = repeat(mean_Stheta, REPEAT);
    Aₛ[:,:]     = repeat(mean_Aₛ, REPEAT);
end

# compute temporal mean of profiles
mean_Aₛ = mean(Aₛ, dims = 1)
mean_Sq = mean(Sq, dims = 1)
mean_Stheta = mean(Stheta, dims = 1)