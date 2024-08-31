using Pandas: read_pickle
using DelimitedFiles


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

Δt = 0.5*dx; # [dimentionless], dt*T_a = 1.7 hours #
T = DAYS*24/Tₐ; # number of hours made dimensionless 
Nt = floor(Int, T/Δt);
