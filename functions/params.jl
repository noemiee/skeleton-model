using Pandas: read_pickle


function smooth_background(field, num_k, alpha, Nx)
    idx0 = Int(floor(Nx/2+1));
    
    rft = fftshift(fft(field));
    rft[1:idx0-num_k-1] .= 0.0;
    rft[idx0+num_k+1:end] .= 0.0;
    rft[idx0+1:end] .= alpha*rft[idx0+1:end];
    rft[1:idx0-1] .= alpha*rft[1:idx0-1];
    Field_smoothed = ifft(ifftshift(rft));
       
    return real.(Field_smoothed)

end




#module params

Hbar = 0.22; #ok (HH)
Qbar = 0.9; #ok (QQ)
Γ =  1.66; #ok (GG)
γ = 1.0/π^(1.0/4.0); # comes from projection on first Hermite modes meridionally
Abar = 0.1331; # comes form projection of Abar on first Hermite mode, Stheta/Hbar/γ
Δa = 0.001;
Sref = 0.0220; # external forcing

# dimensionalization constants
Dₐ = 1500; #[km]
Tₐ = 8; #[h]
Θₐ = 15; #[Kelvin]
Uₐ = 50; #[m/s]  

## SIMULATION specific parameters
YEARS = 34+7;
DAYS = 365*YEARS;

# initialize the domain
L = 26.6666666666666; # [dimensionless], L*D_a = 40000 km
Nx = 64; # number of grid points
dx = L/Nx; 
Δt = 0.5*dx; # [dimentionless], dt*T_a = 1.7 hours #
tlength =  DAYS; # total simulation time in days
T = tlength*24/Tₐ; # [dimensionless] 
Nt = floor(Int, T/Δt);

# domain frequencies
κ=2*π*im*fftshift(fftfreq(Nx,Nx))/L;

# stochastic switch
dostocha = 1;
# repetitions (to reach equilibrium)
num_rep = 1;
# warm pool or homogeneous
background_type = "EN"; # warmpool "WP"; homogeneous "H"; La Nina "LN"; El Nino "EN"; mean "M"

# smoothing parameters for the background
alpha = 1.0;
num_k = Int(floor(Nx/2)); 







#longitudes = [  0.   ,   5.625,  11.25 ,  16.875,  22.5  ,  28.125,  33.75 ,  39.375,  45.   ,  50.625,  56.25 ,  61.875,  67.5  ,  73.125,  78.75 ,  84.375,  90.   ,  95.625, 101.25 , 106.875, 112.5  , 118.125, 123.75 , 129.375, 135.   , 140.625, 146.25 , 151.875, 157.5  , 163.125, 168.75 , 174.375, 180.   , 185.625, 191.25 , 196.875, 202.5  , 208.125, 213.75 , 219.375, 225.   , 230.625, 236.25 , 241.875, 247.5  , 253.125, 258.75 , 264.375, 270.   , 275.625, 281.25 , 286.875, 292.5  , 298.125, 303.75 , 309.375, 315.   , 320.625, 326.25 , 331.875, 337.5  , 343.125, 348.75 , 354.375];
#
#ref_Stheta = [0.04154562, 0.0320926 , 0.0249229 , 0.01752242, 0.01760054, 0.02066146, 0.04039979, 0.05863459, 0.06511504, 0.05947743,  0.05207649, 0.04433111, 0.03706885, 0.03393769, 0.02595959,  0.02234349, 0.01979342, 0.01028133, 0.00931022, 0.00943257,  0.02130727, 0.0218136 , 0.0185833 , 0.01967533, 0.02087887,  0.02535353, 0.025211  , 0.0258603 , 0.02600104, 0.0291941 ,  0.03246797, 0.03656712, 0.04075505, 0.04426148, 0.04824795,  0.05017088, 0.0501721 , 0.05060286, 0.05144824, 0.05245601,  0.05289322, 0.05300512, 0.05340075, 0.05283576, 0.05083652,  0.04873865, 0.0457079 , 0.04504013, 0.0397699 , 0.03002696,  0.01857897, 0.01296339, 0.01646208, 0.02146161, 0.03175237,  0.03743654, 0.04024113, 0.04577519, 0.05150002, 0.04900662,  0.04595611, 0.04492111, 0.0458396 , 0.04497065];
#ref_Sq = [0.02726517, 0.02514457, 0.02165382, 0.0227487 , 0.02275289, 0.0222607 , 0.02004985, 0.02373933, 0.0265001 , 0.03613674, 0.03956223, 0.04211761, 0.04263723, 0.03931252, 0.03945283, 0.04429322, 0.04358685, 0.03994474, 0.03790025, 0.03722263, 0.03769172, 0.03289295, 0.03286652, 0.03674719, 0.03250908, 0.03596094, 0.0367345 , 0.04147149, 0.04095047, 0.04034827, 0.04179791, 0.04296932, 0.0429155 , 0.04354185, 0.04312555, 0.04414219, 0.04284789, 0.04146934, 0.04032892, 0.03919208, 0.03919565, 0.03862465, 0.03772888, 0.03794295, 0.03705076, 0.03598919, 0.03563955, 0.03528401, 0.03480808, 0.03576178, 0.02770517, 0.03140604, 0.03275321, 0.0344637 , 0.03625075, 0.03718675, 0.03706563, 0.03714983, 0.04365295, 0.04222641, 0.03950713, 0.03745733, 0.03063784, 0.02946446];
#ref_HbarA = [0.01406036, 0.0185273 , 0.02320315, 0.0288911 , 0.03160744, 0.03019016, 0.02381693, 0.01916279, 0.01462113, 0.01702146, 0.02244641, 0.02996128, 0.03432718, 0.04281208, 0.04305818, 0.047717  , 0.05589337, 0.05505773, 0.05252596, 0.04411708, 0.0550125 , 0.04585025, 0.04267429, 0.0474713 , 0.05534015, 0.06238737, 0.05972578, 0.05947293, 0.05665611, 0.05429918, 0.05094055, 0.04839706, 0.04604246, 0.04322603, 0.04058331, 0.03736818, 0.03333633, 0.03076613, 0.02871798, 0.02734584, 0.02577249, 0.02475185, 0.02459089, 0.02450548, 0.02273799, 0.02396189, 0.02380518, 0.02491772, 0.02485025, 0.02742422, 0.04534358, 0.05110595, 0.05438311, 0.0489477 , 0.04557988, 0.04597347, 0.03568772, 0.03282235, 0.03408303, 0.02701412, 0.02080564, 0.01829342, 0.01811061, 0.01337537];
#
#ElNino_Sq = [0.03485808, 0.03357372, 0.02782994, 0.03027396, 0.03035573, 0.02888842, 0.02587732, 0.02929485, 0.03433139, 0.0453472 , 0.04796958, 0.04885229, 0.04732569, 0.04337693, 0.04449205, 0.05009883, 0.04897843, 0.04101416, 0.03888658, 0.03959165, 0.03768122, 0.03240462, 0.03518403, 0.03956837, 0.03559972, 0.03927183, 0.04390573, 0.04917553, 0.04779815, 0.04653629, 0.04695374, 0.04881964, 0.04463026, 0.04343506, 0.04360234, 0.04632216, 0.04561594, 0.04613553, 0.04650837, 0.04350419, 0.04455723, 0.04324498, 0.04353053, 0.04460993, 0.04291537, 0.04230839, 0.04061377, 0.03980201, 0.0391851 , 0.04153647, 0.02646615, 0.03543324, 0.03550563, 0.03566788, 0.03878146, 0.0404982 , 0.04010412, 0.04030564, 0.04439997, 0.04225567, 0.03979417, 0.03768464, 0.03491309, 0.03459447];
#ElNino_Stheta = [0.05302188, 0.0399836 , 0.02722046, 0.0134926 , 0.00577522, 0.01013003, 0.03196816, 0.04936257, 0.05888103, 0.05440039, 0.04369241, 0.03568007, 0.03291927, 0.03623655, 0.03541638, 0.04133143, 0.04006312, 0.02818544, 0.02946256, 0.02959027, 0.04424299, 0.04565974, 0.05052347, 0.0475151 , 0.04582552, 0.04725641, 0.04858047, 0.04558472, 0.04215833, 0.03968005, 0.04082813, 0.03976639, 0.0377135 , 0.03797732, 0.04181236, 0.04259597, 0.03784945, 0.04108122, 0.04569841, 0.04865417, 0.05249142, 0.05859304, 0.05867031, 0.05747067, 0.05269442, 0.04881345, 0.04729043, 0.04555841, 0.04105076, 0.03223439, 0.024043  , 0.01563159, 0.02298706, 0.02978971, 0.04571015, 0.05283925, 0.0520417 , 0.05816038, 0.06388509, 0.0618505 , 0.05885364, 0.05662971, 0.05604143, 0.05617347];
#ElNino_HbarA = [0.01568804, 0.02297661, 0.03029768, 0.03045468, 0.03399403, 0.03593604, 0.03140456, 0.03109974, 0.02425493, 0.03094491, 0.03485219, 0.03571092, 0.03868731, 0.04277387, 0.04193643, 0.04272534, 0.04191456, 0.0405206 , 0.03965415, 0.03656723, 0.04815627, 0.02735893, 0.02389219, 0.0273092 , 0.03994915, 0.04906644, 0.04876072, 0.04943306, 0.05379994, 0.05802396, 0.06593647, 0.06640395, 0.0649728 , 0.0684695 , 0.07043636, 0.07188166, 0.06349066, 0.06451233, 0.06024417, 0.06024415, 0.05851117, 0.05714404, 0.05295146, 0.05282901, 0.04812159, 0.04483278, 0.04331806, 0.04136075, 0.0418665 , 0.04289979, 0.05929308, 0.05432805, 0.05564203, 0.04264765, 0.03705004, 0.0396402 , 0.03038207, 0.02736632, 0.02783105, 0.02425455, 0.01917257, 0.01859165, 0.01941458, 0.01302864];
#
#LaNina_Sq = [0.03309182, 0.03106602, 0.02785924, 0.02898068, 0.02896235, 0.02778134, 0.02138993, 0.02479244, 0.03183026, 0.04457253, 0.04662371, 0.04945878, 0.04974927, 0.04651384, 0.04708837, 0.05154306, 0.04844466, 0.04369585, 0.04072143, 0.03948637, 0.03881758, 0.0348689 , 0.03492033, 0.04135198, 0.03608263, 0.04129681, 0.04001316, 0.0421138 , 0.04061919, 0.04032544, 0.04472668, 0.04691949, 0.04829904, 0.04902804, 0.04815716, 0.04746231, 0.04465795, 0.04354803, 0.04219303, 0.04022698, 0.04014235, 0.03902353, 0.03749143, 0.03805176, 0.03714911, 0.0363861 , 0.03718835, 0.03667386, 0.03556389, 0.03797487, 0.02716861, 0.03528281, 0.03643363, 0.03800974, 0.03987365, 0.04160735, 0.04199453, 0.04283852, 0.0499253 , 0.04707233, 0.04303816, 0.03999965, 0.03568165, 0.03484214];
#LaNina_Stheta = [ 0.0420957 ,  0.02807712,  0.0150543 ,  0.00207794,  0.00098759, 0.00248897,  0.03409497,  0.06321969,  0.07203212,  0.06695276,  0.05632127,  0.04586337,  0.03494404,  0.02922719,  0.01993888,  0.01690665,  0.01035884, -0.00290248, -0.00523507, -0.00427787,  0.004453  ,  0.00373526,  0.00079644,  0.00649897,  0.01234832,  0.02340618,  0.02946447,  0.03441455,  0.03835481,  0.04365517,  0.0483617 ,  0.05351431,  0.0586179 ,  0.06202036,  0.06509613,  0.06698069,  0.06775327,  0.06825169,  0.06842936,  0.06842471,  0.06838364,  0.06728693,  0.06736906,  0.06544788,  0.06262897,  0.06000174,  0.05754634,  0.05505549,  0.04786913,  0.0357667 ,  0.01827691,  0.01037303,  0.01259553,  0.01522793,  0.03200489,  0.04209244,  0.04621858,  0.05520975,  0.06476909,  0.06127025,  0.05628623,  0.05261619,  0.05184237,  0.04850439];
#LaNina_HbarA = [0.01769629, 0.02517658, 0.03305585, 0.03776986, 0.03906812, 0.03615473, 0.02516112, 0.0178995 , 0.01031444, 0.01281991, 0.01900407, 0.02962548, 0.03614791, 0.04674147, 0.04916417, 0.05510205, 0.0656146 , 0.06956754, 0.06932959, 0.05655988, 0.07151858, 0.05987336, 0.05507706, 0.06188183, 0.07164916, 0.07622532, 0.06658868, 0.06132896, 0.05505219, 0.04909395, 0.04331977, 0.03848488, 0.0354518 , 0.03225312, 0.02973777, 0.02639998, 0.02436148, 0.02281916, 0.02328615, 0.02282624, 0.02105112, 0.02063064, 0.0222195 , 0.02245985, 0.02138967, 0.02186662, 0.02375302, 0.02408632, 0.02420218, 0.02932977, 0.05994108, 0.06342795, 0.06623351, 0.05887092, 0.05286814, 0.05286424, 0.03901672, 0.03676776, 0.04148018, 0.0350502 , 0.02766202, 0.02368158, 0.02464221, 0.01738163];

#load longitudes
longitudes = read_pickle("./functions/longitudes.pckl");

# select the profile for simulations
if background_type == "WP"
    Stheta = Sref*(1.0 .-2*0.3 .*cos.(2*pi/L.*(0:Nx-1).*dx));
    Stheta = circshift(Stheta,-4);
    Sq = Stheta;
    Aₛ = 1.0/ Hbar .*(Sq .- Qbar*Stheta)/(1.0-Qbar);
elseif background_type == "H"
    Stheta = Sref*ones(64);
    Sq = Stheta;
    Aₛ = 1.0/ Hbar .*(Sq .- Qbar*Stheta)/(1.0-Qbar);
    
else
    # load profiles
    ref_Stheta = read_pickle("./functions/Stheta.pckl");
    ref_Sq = read_pickle("./functions/Sq.pckl");
    ref_HbarA = read_pickle("./functions/HbarA0.pckl");
        
    #smoothed_ref_Stheta = smooth_background(ref_Stheta, num_k, alpha, Nx);
    #smoothed_ref_A = smooth_background(ref_HbarA./Hbar, num_k, alpha, Nx);
    #Stheta = smoothed_ref_Stheta;
    #Aₛ = smoothed_ref_A .- mean(smoothed_ref_A) .+ mean(Stheta)./Hbar; 
    #Sq = Hbar.*Aₛ .- Qbar.*(Hbar.*Aₛ .- Stheta);
    
    smoothed_ref_Sq = smooth_background(ref_Sq, num_k, alpha, Nx);
    smoothed_ref_A = smooth_background(ref_HbarA./Hbar, num_k, alpha, Nx);
    Sq = smoothed_ref_Sq;
    Aₛ = smoothed_ref_A .- mean(smoothed_ref_A) .+ mean(Sq)./Hbar; 
    Stheta = -(1-Qbar)/Qbar*Hbar.*Aₛ .+ 1/Qbar.*Sq; 
end;