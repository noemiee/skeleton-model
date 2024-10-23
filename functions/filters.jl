module filters

export filter
export lanczosfilter

#using Statistics
#using FFTW


function lanczosfilter(x,dT,Cf,M,pass)
    
#  LANCZOSFILTER   Filters a time series via Lanczos method (cosine filter).  
#  [Y,coef,window,Cx,Ff] = LANCZOSFILTER(X,dT,Cf,M,pass) Filters the time
#  series via the Lanczos filter in the frequency space (FFT), where
#
#   INPUTS:
#      X    - Time series
#      dT   - Sampling interval       (default: 1)
#      Cf   - Cut-off frequency       (default: half Nyquist)
#      M    - Number of coefficients  (default: 100)
#      pass - Low or high-pass filter (default: 'low')
#
#   OUTPUTS:
#      Y      - Filtered time series
#      coef   - Coefficients of the time window (cosine)
#      window - Frequency window (aprox. ones for Ff lower(greater) than Fc 
#               if low(high)-pass filter and zeros otherwise)
#      Cx     - Complex Fourier Transform of X for Ff frequencies
#      Ff     - Fourier frequencies, from 0 to the Nyquist frequency.
#  
#   The program removes from the time series the frequencies greater than   
#   the cut off frequency if "pass" is 'low', i.e., low-pass filter .
#   Otherwise, if pass is 'high', frequencies from zero to Cf are removed,
#   i.e., a high-pass filter. Units of the cut-off frequency, [Cf], must be
#   [dT]^{-1}. 
#  
#   In consequence, when used as a low-pass the time series is smoothed   
#   like a cosine filter in time space with M coefficients where greater is
#   better.    
    
    Nd =  size(x)[1];
    
    Nf = 1/(2*dT); # Nyquist frequency

    # Normalize the cut off frequency with the Nyquist frequency:
    Cf = Cf/Nf;

    # Lanczos coeficients:
    hkcs = Cf*[1; sin.(pi*(1:M)*Cf)./(pi*(1:M)*Cf)];
    sigma = [1; sin.(pi*(1:M)/M)./(pi*(1:M)/M)];
    hkB = hkcs.*sigma; # defines the low-pass filter
    
    hkA = -hkB; hkA[1] = hkA[1]+1; # defines the high-pass filter   
    
    # select the desired filter
    if pass == "high"
        coef = hkA;
    else
        coef = hkB;
    end;
    
    # create Window of cosine filter in frequency space.
    Ff = fftshift(fftfreq(Nd,Nd)); 
    window = zeros(length(Ff));
    
    for i = 1:length(Ff)
        window[i] = coef[1] + 2*sum(coef[2:end] .* cos.(2*pi*(1:length(coef)-1)*Ff[i]/Nd)); # implementation of (7)
    end
    
    # Filtering in frequency space is multiplication, (convolution in time space).
    Cx  = fftshift(fft(x)); 
    CxH = zeros(Nd)*1im;
    CxH = Cx.*window;
    y = real(ifft(ifftshift(CxH)));
    
    return y,coef,window,Cx,Ff
    
end


function filter(month_array, u, theta, q, ha, Qbar, Hbar, γ, Tₐ, times, Nt, Nx, NUM_REPEAT)
    
    """ Performs the 20-100 days filtering method described in Stachnik et al. 2015

    INPUTS: 
    - month_array: array of months corresponding to each time step of simulation (based on forcing)
    - u: zonal wind
    - theta: potential temperature
    - q: specific humidity
    - ha: convective activity
    - Qbar: mean background vertical moisture gradient
    - Hbar: scaling constant for convective activity
    - γ: parameter, comes from projection on first Hermite modes meridionally
    - Tₐ: dimensionalization constant for time
    - times: simulation time [non-dimensional] corresponding to the data provided, separated by Δt 
    - Nt: number of time steps
    - Nx: number of grid points
    - NUM_REPEAT: number of independent simulation runs
    
    """
    
    #step 1: daily mean data 
    sim_day = floor.(times*Tₐ/24) # simulation day corresponding to each time step 
    u_daily = zeros(Nx, NUM_REPEAT, length(unique(sim_day)))
    ha_daily = zeros(Nx, NUM_REPEAT, length(unique(sim_day)))
    q_daily = zeros(Nx, NUM_REPEAT, length(unique(sim_day)))
    theta_daily = zeros(Nx, NUM_REPEAT, length(unique(sim_day)))
    months_tracked = zeros(length(unique(sim_day))) # for each simulation day, the corresponding month
    
    i=1
    for d in unique(sim_day)
        idxs = findall(sim_day .== d);
        months_tracked[i] = mean(month_array[idxs]); 
        i+=1;
    end
    
    for nr in 1:NUM_REPEAT
        i=1
        for d in unique(sim_day)
            idxs = findall(sim_day .== d)
            u_daily[:,nr,i] = mean(u[:,nr,idxs], dims = 2);
            ha_daily[:,nr,i] = mean(ha[:,nr,idxs], dims = 2);
            q_daily[:,nr,i] = mean(q[:,nr,idxs], dims = 2);
            theta_daily[:,nr,i] = mean(theta[:,nr,idxs], dims = 2);
            i+=1
        end
    end
    
    sim_day = unique(sim_day);
    Nd = length(sim_day); 
    println("STEP: daily means computed")
    
    # step 2: remove the long term mean OK
    #u_daily = u_daily .- mean(u_daily, dims = 3);
    #ha_daily = ha_daily .- mean(ha_daily, dims = 3);
    
    u_daily_anom = zeros(size(u_daily));
    ha_daily_anom = zeros(size(ha_daily));
    q_daily_anom = zeros(size(q_daily));
    theta_daily_anom = zeros(size(theta_daily));
    
    
    # step 3: remove mean and first three harmonics of the annual cycle (ie periods 12 months, 6 months and 4 months, corresponding to 1,2 and 3 cycle per year)
    for nr in 1:NUM_REPEAT
        for pt in 1:64 # filter each spatial point separately
            uH = fftshift(fft(u_daily[pt,nr,:]));
            haH = fftshift(fft(ha_daily[pt,nr,:]));
            qH = fftshift(fft(q_daily[pt,nr,:]));
            thetaH = fftshift(fft(theta_daily[pt,nr,:]));
            
            freqs = fftshift(fftfreq(Nd, Nd))/Nd;
            idx1 = findall(x-> abs(x) <= 3/365.0, freqs);
            uH[idx1].=0;
            haH[idx1].=0;
            qH[idx1].=0;
            thetaH[idx1].=0;
        
            u_daily_anom[pt,nr,:] = real(ifft(ifftshift(uH)));
            ha_daily_anom[pt,nr,:] = real(ifft(ifftshift(haH)));
            q_daily_anom[pt,nr,:] = real(ifft(ifftshift(qH)));
            theta_daily_anom[pt,nr,:] = real(ifft(ifftshift(thetaH)));
        end
    end
    println("STEP: removed mean and first three harmonics of the annual cycle")
    
    # step 4: apply Lanczos cosine filder with cutoff 20–100 days and 201 weights 
    M = 100 # that is 2M+1=201 weights
    fl = 1/20.0;
    fh = 1/100.0;
    
    yu = zeros(size(u_daily_anom));
    yha = zeros(size(ha_daily_anom ));
    yq = zeros(size(q_daily_anom ));
    ytheta = zeros(size(theta_daily_anom ));
    
    for nr in 1:NUM_REPEAT
        for pt in 1:64 # filter each spatial point separately
            # first apply low pass filter
            yu_tmp, coef, window, Cx, Ff = lanczosfilter(u_daily_anom[pt,nr, :], 1, fl, M, "low")
            # then apply high pass filter
            yu[pt,nr,:], coef, window, Cx, Ff = lanczosfilter(yu_tmp,1,fh,M,"high")
            
            # low pass filter for u
            yha_tmp, coef, window, Cx, Ff = lanczosfilter(ha_daily_anom[pt,nr, :], 1, fl, M, "low")
            # high pass filter for ha
            yha[pt,nr,:], coef, window, Cx, Ff = lanczosfilter(yha_tmp,1,fh,M,"high")
            
            # first apply low pass filter
            yq_tmp, coef, window, Cx, Ff = lanczosfilter(q_daily_anom[pt,nr, :], 1, fl, M, "low")
            # then apply high pass filter
            yq[pt,nr,:], coef, window, Cx, Ff = lanczosfilter(yq_tmp,1,fh,M,"high")
            
            # first apply low pass filter
            ytheta_tmp, coef, window, Cx, Ff = lanczosfilter(theta_daily_anom[pt,nr, :], 1, fl, M, "low")
            # then apply high pass filter
            ytheta[pt,nr,:], coef, window, Cx, Ff = lanczosfilter(ytheta_tmp,1,fh,M,"high")
        end
    end
    
    u_daily_filtered = yu;
    ha_daily_filtered = yha;
    q_daily_filtered = yq;
    theta_daily_filtered = ytheta;
    println("STEP: applied Lanczos cosine filder with cutoff 20–100 days and 201 weights")
    
    return sim_day, months_tracked, u_daily, u_daily_filtered, ha_daily, ha_daily_filtered, q_daily, q_daily_filtered, theta_daily, theta_daily_filtered
  
end

end