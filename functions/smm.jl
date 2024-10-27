module smm

export SMM

using LinearAlgebra
using Statistics



function SMM(sim_days, u_daily_filtered, ha_daily_filtered)
    """
       Implementation of the Skeleton Multivariate MJO index, as presented in [Stachnik et al., 2015](https://doi.org/10.1002/2015JD023916) based on the Real-Time Multivariate MJO index (RMM) by [Wheeler and Hendon (2004)](https://doi.org/10.1175/1520-0493(2004)132%3C1917:AARMMI%3E2.0.CO;2).
    """
    # INPUT: daily data, filtered
    Nd = size(u_daily_filtered)[end];
    
    # step 1: normalize each anomaly field by its (global) standard deviation
      
    tmp_u = zeros(size(u_daily_filtered)[1]);
    tmp_ha = zeros(size(ha_daily_filtered)[1]);
    
    for pt in 1:64
        tmp_u[pt] = var(u_daily_filtered[pt,:]);
        tmp_ha[pt] = var(ha_daily_filtered[pt,:]);
    end
    
    
    u_daily_filtered_normalized = u_daily_filtered ./ sqrt(mean(tmp_u));
    ha_daily_filtered_normalized = ha_daily_filtered ./ sqrt(mean(tmp_ha));
    
    # step 2: compute multivariate PCA of data
    state = vcat(u_daily_filtered_normalized, -ha_daily_filtered_normalized) # ~ wind , OLR


    # 2.1 compute corelation matrix
    CCM = state*transpose(state);#/Nd;
    
    # 2.2 compute SVD
    F = svd(CCM);
    V, EV, U = F; 
      
    # 2.3 first two principal components
    PC1 = -U[:,1];
    PC2 = U[:,2];
    
    # 2.4 variance explained
    var1 = EV[1]/sum(EV);
    var2 = EV[2]/sum(EV);
    
    # the first two PC are:
    # PC1u = PC1[1:64];
    # PC1ha = PC1[65:end];
    # PC2u = PC2[1:64];
    # PC2ha = PC2[65:end];

        
    # step 3: projection onto the principal components
    SMM1 = zeros(length(sim_days));
    SMM2 = zeros(length(sim_days));
    for i = 1:length(sim_days)
        SMM1[i] = dot(state[:,i], PC1);
        SMM2[i] = dot(state[:,i], PC2);
    end
        
    SMM1 = (SMM1 .- mean(SMM1))/std(SMM1);
    SMM2 = (SMM2 .- mean(SMM2))/std(SMM2);
    
    amplitude = sqrt.(SMM1.^2 .+ SMM2.^2);
    phase = angle.( SMM1 .+ SMM2.*1im);
    
    return var1, var2, PC1, PC2, amplitude, phase, SMM1, SMM2
   
end



end


