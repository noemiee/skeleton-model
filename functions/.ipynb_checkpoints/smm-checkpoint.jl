module smm

export SMM

using LinearAlgebra
using Statistics
using JLD
using Interpolations

using Pandas: read_pickle


function SMM(sim_days, u_daily_filtered, ha_daily_filtered, background_type)
    # INPUT: daily data, filtered
    Nd = size(u_daily_filtered)[end];
    
    # step 1: normalize each anomaly field by its (global) standard deviation
    #tmp_u = zeros(size(u_daily_filtered));
    #tmp_ha = zeros(size(ha_daily_filtered));
    
    tmp_u = zeros(size(u_daily_filtered)[1]);
    tmp_ha = zeros(size(ha_daily_filtered)[1]);
    
    for pt in 1:64
        #tmp_u[pt,:] = (u_daily_filtered[pt,:].-mean(u_daily_filtered))./std(u_daily_filtered)
        #tmp_ha[pt,:] = (ha_daily_filtered[pt,:].-mean(ha_daily_filtered))./std(ha_daily_filtered)
        tmp_u[pt] = var(u_daily_filtered[pt,:]);
        tmp_ha[pt] = var(ha_daily_filtered[pt,:]);
    end
    
    #u_daily_filtered_normalized = tmp_u;
    #ha_daily_filtered_normalized = tmp_ha;
    
    u_daily_filtered_normalized = u_daily_filtered ./ sqrt(mean(tmp_u));
    ha_daily_filtered_normalized = ha_daily_filtered ./ sqrt(mean(tmp_ha));
    
    # step 2: compute multivariate PCA of data
    state = vcat(u_daily_filtered_normalized, -ha_daily_filtered_normalized) # ~ wind , OLR
    
    if background_type in ["M", "H", "WP"]
        # 2.1 compute corelation matrix
        CCM = state*transpose(state);#/Nd;
        
        # 2.2 compute SVD
        F = svd(CCM);
        V, EV, U = F; 
        #U = transpose(U);
        #U = eigvecs(CCM);
        #EV = eigvals(CCM);
        
        
        # 2.3 first two principal components
        PC1 = -U[:,1];
        PC2 = U[:,2];
        
        # 2.4 variance explained
        var1 = EV[1]/sum(EV);
        var2 = EV[2]/sum(EV);
        
        # the first two PC are:
        PC1u = PC1[1:64];
        PC1ha = PC1[65:end];
        PC2u = PC2[1:64];
        PC2ha = PC2[65:end];
        
        # save the principal components
        file_name = "./PCs_"*background_type*".jld";
        save(file_name, "PC1", PC1, "PC2", PC2);
    
        
    else
        PC1 = load("./PCs_M.jld", "PC1");
        PC2 = load("./PCs_M.jld", "PC2");
        # the first two PC are:
        PC1u = PC1[1:64];
        PC1ha = PC1[65:end];
        PC2u = PC2[1:64];
        PC2ha = PC2[65:end];  
        
    end    
    #state = vcat(u_daily_filtered, -ha_daily_filtered) # ~ wind , OLR
    
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
    if background_type in ["M"]
        return var1, var2, PC1, PC2, amplitude, phase, SMM1, SMM2
    else
        return PC1, PC2, amplitude, phase, SMM1, SMM2
    end
end



end


