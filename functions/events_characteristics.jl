module events


#using LinearAlgebra
#using Statistics
#using FFTW
#using PyPlot

function circleShape(h,k,r)
    θ = LinRange(0, 2*π, 500);
    h .+ r*sin.(θ), k .+ r*cos.(θ)
end




function find_zero_and_nonzero_sequences(arr)
    """ 
        Utility function for the function 'identify'.

        Detect contiguous subsequences of zero and non-zero elements in the array arr.
        Returns 
        - zero_sequences: a list of tuples representing the start and end indices of sequences where the array contains consecutive zeros
        - nonzero_sequences: a list of tuples representing the start and end indices of sequences where the array contains consecutive non-zero values

        (Zero sequences identify 'backward'/ westward movement of the MJO, while non-zero sequences identify 'forward'/ eastward movement of the MJO.)
    """
    zero_sequences = []
    nonzero_sequences = []
    start_index = 1
    in_zero_sequence = false
    in_nonzero_sequence = false

    for (index, value) in enumerate(arr)
        if value == 0
            if in_nonzero_sequence
                push!(nonzero_sequences, (start_index, index - 1))
                in_nonzero_sequence = false
            end

            if !in_zero_sequence
                start_index = index
                in_zero_sequence = true
            end
        else
            if in_zero_sequence
                push!(zero_sequences, (start_index, index - 1))
                in_zero_sequence = false
            end

            if !in_nonzero_sequence
                start_index = index
                in_nonzero_sequence = true
            end
        end
    end

    if in_zero_sequence
        push!(zero_sequences, (start_index, length(arr)))
    end

    if in_nonzero_sequence
        push!(nonzero_sequences, (start_index, length(arr)))
    end

    return zero_sequences, nonzero_sequences
end


function adjust_sequences(zero_sequences, nonzero_sequences)
    """ 
        Utility function for the function 'identify'.
        
        Adjust zero and non-zero sequences. If the length of a zero sequence is smaller than 4 (i.e. MJO is moving backward for less than 4 consecutive days), the sequence is reclassified as a non-zero sequence.
    """

    adjusted_zero_sequences = Tuple{Int, Int}[];
    adjusted_nonzero_sequences = nonzero_sequences;

    for (start, stop) in zero_sequences
        if (stop - start + 1) > 4
            push!(adjusted_zero_sequences, (start, stop))
        else
            push!(adjusted_nonzero_sequences, (start, stop))
        end
    end

    # Initialize an empty array for the modified indices
    modified_indices = []
    sorted_adjusted_nonzero_sequences = sort(adjusted_nonzero_sequences)
    # Iterate through the sorted indices and merge consecutive sequences
    start, stop = sorted_adjusted_nonzero_sequences[1]
    for i in 2:length(sorted_adjusted_nonzero_sequences)
        current_start, current_stop = sorted_adjusted_nonzero_sequences[i]
        if current_start == stop + 1
            stop = current_stop
        else
            push!(modified_indices, (start, stop))
            start, stop = current_start, current_stop
        end
    end
    
    # Add the last modified sequence
    push!(modified_indices, (start, stop))
    
    

    return adjusted_zero_sequences, modified_indices
end



function identify(amplitude, phase)
    """ 
        Given the amplitude and phase of the SMM index, identify MJO events according to the criteria described in [Stachnik et al., 2015](https://doi.org/10.1002/2015JD023916)
        Returns a sequence of tuples representing the start and end indices (simulation days) of each event.
    """

    # step 1: identify days with sufficient amplitude 
    indices1 = findall(x->x>=1, amplitude); # days where amplitude is >=1
    
    ## step 2: constrain this to the days which are not part of a triplet of backward steps (ie westward movement)
    #indices1 = []; 
    #steps = phase[2:end]-phase[1:end-1];
    #steps = (sort(hcat(steps, steps.+2*pi, steps.-2*pi), by = abs, dims=2))[:,1]; # difference on the circle
##
    #for i in indices11[3:end-2]
    #    if (steps[i]<0) # if it is a backward step we have to check that it is not part of a triplet of backward step, if not: store it
    #        if !( ((steps[i+1]<0) & (steps[i+2]<0)) | ((steps[i-1]<0) & (steps[i+1]<0)) | ((steps[i-2]<0) & (steps[i-1]<0)) ) 
    #            append!(indices1,[i])
    #        end
    #    else # if it is a forward step: store it
    #        append!(indices1,[i])
    #    end
    #end
    #
    ## boundaries
    #if steps[1]>=0
    #    append!([indices11[1]], indices1)
    #end
    #if steps[2]>=0
    #    append!([indices11[2]], indices1)
    #end
    #if steps[end-1]>=0
    #    append!(indices1, [indices11[end-1]])
    #end
    #if steps[end]>=0
    #    append!(indices1, [indices11[end]])
    #end
    
    
    # step 3: find continuous sequences in the array of days where amplitude is >=1
    diff = indices1[2:end]-indices1[1:end-1]; 
    indices2 = findall(x->x>2, diff); # find "cuts" 
    num_events = length(indices2)+1; # number of events is the number of cuts + 1
    
    # step 4: find the initial and final day of each event (events of amplitude >1 and eastward movement)
    seq = zeros((num_events,2));
    seq[1,1]=indices1[1]trailer;
    i=1;trailer
    for cut in indices2[1:end]
        seq[i,2] = indices1[cut];
        seq[i+1,1] = indices1[indices2[i]+1];
        i+=1;
    end
    seq[end,2] = indices1[end]; 
    
    
    # step 5: only keep events which propagate through 4 phases of the phase diagram at least
    num_el_deleted = 0;
    final_seq = copy(seq);
    for i in 1:num_events
        d1 = Int(seq[i,1]);
        d2 = Int(seq[i,2]);
        steps = (phase[d1+1:d2] .- phase[d1:d2-1]);
        diffs = (sort(hcat(steps, steps.+2*pi, steps.-2*pi), by = abs, dims=2))[:,1];
        if sum(diffs)< 3*π/4
            # delete event from matrix
            final_seq = final_seq[1:end .!= i-num_el_deleted,:]
            num_el_deleted += 1;
        end
    end
    
    
    # step 6: check that no event propagates westward more than one phase
    final_final_seq = [];
    for fsi in 1:size(final_seq)[1]
        fs = final_seq[fsi,:]
        d1 = Int(fs[1]);
        d2 = Int(fs[2]);
        progression = phase[d1:d2];
        sub_seq = [1]; # log of indices to track forward/backward motion (backward motion indicated by 0s) 
        for (i, x) in enumerate(progression[1:end-1])
            if i>1
                if((x<=π/4.0)&(x>0)) # phase 5
                    x1 = progression[i+1];
                    if !((x1<=π/2.0)&(x1>0)) # if we are in phase 5, the next phase must be either 6 or 5
                        j = i;
                        # test the previous step -> aim is to identify the first where the RMM began to retrograde (must find first clockwise steps in phase 5)
                        step = progression[j] - progression[j-1];
                        diff = (sort(hcat(step, step+2*pi, step-2*pi), by = abs, dims=2))[:,1][1];
                        count=0;
                        while (j>2) & (progression[j-1]<=π/4.0) & (progression[j-1]>0.0) &  (diff<=0)  # j>2 & still in phase 5 & clockwise progression  
                            count+=1;
                            j=j-1
                            step = progression[j] - progression[j-1];
                            diff = (sort(hcat(step, step+2*pi, step-2*pi), by = abs, dims=2))[:,1][1];
                        end # j should be the index of the last clockwise element
                        sub_seq[end-count+1:end].=0; # retrograde steps are zeroed
                        push!(sub_seq,0)
                    else
                        push!(sub_seq, i)
                    end
                    
                elseif((x<=π/2.0)&(x>π/4.0)) # phase 6
                    x1 = progression[i+1];
                    if !((x1<=3*π/4.0)&(x1>π/4.0)) # if we are in phase 6, the next phase must be either 7 or 6
                        j = i;
                        # test the previous step -> aim is to identify the first where the RMM began to retrograde (must find first clockwise steps in phase 6)
                        step = progression[j] - progression[j-1];
                        diff = (sort(hcat(step, step+2*pi, step-2*pi), by = abs, dims=2))[:,1][1];
                        count=0;
                        while (j>2) & (progression[j-1]<=π/2.0) & (progression[j-1]>π/4.0) & (diff<=0)   # j>1 & still in phase 6 & clockwise progression  
                            count+=1;
                            j=j-1
                            step = progression[j] - progression[j-1];
                            diff = (sort(hcat(step, step+2*pi, step-2*pi), by = abs, dims=2))[:,1][1];
                        end
                        sub_seq[end-count+1:end].=0;
                        push!(sub_seq,0)
                    else
                        push!(sub_seq, i)
                    end
                    
                elseif((x<=3*π/4.0)&(x>π/2.0)) # phase 7
                    x1 = progression[i+1];
                    if !((x1<=π)&(x1>π/2.0)) # if we are in phase 7, the next phase must be either 8 or 7
                        j = i;
                        # test the previous step -> aim is to identify the first where the RMM began to retrograde (must find first clockwise steps in phase 7)
                        step = progression[j] - progression[j-1];
                        diff = (sort(hcat(step, step+2*pi, step-2*pi), by = abs, dims=2))[:,1][1];
                        count=0;
                        while (j>2) & (progression[j-1]<=3*π/4.0) & (progression[j-1]>π/2.0) & (diff<=0)   # j>1 & still in phase 7 & clockwise progression  
                            count+=1;
                            j=j-1
                            step = progression[j] - progression[j-1];
                            diff = (sort(hcat(step, step+2*pi, step-2*pi), by = abs, dims=2))[:,1][1];
                        end
                        sub_seq[end-count+1:end].=0;
                        push!(sub_seq,0) 
                    else
                        push!(sub_seq, i)
                    end
                    
                elseif((x<=π)&(x>3*π/4.0)) # phase 8
                    x1 = progression[i+1]; 
                    if !(((x1<=π)&(x1>3*π/4.0)) | ((x1<=-3*π/4.0)&(x1>-π))) # if we are in phase 8, the next phase must be either 1 or 8
                        j = i;
                        # test the previous step -> aim is to identify the first where the RMM began to retrograde (must find first clockwise steps in phase 8)
                        step = progression[j] - progression[j-1];
                        diff = (sort(hcat(step, step+2*pi, step-2*pi), by = abs, dims=2))[:,1][1];
                        count=0;
                        while (j>2) & (progression[j-1]<=π) & (progression[j-1]>3*π/4.0)  & (diff<=0)   # j>1 & still in phase 8 & clockwise progression  
                            count+=1;
                            j=j-1
                            step = progression[j] - progression[j-1];
                            diff = (sort(hcat(step, step+2*pi, step-2*pi), by = abs, dims=2))[:,1][1];
                        end
                        sub_seq[end-count+1:end].=0;
                        push!(sub_seq,0)
                    else
                        push!(sub_seq, i)
                    end
    
                elseif((x<=-3*π/4.0)&(x>-π))
                    x1 = progression[i+1];
                    if !((x1<=-π/2.0)&(x1>-π))
                        j = i;
                        # test the previous step -> aim is to identify the first where the RMM began to retrograde 
                        step = progression[j] - progression[j-1];
                        diff = (sort(hcat(step, step+2*pi, step-2*pi), by = abs, dims=2))[:,1][1];
                        count=0;
                        while (j>2) & (progression[j-1]<=-3*π/4.0) & (progression[j-1]>-π) & (diff<=0)   #   
                            count+=1;
                            j=j-1
                            step = progression[j] - progression[j-1];
                            diff = (sort(hcat(step, step+2*pi, step-2*pi), by = abs, dims=2))[:,1][1];
                        end
                        sub_seq[end-count+1:end].=0;
                        push!(sub_seq,0)
                    else
                        push!(sub_seq, i)
                    end
                    
                elseif((x<=-π/2.0)&(x>-3*π/4.0))
                    x1 = progression[i+1];
                    if !((x1<=-π/4.0)&(x1>-3*π/4.0))
                        j = i;
                        # test the previous step -> aim is to identify the first where the RMM began to retrograde 
                        step = progression[j] - progression[j-1];
                        diff = (sort(hcat(step, step+2*pi, step-2*pi), by = abs, dims=2))[:,1][1];
                        count=0;
                        while (j>2) & (progression[j-1]<=-π/2.0) & (progression[j-1]>-3*π/4.0) &  (diff<=0)   #   
                            count+=1;
                            j=j-1
                            step = progression[j] - progression[j-1];
                            diff = (sort(hcat(step, step+2*pi, step-2*pi), by = abs, dims=2))[:,1][1];
                        end
                        sub_seq[end-count+1:end].=0;
                        push!(sub_seq,0)
                    else
                        push!(sub_seq, i)
                    end
                     
                elseif((x<=-π/4.0)&(x>-π/2.0))
                    x1 = progression[i+1];
                    if !((x1<=0)&(x1>-π/2.0))
                        j = i;
                        # test the previous step -> aim is to identify the first where the RMM began to retrograde 
                        step = progression[j] - progression[j-1];
                        diff = (sort(hcat(step, step+2*pi, step-2*pi), by = abs, dims=2))[:,1][1];
                        count=0;
                        while (j>2) & (progression[j-1]<=-π/4.0) & (progression[j-1]>-π/2.0) &  (diff<=0)   #   
                            count+=1;
                            j=j-1
                            step = progression[j] - progression[j-1];
                            diff = (sort(hcat(step, step+2*pi, step-2*pi), by = abs, dims=2))[:,1][1];
                        end
                        sub_seq[end-count+1:end].=0;
                        push!(sub_seq,0)
                    else
                        push!(sub_seq, i)
                    end
                    
                elseif((x<=0)&(x>-π/4.0)) # phase 4
                    x1 = progression[i+1];
                    if !((x1<=π/4.0)&(x1>-π/4.0)) # next must be phase 4 or 5
                        j = i;
                        # test the previous step -> aim is to identify the first where the RMM began to retrograde (must find first clockwise steps in phase 4)
                        step = progression[j] - progression[j-1];
                        diff = (sort(hcat(step, step+2*pi, step-2*pi), by = abs, dims=2))[:,1][1];
                        count=0;
                        while (j>2) & (progression[j-1]<=0) & (progression[j-1]>-π/4.0) &  (diff<=0)   # j>1 & still in phase 4 & clockwise progression  
                            count+=1;
                            j=j-1
                            step = progression[j] - progression[j-1];
                            diff = (sort(hcat(step, step+2*pi, step-2*pi), by = abs, dims=2))[:,1][1];
                        end
                        sub_seq[end-count+1:end].=0;
                        push!(sub_seq,0) 
                    else
                        push!(sub_seq, i)
                    end
                end
            end
        end
        push!(sub_seq,d2-d1+1)
        if 0 in sub_seq 
            # find all continuous sequences of zeros (indicating retrograde motion)
            zero_sequences, nonzero_sequences = find_zero_and_nonzero_sequences(sub_seq);
            # adjust so that only 4 steps are taken backward 
            adjusted_zero_sequences, adjusted_nonzero_sequences = adjust_sequences(zero_sequences, nonzero_sequences)
            days_in_original_seq = collect(d1:d2)
            
            if length(adjusted_nonzero_sequences)!=0
                for s in adjusted_nonzero_sequences
                    i1 = days_in_original_seq[sub_seq[s[1]]]
                    i2 = days_in_original_seq[sub_seq[s[2]]]
                    if (i2-i1)>10 # consider only events lasting more than 10 days
                        steps = (phase[i1+1:i2] .- phase[i1:i2-1]);
                        diffs = (sort(hcat(steps, steps.+2*pi, steps.-2*pi), by = abs, dims=2))[:,1];
                        if sum(diffs)> 3*π/4
                            push!(final_final_seq, (i1,i2));
                        end
                    end
                end
            else
                push!(final_final_seq, (d1,d2));
                #println("append : ", i1:i2)
            end
        else
            push!(final_final_seq, (d1,d2));
        end
    end
    
    num_events = size(final_final_seq)[1];
    final_seq = zeros(num_events,2)
    for i in 1:num_events
        final_seq[i,1] = final_final_seq[i][1]
        final_seq[i,2] = final_final_seq[i][2]
    end
        
    return final_seq;
end




function div_seq_ENSO_phases(seq, years_tracked, months_tracked, ENSO_ym, ENSO_phases)

    """
        Classifies MJO events based on the corresponding ENSO phases. 
        If an MJO event occurs over more than one ENSO phase, the dominant phase is chosen. 

        ### Inputs:
        - 'seq': sequence of tupples (start_idx, end_idx) representing the start and end indices of MJO events.
        - 'years_tracked', 'months_tracked': lists mapping each index to a year, month
        - 'ENSO_ym': list of dates for which ENSO phases are tracked
        - 'ENSO_phases': list of ENSO phases corresponding to the dates in 'ENSO_ym'

        ### Outputs:
        - 'seqN': Subset of MJO events that occurred predominantly during Neutral ENSO phases.
        - 'seqEN': Subset of MJO events that occurred predominantly during El Niño phases.
        - 'seqLN': Subset of MJO events that occurred predominantly during La Niña phases.
        - 'frac': fraction of time each MJO event spent in Neutral, El Niño, and La Niña phases, respectively. 

    """

    sep = zeros(size(seq)[1]); # contains the phase of ENSO for each event (as -1/0/+1 for LN/N/EN)
    frac = zeros(size(seq)[1], 3); 
    # for each MJO event
    for i in 1:size(seq)[1]
        y1 = years_tracked[Int(seq[i,1])];
        y2 = years_tracked[Int(seq[i,2])];
        m1 = months_tracked[Int(seq[i,1])];
        m2 = months_tracked[Int(seq[i,2])];
        
        idx_idx1 = findall(all(ENSO_ym[:,:] .== [y1 m1], dims=2));
        idx_idx2 = findall(all(ENSO_ym[:,:] .== [y2 m2], dims=2));

        sub_phases = ENSO_phases[idx_idx1[1][1]:idx_idx2[1][1]];
        
        nN = length(findall(x->x=="N", sub_phases));
        nEN = length(findall(x->x=="EN", sub_phases));
        nLN = length(findall(x->x=="LN", sub_phases));
        if (nN >= nEN) & (nN >= nLN) # mainly neutral
            sep[i] = 0;
        elseif (nEN > nN) & (nEN > nLN)
            sep[i] = 1;
        elseif (nLN > nN) & (nLN > nEN)
            sep[i] = -1;
        end
        frac[i,:] = [nN/length(sub_phases), nEN/length(sub_phases), nLN/length(sub_phases)];
    end
    
    seqN = seq[findall(x->x==0, sep),:];
    seqEN = seq[findall(x->x==1, sep),:];
    seqLN = seq[findall(x->x==-1, sep),:];
    
    return seqN, seqEN, seqLN, frac        
end





function number(seq)
    """ 
        Inputs: 
        - 'seq': sequence of tupples (start_idx, end_idx) representing the start and end indices of MJO events
        Outputs:
        - total number of MJO events
    """
    return size(seq)[1];
end


function seasonality(seq, months_tracked)
    """ 
        Count number of MJO events occuring during each month of the year. 
    """
    seasonal_var_local = zeros(12);
    num_events = size(seq)[1];
    for i in 1:num_events
        m1 = Int(seq[i,1]);
        m2 = Int(seq[i,2]);
        event_months = months_tracked[m1:m2];
        for em in unique(event_months)
            seasonal_var_local[Int(em)]+=1;
        end
    end
    return seasonal_var_local
end


function total_angle(seq, phase)
    """ 
        Inputs:
        - 'seq': sequence of tupples (start_idx, end_idx) representing the start and end indices of MJO events
        - 'phase': phase of the Skeleton Multivariate MJO index for each simulation day
        Output: total angle covered by each MJO event in the (SMM1,SMM2) diagram.
    """
    num_events = size(seq)[1];
    distances = zeros(num_events); 
    for i in 1:num_events;
        
        d1 = Int(seq[i,1]);
        d2 = Int(seq[i,2]);
        
        steps = (phase[d1+1:d2] .- phase[d1:d2-1]);
        diffs = (sort(hcat(steps, steps.+2*pi, steps.-2*pi), by = abs, dims=2))[:,1];
        distances[i] = sum(diffs);
    end
    return distances    
end


function duration(seq)
    """ 
        Output: lifetime (in days) of each MJO event in the sequence.
    """
    durations = zeros(size(seq)[1]);
    for i in 1:size(seq)[1]
        durations[i] = Int(seq[i,2])-Int(seq[i,1])+1;
    end
    return durations;
end

function max_amplitude(seq, amplitude, phase)
    """ 
        Output: maximum amplitude of the SMM index and its corresponding phase for each MJO event in the sequence.
    """
    max_amplitude = zeros(size(seq)[1]);
    max_ampl_location = zeros(size(seq)[1]); 
    
    
    for i in 1:size(seq)[1]
        max_amplitude[i], maxidx = findmax(amplitude[Int(seq[i,1]):Int(seq[i,2])]);
        max_ampl_location[i] = phase[Int(seq[i,1]):Int(seq[i,2])][maxidx]; 
    end
    
    return max_amplitude, max_ampl_location;
end



function phase_space_diagram(SMM1, SMM2, ii1, ii2, amplitude, f1, f2, jj1, jj2, L1, L2, ft)

    """ 
        Plot (SMM1,SMM2) diagram, illustrating the MJO events.
    """    
    fig, axs = subplots(1,1, figsize = (6,6))
    
    #plot a circle of radius 1
    Circle = circleShape(0,0,1);
    axs.plot(Circle[1], Circle[2], color = "black", linewidth = 1)
    
    # plot phases
    axs.plot(cos(π/4):0.01:4, sin(π/4):0.01:4, "--", linewidth =1, color = "black")
    axs.plot(-4:0.01:cos(5*π/4), -4:0.01:sin(5*π/4), "--", linewidth =1, color = "black")
    axs.plot(-4:0.01:cos(3*π/4), 4:-0.01:sin(3*π/4), "--", linewidth =1, color = "black")
    axs.plot(cos(-π/4):0.01:4, sin(-π/4):-0.01:-4, "--", linewidth =1, color = "black")
    axs.plot(zeros(length(1:0.01:4)), 1:0.01:4, "--", linewidth =1, color = "black")
    axs.plot(zeros(length(1:0.01:4)), -1:-0.01:-4, "--", linewidth =1, color = "black")
    axs.plot(-4:0.01:-1, zeros(length(1:0.01:4)), "--", linewidth =1, color = "black")
    axs.plot(1:0.01:4, zeros(length(1:0.01:4)), "--", linewidth =1, color = "black")
    
    
    text(-1.4, -2.8, "2", fontsize=14, color = "grey")
    text(1.4, -2.8, "3", fontsize=14, color = "grey")
    text(-0.7, -3.9, "Indian Ocean", fontsize=12, color = "grey")
    text(2.6, -1, "4", fontsize=14, color = "grey")
    text(2.6, 1, "5", fontsize=14, color = "grey")
    text(3.5, -1, "Maritime Continent", fontsize=12, color = "grey", rotation="vertical")
    text(1.4, 2.6, "6", fontsize=14, color = "grey")
    text(-1.4, 2.6, "7", fontsize=14, color = "grey")
    text(-0.7, 3.2, "Pacific Ocean", fontsize=12, color = "grey")
    text(-2.7, 1, "8", fontsize=14, color = "grey")
    text(-2.7, -1, "1", fontsize=14, color = "grey")
    text(-4.0, -1, "West. Hem. & Africa", fontsize=12, color = "grey", rotation="vertical")
    
    
    axs.set_xlim([-3,3])
    axs.set_ylim([-3,3])
    
    # plot event
    axs.scatter(SMM1[ii1:jj1], SMM2[ii1:jj1], s=30, color = "skyblue", alpha = 0.6)
    axs.scatter(SMM1[jj2:ii2], SMM2[jj2:ii2], s=30, color = "skyblue", alpha = 0.6)
    axs.scatter(SMM1[jj1:jj2], SMM2[jj1:jj2], s=30, color = "steelblue", alpha = 1)
    
    axs.text(SMM1[f1], SMM2[f1], L1, fontsize=ft)
    axs.text(SMM1[f2], SMM2[f2], L2, fontsize=ft)
   
    #axs.scatter(SMM1[f1], SMM2[f1], s=30, color = "tomato")
    #axs.scatter(SMM1[f2], SMM2[f2], s=30, color = "black")
    #axs.scatter(SMM1[ii1+1], SMM2[ii1+1], s=12, color = "yellow")
    #axs.scatter(SMM1[ii1+2], SMM2[ii1+2], s=12, color = "blue")

    #m, j = findmax(amplitude[ii1:ii2]);
    #axs.scatter(SMM1[ii1-1+j], SMM2[ii1-1+j], s=12, color = "red")
    
    axs.set_xlabel("RMM1", fontsize = 10)
    axs.set_ylabel("RMM2", fontsize = 10)
    
    savefig("RMM.png", dpi = 400, bbox_inches="tight")
end





function starting_phase(seq, phase)
    """ 
        Output: number of events starting in each phase
    """

    sp = zeros(8);
    num_events = size(seq)[1];
    for i in 1:num_events
        d1 = Int(seq[i,1]);
        p = phase[d1];
        if (p<=π/4.0)&(p>0) # phase 5 
            sp[5] += 1; 
        elseif (p<=π/2.0)&(p>π/4.0) # phase 6
            sp[6] += 1;
        elseif (p<=3*π/4.0)&(p>π/2.0) # phase 7
            sp[7] += 1; 
        elseif (p<=π)&(p>3*π/4.0) # phase 8
            sp[8] += 1; 
        elseif (p<=-3*π/4.0)&(p>-π) # phase 1
            sp[1] += 1; 
        elseif (p<=-π/2.0)&(p>-3*π/4.0) # phase 2 
            sp[2] += 1; 
        elseif (p<=-π/4.0)&(p>-π/2.0) # phase 3
            sp[3] += 1; 
        elseif (p<=0)&(p>-π/4.0) # phase 4 
            sp[4] += 1; 
        end
    end
    return sp
end

function ending_phase(seq, phase)
    """ 
        Output: number of events stoping in each phase
    """
    
    ep = zeros(8);
    num_events = size(seq)[1];
    for i in 1:num_events
        d2 = Int(seq[i,2]);
        p = phase[d2];
        if (p<=π/4.0)&(p>0) # phase 5 
            ep[5] += 1; 
        elseif (p<=π/2.0)&(p>π/4.0) # phase 6
            ep[6] += 1;
        elseif (p<=3*π/4.0)&(p>π/2.0) # phase 7
            ep[7] += 1; 
        elseif (p<=π)&(p>3*π/4.0) # phase 8
            ep[8] += 1; 
        elseif (p<=-3*π/4.0)&(p>-π) # phase 1
            ep[1] += 1; 
        elseif (p<=-π/2.0)&(p>-3*π/4.0) # phase 2 
            ep[2] += 1; 
        elseif (p<=-π/4.0)&(p>-π/2.0) # phase 3
            ep[3] += 1; 
        elseif (p<=0)&(p>-π/4.0) # phase 4 
            ep[4] += 1; 
        end
    end
    return ep
end



end