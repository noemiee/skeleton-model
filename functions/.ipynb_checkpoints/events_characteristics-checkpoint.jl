module events

#export identify
using PyPlot
using FFTW



function circleShape(h,k,r)
    θ = LinRange(0, 2*π, 500);
    h .+ r*sin.(θ), k .+ r*cos.(θ)
end


function consecutive_groups(array)
    groups = []
    j = 0
    for i=1:length(array)-1
        if array[i]+1 != array[i+1]
            push!(groups, array[j+1:i])
            j = i
        end
    end
    push!(groups, array[j+1:end])
    return groups;
end


function groups(phase_name, phase_values)
   
    if phase_name == "p5"
        p5 = findall(x->((x<=π/4.0)&(x>0)), phase_values);
        return consecutive_groups(p5);
    elseif phase_name == "p6"
        p6 = findall(x->((x<=π/2.0)&(x>π/4.0)), phase_values);
        return consecutive_groups(p6);
    elseif phase_name == "p7"
        p7 = findall(x->((x<=3*π/4.0)&(x>π/2.0)), phase_values);
        return consecutive_groups(p7);
    elseif phase_name == "p8"
        p8 = findall(x->((x<=π)&(x>3*π/4.0)), phase_values);
        return consecutive_groups(p8);
    elseif phase_name == "p1"
        p1 = findall(x->((x<=-3*π/4.0)&(x>-π)), phase_values);
        return consecutive_groups(p1);
    elseif phase_name == "p2"
        p2 = findall(x->((x<=-π/2.0)&(x>-3*π/4.0)), phase_values);
        return consecutive_groups(p2);
    elseif phase_name == "p3"
        p3 = findall(x->((x<=-π/4.0)&(x>-π/2.0)), phase_values);
        return consecutive_groups(p3);
    elseif phase_name == "p4"
        p4 = findall(x->((x<=0)&(x>-π/4.0)), phase_values);
        return consecutive_groups(p4);
    end  
    
end


function find_zero_and_nonzero_sequences(arr)
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
    
    
    # step 3: find continuous sequences in the array of days
    diff = indices1[2:end]-indices1[1:end-1]; 
    indices2 = findall(x->x>2, diff); # find "cuts" 
    num_events = length(indices2)+1; # number of events is the number of cuts + 1
    
    # step 4: find the initial and final day of each event (events of amplitude >1 and eastward movement)
    seq = zeros((num_events,2));
    seq[1,1]=indices1[1];
    i=1;
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
        sub_seq = [1];
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
                        sub_seq[end-count+1:end].=0;
                        push!(sub_seq,0)
                    else
                        push!(sub_seq, i)
                    end
                    
                elseif((x<=π/2.0)&(x>π/4.0))
                    x1 = progression[i+1];
                    if !((x1<=3*π/4.0)&(x1>π/4.0)) # if we are in phase 6, the next phase must be either 7 or 6
                        j = i;
                        # test the previous step -> aim is to identify the first where the RMM began to retrograde (must find first clockwise steps in phase 5)
                        step = progression[j] - progression[j-1];
                        diff = (sort(hcat(step, step+2*pi, step-2*pi), by = abs, dims=2))[:,1][1];
                        count=0;
                        while (j>2) & (progression[j-1]<=π/2.0) & (progression[j-1]>π/4.0) & (diff<=0)   # j>1 & still in phase 5 & clockwise progression  
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
                    
                elseif((x<=3*π/4.0)&(x>π/2.0))
                    x1 = progression[i+1];
                    if !((x1<=π)&(x1>π/2.0)) # if we are in phase 7, the next phase must be either 8 or 7
                        j = i;
                        # test the previous step -> aim is to identify the first where the RMM began to retrograde (must find first clockwise steps in phase 5)
                        step = progression[j] - progression[j-1];
                        diff = (sort(hcat(step, step+2*pi, step-2*pi), by = abs, dims=2))[:,1][1];
                        count=0;
                        while (j>2) & (progression[j-1]<=3*π/4.0) & (progression[j-1]>π/2.0) & (diff<=0)   # j>1 & still in phase 5 & clockwise progression  
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
                    
                elseif((x<=π)&(x>3*π/4.0))
                    x1 = progression[i+1]; 
                    if !(((x1<=π)&(x1>3*π/4.0)) | ((x1<=-3*π/4.0)&(x1>-π))) # if we are in phase 8, the next phase must be either 1 or 8
                        j = i;
                        # test the previous step -> aim is to identify the first where the RMM began to retrograde (must find first clockwise steps in phase 5)
                        step = progression[j] - progression[j-1];
                        diff = (sort(hcat(step, step+2*pi, step-2*pi), by = abs, dims=2))[:,1][1];
                        count=0;
                        while (j>2) & (progression[j-1]<=π) & (progression[j-1]>3*π/4.0)  & (diff<=0)   # j>1 & still in phase 5 & clockwise progression  
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
                        # test the previous step -> aim is to identify the first where the RMM began to retrograde (must find first clockwise steps in phase 5)
                        step = progression[j] - progression[j-1];
                        diff = (sort(hcat(step, step+2*pi, step-2*pi), by = abs, dims=2))[:,1][1];
                        count=0;
                        while (j>2) & (progression[j-1]<=-3*π/4.0) & (progression[j-1]>-π) & (diff<=0)   # j>1 & still in phase 5 & clockwise progression  
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
                        # test the previous step -> aim is to identify the first where the RMM began to retrograde (must find first clockwise steps in phase 5)
                        step = progression[j] - progression[j-1];
                        diff = (sort(hcat(step, step+2*pi, step-2*pi), by = abs, dims=2))[:,1][1];
                        count=0;
                        while (j>2) & (progression[j-1]<=-π/2.0) & (progression[j-1]>-3*π/4.0) &  (diff<=0)   # j>1 & still in phase 5 & clockwise progression  
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
                        # test the previous step -> aim is to identify the first where the RMM began to retrograde (must find first clockwise steps in phase 5)
                        step = progression[j] - progression[j-1];
                        diff = (sort(hcat(step, step+2*pi, step-2*pi), by = abs, dims=2))[:,1][1];
                        count=0;
                        while (j>2) & (progression[j-1]<=-π/4.0) & (progression[j-1]>-π/2.0) &  (diff<=0)   # j>1 & still in phase 5 & clockwise progression  
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
                        # test the previous step -> aim is to identify the first where the RMM began to retrograde (must find first clockwise steps in phase 5)
                        step = progression[j] - progression[j-1];
                        diff = (sort(hcat(step, step+2*pi, step-2*pi), by = abs, dims=2))[:,1][1];
                        count=0;
                        while (j>2) & (progression[j-1]<=0) & (progression[j-1]>-π/4.0) &  (diff<=0)   # j>1 & still in phase 5 & clockwise progression  
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
        if 0 in sub_seq # something has to be done
            #println("original: ", d1:d2)
            #println(sub_seq)
            # find all continuous sequences of zeros
            zero_sequences, nonzero_sequences = find_zero_and_nonzero_sequences(sub_seq);
            # adjust so that if only 4 steps are taken backward 
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
                            #println("append : ", i1:i2)
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




#function div_seq_ENSO_phases(seq, years_tracked, months_tracked, ENSO_ym, ENSO_phases)
#    sep = zeros(size(seq)[1]); # contains the phase of ENSO for each event (as -1/0/+1 for LN/N/EN)
#    frac = zeros(size(seq)[1], 3); 
#    # for each MJO event
#    for i in 1:size(seq)[1]
#        y1 = years_tracked[Int(seq[i,1])];
#        y2 = years_tracked[Int(seq[i,2])];
#        m1 = months_tracked[Int(seq[i,1])];
#        m2 = months_tracked[Int(seq[i,2])];
#        
#        idx_idx1 = findall(all(ENSO_ym[:,:] .== [y1 m1], dims=2));
#        idx_idx2 = findall(all(ENSO_ym[:,:] .== [y2 m2], dims=2));
#
#        sub_phases = ENSO_phases[idx_idx1[1][1]:idx_idx2[1][1]];
#        
#        nN = length(findall(x->x=="N", sub_phases));
#        nEN = length(findall(x->x=="EN", sub_phases));
#        nLN = length(findall(x->x=="LN", sub_phases));
#        if (nN >= nEN) & (nN >= nLN) # mainly neutral
#            sep[i] = 0;
#        elseif (nEN > nN) & (nEN > nLN)
#            sep[i] = 1;
#        elseif (nLN > nN) & (nLN > nEN)
#            sep[i] = -1;
#        else
#            sep[i] = 1000;
#        end
#        frac[i,:] = [nN/length(sub_phases), nEN/length(sub_phases), nLN/length(sub_phases)];
#    end
#    
#    seqN = seq[findall(x->x==0, sep),:];
#    seqEN = seq[findall(x->x==1, sep),:];
#    seqLN = seq[findall(x->x==-1, sep),:];
#    
#    return seqN, seqEN, seqLN, frac
#            
#end


function div_seq_ENSO_phases(seq, years_tracked, months_tracked, ENSO_ym, ENSO_phases)
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

function div_seq_PDO_phases(seq, sim_day, PDO_phases_monthly)
    sep = zeros(size(seq)[1]);
    frac = zeros(size(seq)[1], 2);
    # for each MJO event
    for i in 1:size(seq)[1]
        d1 = sim_day[Int(seq[i,1])];
        d2 = sim_day[Int(seq[i,2])];
        m1 = floor(Int, d1/30)+1;
        m2 = floor(Int, d2/30)+1;
        sub_phases = PDO_phases_monthly[m1:m2]
        nPPDO = length(findall(x->x=="PPDO", sub_phases));
        nNPDO = length(findall(x->x=="NPDO", sub_phases));
        if (nPPDO > nNPDO)
            sep[i] = 1;
        else
            sep[i] = -1;
        end
        frac[i,:] = [nPPDO/length(sub_phases), nNPDO/length(sub_phases)];
    end
    
    seqPPDO = seq[findall(x->x==1, sep),:];
    seqNPDO = seq[findall(x->x==-1, sep),:];
    
    return seqPPDO, seqNPDO, frac
            
end



function number(seq)
    return size(seq)[1];
end

function total_angle(seq, phase)
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


function num_sub_events(seq, phase)
    num_events = size(seq)[1];
    n_sub_events = zeros(num_events); 
    for i in 1:num_events;
        
        d1 = Int(seq[i,1]);
        d2 = Int(seq[i,2]);
        
        phase_values = phase[d1:d2];
        
        groups = groups_primary_continuing(phase_values);
    
        n_sub_events[i] = size(groups)[1];
    end
    return n_sub_events    
end




function duration(seq)
    durations = zeros(size(seq)[1]);
    for i in 1:size(seq)[1]
        durations[i] = Int(seq[i,2])-Int(seq[i,1])+1;
    end
    return durations;
end

function max_amplitude(seq, amplitude, phase)
    max_amplitude = zeros(size(seq)[1]);
    max_ampl_location = zeros(size(seq)[1]); 
    
    
    for i in 1:size(seq)[1]
        max_amplitude[i], maxidx = findmax(amplitude[Int(seq[i,1]):Int(seq[i,2])]);
        max_ampl_location[i] = phase[Int(seq[i,1]):Int(seq[i,2])][maxidx]; 
    end
    
    return max_amplitude, max_ampl_location;
end

function phase_space_diagram(SMM1, SMM2, ii1, ii2, amplitude, f1, f2, jj1, jj2, L1, L2, ft)
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


function composites(field, seq, phase, phase_name, filter, maxk)
    Nx = size(field)[1];
    nts = size(field)[end]; 
    
    if filter == true
        # filter spaially the subset
        field_spatially_filtered = zeros(size(field))
        for ts in 1:nts # filter each time step separately
            fH = fftshift(fft(field[:,ts]));   
            freqs = fftshift(fftfreq(Nx, Nx));
            idx1 = findall(x-> abs(x) >maxk , freqs);
            fH[idx1].=0;
            
            field_spatially_filtered[:,ts] = real(ifft(ifftshift(fH)));
        end
    field = field_spatially_filtered;
    end
    
    compo = zeros(Nx);
    count = 0;
    for i in 1:size(seq)[1]
        
        if phase_name == "p5"
            p = findall(x->((x<=π/4.0)&(x>0)), phase[Int(seq[i,1]):Int(seq[i,2])]);
        elseif phase_name == "p6"
            p = findall(x->((x<=π/2.0)&(x>π/4.0)), phase[Int(seq[i,1]):Int(seq[i,2])]);
        elseif phase_name == "p7"
            p = findall(x->((x<=3*π/4.0)&(x>π/2.0)), phase[Int(seq[i,1]):Int(seq[i,2])]);
        elseif phase_name == "p8"
            p = findall(x->((x<=π)&(x>3*π/4.0)), phase[Int(seq[i,1]):Int(seq[i,2])]);
        elseif phase_name == "p1"
            p = findall(x->((x<=-3*π/4.0)&(x>-π)), phase[Int(seq[i,1]):Int(seq[i,2])]);
        elseif phase_name == "p2"
            p = findall(x->((x<=-π/2.0)&(x>-3*π/4.0)), phase[Int(seq[i,1]):Int(seq[i,2])]);
        elseif phase_name == "p3"
            p = findall(x->((x<=-π/4.0)&(x>-π/2.0)), phase[Int(seq[i,1]):Int(seq[i,2])]);
        elseif phase_name == "p4"
            p = findall(x->((x<=0)&(x>-π/4.0)), phase[Int(seq[i,1]):Int(seq[i,2])]);
        else
            println("phase not recognized");
        end
        
        for pidx in p
            compo += field[:,Int(seq[i,1]):Int(seq[i,2])][:,pidx];
            count += 1.0;
        end
    end
    println(count)
    return compo./count; 
end

function ndays_per_phase(seq, phase, typ)
    ndpp = zeros(8);
    is = zeros(8);
    for i in 1:size(seq)[1]
        # steps in phase 1
        
        p5 = findall(x->((x<=π/4.0)&(x>0)), phase[Int(seq[i,1]):Int(seq[i,2])]);
        p6 = findall(x->((x<=π/2.0)&(x>π/4.0)), phase[Int(seq[i,1]):Int(seq[i,2])]);
        p7 = findall(x->((x<=3*π/4.0)&(x>π/2.0)), phase[Int(seq[i,1]):Int(seq[i,2])]);
        p8 = findall(x->((x<=π)&(x>3*π/4.0)), phase[Int(seq[i,1]):Int(seq[i,2])]);
        p1 = findall(x->((x<=-3*π/4.0)&(x>-π)), phase[Int(seq[i,1]):Int(seq[i,2])]);
        p2 = findall(x->((x<=-π/2.0)&(x>-3*π/4.0)), phase[Int(seq[i,1]):Int(seq[i,2])]);
        p3 = findall(x->((x<=-π/4.0)&(x>-π/2.0)), phase[Int(seq[i,1]):Int(seq[i,2])]);
        p4 = findall(x->((x<=0)&(x>-π/4.0)), phase[Int(seq[i,1]):Int(seq[i,2])]);
        
        for (i, pp) in enumerate([p1,p2,p3,p4,p5,p6,p7,p8])
            if length(pp) !=0
                is[i]+=1;
                ndpp[i]+=length(pp);
            end
        end
    
    end
    if typ == "total"
        return ndpp
    elseif typ == "mean"
        return ndpp./is
    end
end


function percentage_lifetime_in_phase(seq, phase)
    plip = zeros(8);
    is = zeros(8);
    for i in 1:size(seq)[1]
        # steps in phase 1
        
        p5 = findall(x->((x<=π/4.0)&(x>0)), phase[Int(seq[i,1]):Int(seq[i,2])]);
        p6 = findall(x->((x<=π/2.0)&(x>π/4.0)), phase[Int(seq[i,1]):Int(seq[i,2])]);
        p7 = findall(x->((x<=3*π/4.0)&(x>π/2.0)), phase[Int(seq[i,1]):Int(seq[i,2])]);
        p8 = findall(x->((x<=π)&(x>3*π/4.0)), phase[Int(seq[i,1]):Int(seq[i,2])]);
        p1 = findall(x->((x<=-3*π/4.0)&(x>-π)), phase[Int(seq[i,1]):Int(seq[i,2])]);
        p2 = findall(x->((x<=-π/2.0)&(x>-3*π/4.0)), phase[Int(seq[i,1]):Int(seq[i,2])]);
        p3 = findall(x->((x<=-π/4.0)&(x>-π/2.0)), phase[Int(seq[i,1]):Int(seq[i,2])]);
        p4 = findall(x->((x<=0)&(x>-π/4.0)), phase[Int(seq[i,1]):Int(seq[i,2])]);
        
        # compute duration of event
        duration = length(p1)+length(p2)+length(p3)+length(p4)+length(p5)+length(p6)+length(p7)+length(p8);
        for (i, pp) in enumerate([p1,p2,p3,p4,p5,p6,p7,p8])
            if length(pp) !=0
                is[i] += 1;
                plip[i]+=length(pp)/duration*100;
            end
        end
       
    end
    return plip./is    
    
end




function lifetime_distribution_per_phase(seq, phase)
    d1 = [];
    d2 = [];
    d3 = [];
    d4 = [];
    d5 = [];
    d6 = [];
    d7 = [];
    d8 = [];
    
    for i in 1:size(seq)[1]
        # steps in phase 1
        
        p5 = findall(x->((x<=π/4.0)&(x>0)), phase[Int(seq[i,1]):Int(seq[i,2])]);
        p6 = findall(x->((x<=π/2.0)&(x>π/4.0)), phase[Int(seq[i,1]):Int(seq[i,2])]);
        p7 = findall(x->((x<=3*π/4.0)&(x>π/2.0)), phase[Int(seq[i,1]):Int(seq[i,2])]);
        p8 = findall(x->((x<=π)&(x>3*π/4.0)), phase[Int(seq[i,1]):Int(seq[i,2])]);
        p1 = findall(x->((x<=-3*π/4.0)&(x>-π)), phase[Int(seq[i,1]):Int(seq[i,2])]);
        p2 = findall(x->((x<=-π/2.0)&(x>-3*π/4.0)), phase[Int(seq[i,1]):Int(seq[i,2])]);
        p3 = findall(x->((x<=-π/4.0)&(x>-π/2.0)), phase[Int(seq[i,1]):Int(seq[i,2])]);
        p4 = findall(x->((x<=0)&(x>-π/4.0)), phase[Int(seq[i,1]):Int(seq[i,2])]);
        
        if length(p1) !=0
            push!(d1, length(p1))
        end
        if length(p2) !=0
            push!(d2, length(p2))
        end
        if length(p3) !=0
            push!(d3, length(p3))
        end
        if length(p4) !=0
            push!(d4, length(p4))
        end
        if length(p5) !=0
            push!(d5, length(p5))
        end
        if length(p6) !=0
            push!(d6, length(p6))
        end
        if length(p7) !=0
            push!(d7, length(p7))
        end
        if length(p8) !=0
            push!(d8, length(p8))
        end
    end
    
    return d1, d2, d3, d4, d5, d6, d7, d8    
    
end


function starting_phase(seq, phase)
    # number of events starting in each phase
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
    # number of events stoping in each phase
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

function sep(seq, phase)
     # number of events stoping in each phase
    num_events = size(seq)[1];
    sp = zeros(num_events);
    ep = zeros(num_events);
    
    for i in 1:num_events
        d1 = Int(seq[i,1]);
        d2 = Int(seq[i,2]);
        p1 = phase[d1];
        p2 = phase[d2];
        
        if (p1<=π/4.0)&(p1>0) # phase 5 
            sp[i] = 5; 
        elseif (p1<=π/2.0)&(p1>π/4.0) # phase 6
            sp[i] = 6;
        elseif (p1<=3*π/4.0)&(p1>π/2.0) # phase 7
            sp[i] = 7; 
        elseif (p1<=π)&(p1>3*π/4.0) # phase 8
            sp[i] = 8; 
        elseif (p1<=-3*π/4.0)&(p1>-π) # phase 1
            sp[i] = 1; 
        elseif (p1<=-π/2.0)&(p1>-3*π/4.0) # phase 2 
            sp[i] = 2; 
        elseif (p1<=-π/4.0)&(p1>-π/2.0) # phase 3
            sp[i] = 3; 
        elseif (p1<=0)&(p1>-π/4.0) # phase 4 
            sp[i] = 4; 
        end
        
        if (p2<=π/4.0)&(p2>0) # phase 5 
            ep[i] = 5; 
        elseif (p2<=π/2.0)&(p2>π/4.0) # phase 6
            ep[i] = 6;
        elseif (p2<=3*π/4.0)&(p2>π/2.0) # phase 7
            ep[i] = 7; 
        elseif (p2<=π)&(p2>3*π/4.0) # phase 8
            ep[i] = 8; 
        elseif (p2<=-3*π/4.0)&(p2>-π) # phase 1
            ep[i] = 1; 
        elseif (p2<=-π/2.0)&(p2>-3*π/4.0) # phase 2 
            ep[i] = 2; 
        elseif (p2<=-π/4.0)&(p2>-π/2.0) # phase 3
            ep[i] = 3; 
        elseif (p2<=0)&(p2>-π/4.0) # phase 4 
            ep[i] = 4; 
        end
        
    end
    return collect(zip(sp,ep))
end


function track_minimum(field, maxk, phase)
    Nx = size(field)[1];
    nts = size(field)[end]; 
    
    
    # filter spaially the subset
    field_spatially_filtered = zeros(size(field))
    for ts in 1:nts # filter each time step separately
        fH = fftshift(fft(field[:,ts]));   
        freqs = fftshift(fftfreq(Nx, Nx));
        idx1 = findall(x-> abs(x) >maxk , freqs);
        fH[idx1].=0;
        
        field_spatially_filtered[:,ts] = real(ifft(ifftshift(fH)));
    end
    
    minima = zeros(nts);
    
    # searching boxes
    Box_limits = Dict("p1"=>(7,26), "p2"=>(7,26), "p3"=>(7,26), "p4"=>(7,38), "p5"=>(14,43), "p6"=>(14,43), "p7"=>(14,43), "p8"=>(14,43));
    
    for i in 1:nts
        p = phase[i];
        if (p<=π/4.0)&(p>0) # phase 5 [search in 70ºE - 120ºW <--> spatial indices 14-43]
            j1 = Box_limits["p5"][1];
            j2 = Box_limits["p5"][2];
            minima[i] = j1 - 1 + argmin(field_spatially_filtered[j1:j2,i]);
        elseif (p<=π/2.0)&(p>π/4.0) # phase 6 [search in 70ºE - 120ºW <--> spatial indices 14-43]
            j1 = Box_limits["p6"][1];
            j2 = Box_limits["p6"][2];
            minima[i] = j1 - 1 + argmin(field_spatially_filtered[j1:j2,i]);
        elseif (p<=3*π/4.0)&(p>π/2.0) # phase 7 [search in 70ºE - 120ºW <--> spatial indices 14-43]
            j1 = Box_limits["p7"][1];
            j2 = Box_limits["p7"][2];
            minima[i] = j1 - 1 + argmin(field_spatially_filtered[j1:j2,i]);
        elseif (p<=π)&(p>3*π/4.0) # phase 8 [search in 70ºE - 120ºW <--> spatial indices 14-43]
            j1 = Box_limits["p8"][1];
            j2 = Box_limits["p8"][2];
            minima[i] = j1 - 1 + argmin(field_spatially_filtered[j1:j2,i]);
        elseif (p<=-3*π/4.0)&(p>-π) # phase 1 [search in 30ºE - 140ºE <--> spatial indices 7-26]
            j1 = Box_limits["p1"][1];
            j2 = Box_limits["p1"][2];
            minima[i] = j1 - 1 + argmin(field_spatially_filtered[j1:j2,i]);
        elseif (p<=-π/2.0)&(p>-3*π/4.0) # phase 2 [search in 30ºE - 140ºE <--> spatial indices 7-26]
            j1 = Box_limits["p2"][1];
            j2 = Box_limits["p2"][2];
            minima[i] = j1 - 1 + argmin(field_spatially_filtered[j1:j2,i]);
        elseif (p<=-π/4.0)&(p>-π/2.0) # phase 3 [search in 30ºE - 140ºE <--> spatial indices 7-26]
            j1 = Box_limits["p3"][1];
            j2 = Box_limits["p3"][2];
            minima[i] = j1 - 1 + argmin(field_spatially_filtered[j1:j2,i]);
        elseif (p<=0)&(p>-π/4.0) # phase 4 [search in 30ºE - 150ºW <--> spatial indices 7-38]
            j1 = Box_limits["p4"][1];
            j2 = Box_limits["p4"][2];
            minima[i] = j1 - 1 + argmin(field_spatially_filtered[j1:j2,i]);
        end
        
    end
    
    return field_spatially_filtered, minima
end





function groups_primary_continuing(phase_values)
    groups = [];
    nts = length(phase_values);
        
    # by definition primary events end in phase 8 (on the last date belonging to phase 8), then secondary events start etc
    p8 = findall(x->((x<=π)&(x>3*π/4.0)), phase_values);
    
    if length(p8) != 0 # if the event goes through phase 8 at least once
        cg = consecutive_groups(p8);
        push!(groups, 1:cg[1][end]) # day 1 to last date of the first time we go through phase 8       
        if length(cg)>1 # if the event goes through phase 8 more than once
            for i in 2:length(cg)
                push!(groups, cg[i-1][end]+1:cg[i][end])
            end
        end
        if cg[end][end]!=nts # if the last day is not in phase 8
            push!(groups, cg[end][end]+1:nts)
        end
    else # all phase values are part of one primary event
        push!(groups, 1:nts)
    end
        
    return groups
end

function remove18(days, phase)
    p1 = findall(x->((x<=-3*π/4.0)&(x>-π)), phase);
    p8 = findall(x->((x<=π)&(x>3*π/4.0)), phase);
    to_remove = vcat(p1,p8);
    
    if length(to_remove)!=0
        d1 = collect(days)
        deleteat!(d1, to_remove);
        if length(d1)!=0
            return d1[1]:d1[end]
        else 
            return[]
        end
    else
        return days
    end
end

            

function speed(seq, field, maxk, phase, Nx)
    speeds = []; # store all speed values
    sub_seqs = []; # store days that have been used to compute speeds (used for checking)
    main_seq = [];
    
    lon = range(0, stop=360-360/Nx, length=Nx);
    num_events = size(seq)[1];
    
    for s in 1:num_events # for each event of the sequence
        d1 = Int(seq[s,1]); # first day of event
        d2 = Int(seq[s,2]); # last day
        field_spatially_filtered, minima = track_minimum(field[:,d1:d2], maxk, phase[d1:d2]);
        
        # group dates in primary and continuing events
        g = groups_primary_continuing(phase[d1:d2]); # separate dates belonging to primary, secondary events, etc.
        for sub_g in g
            sub_g = remove18(sub_g, phase[d1:d2][sub_g]);
            nd = length(sub_g);
            if nd >7 # if the event is lasting more than 7 days
                loc_array = lon[Int.(minima[sub_g])]; # location array
                dd = 1:nd; # dates array
                MM = hcat(ones(nd,1), dd);
                
                if all(x -> x==loc_array[1], loc_array) # if the event is completely stationary do not count it
                    push!(speeds,0);
                    push!(sub_seqs, (d1:d2)[sub_g]);
                    push!(main_seq, s);
                #elseif any(x -> x==loc_array[1], loc_array) # then test if at least one point is different, means the wave is only partly stationary 
                #    to_keep = unique(j -> loc_array[j], eachindex(loc_array)) # compute speed with the "moving" points, using linear regression
                #    if length(to_keep)>3 # at least 3 moving points
                #        a, b = (MM'*MM) \(MM'*dd); # linear regression
                #        push!(speeds,(1/b)*111*10^3/(60*60*24));
                #        push!(sub_seqs, (d1:d2)[sub_g]);
                #        push!(main_seq, s);
                #        #a, b = (MM[to_keep,:]'*MM[to_keep,:]) \(MM[to_keep,:]'*dd[to_keep]);
                #        #push!(speeds,(1/b)*111*10^3/(60*60*24)); #[m/s]
                #        #push!(sub_seqs, (d1:d2)[sub_g[to_keep]]);
                #    end
                else # the wave is completely non-stationary
                    a, b = (MM'*MM) \(MM'*loc_array); # linear regression
                    if b>0
                        push!(speeds,(b)*111*10^3/(60*60*24));
                        push!(sub_seqs, (d1:d2)[sub_g]);
                        push!(main_seq, s);
                    end
                end
            end
        end
    end
    return speeds, sub_seqs, main_seq
end



function western_most_point(seq, field, maxk, phase, Nx)
    wmp = []; # store all speed values
    
    lon = range(0, stop=360-360/Nx, length=Nx);
    num_events = size(seq)[1];
    
    for s in 1:num_events # for each event of the sequence
        d1 = Int(seq[s,1]); # first day of event
        d2 = Int(seq[s,2]); # last day
        field_spatially_filtered, minima = track_minimum(field[:,d1:d2], maxk, phase[d1:d2]);
        
        # group dates in primary and continuing events
        g = groups_primary_continuing(phase[d1:d2]); # separate dates belonging to primary, secondary events, etc.
        if length(g)>0
            sub_g = g[1];
            loc_array = lon[Int.(minima[sub_g])]; # location array
            push!(wmp, loc_array[1]);
        end
    end
    return wmp
end


function eastern_most_point(seq, field, maxk, phase, Nx)
    emp = []; # store all speed values
    
    lon = range(0, stop=360-360/Nx, length=Nx);
    num_events = size(seq)[1];
    
    for s in 1:num_events # for each event of the sequence
        d1 = Int(seq[s,1]); # first day of event
        d2 = Int(seq[s,2]); # last day
        field_spatially_filtered, minima = track_minimum(field[:,d1:d2], maxk, phase[d1:d2]);
        
        # group dates in primary and continuing events
        g = groups_primary_continuing(phase[d1:d2]); # separate dates belonging to primary, secondary events, etc.
        if length(g)>0
            sub_g = g[end];
            loc_array = lon[Int.(minima[sub_g])]; # location array
            push!(emp, loc_array[end]);
        end
    end
    return emp
end




function global_minimum(seq, field, maxk, phase, Nx)
    num_events = size(seq)[1];
    location = [];
    value = [];
    
    lon = range(0, stop=360-360/Nx, length=Nx);
    
    for s in 1:num_events # for each event of the sequence
        d1 = Int(seq[s,1]); # first day of event
        d2 = Int(seq[s,2]); # last day
        
        field_spatially_filtered, minima = track_minimum(field[:,d1:d2], maxk, phase[d1:d2]);
        
        # group dates in primary and continuing events
        g = groups_primary_continuing(phase[d1:d2]); # separate dates belonging to primary, secondary events, etc.
        for sub_g in g
            sub_g = remove18(sub_g, phase[d1:d2][sub_g]);
            nd = length(sub_g);
            
            if nd > 7
                sub_field = field_spatially_filtered[:,sub_g];
                min_values = sub_field[CartesianIndex.(zip(Int.(minima[sub_g]),1:64))];
                
                time_of_minimum = argmin(min_values);
                push!(location, lon[Int(minima[sub_g][time_of_minimum])]);
                push!(value, sub_field[Int(minima[sub_g][time_of_minimum]), time_of_minimum]); 
            end
        end
    end
    return location, value;
end


function speed_per_phase(seq, field, maxk, phase, Nx)
    
    speeds = zeros(8); # store the mean speed in each phase
    count = zeros(8); # number of occurence in each phase 
    
    lon = range(0, stop=360-360/Nx, length=Nx);
    
    num_events = size(seq)[1];
    for s in 1:num_events # for each event of the sequence
        d1 = Int(seq[s,1]); # first day of event
        d2 = Int(seq[s,2]); # last day
        field_spatially_filtered, minima = track_minimum(field[:,d1:d2], maxk, phase[d1:d2]);
        
        # dates belonging to each phase 
        for (i, phase_name) = enumerate(["p1", "p2", "p3", "p4", "p5", "p6", "p7", "p8"])
            g = groups(phase_name, phase[d1:d2])
            
            if length(g)!=0
                for sub_g in g
                    nd = length(sub_g);
                    if nd >1
                        loc_array = lon[Int.(minima[sub_g])]; # location array
                        MM = hcat(ones(nd,1), loc_array);
                        dd = 1:nd; # dates array
                        if all(x -> x==loc_array[1], loc_array) # if all locations are equal, wave is stationary (speed = 0 m/s)
                            #speeds[i]+=0; #[m/s]
                            #count[i]+=1;
                            #println("event number ",s, " , phase ", phase_name)    
                            #print(loc_array)
                        elseif any(x -> x!=loc_array[1], loc_array) # then test if at least one point is different, means the wave is only partly stationary 
                            to_keep = unique(j -> loc_array[j], eachindex(loc_array)) # compute speed with the "moving" points, using linear regression
                            if length(to_keep)>2
                                a, b = (MM[to_keep,:]'*MM[to_keep,:]) \(MM[to_keep,:]'*dd[to_keep]);
                                speeds[i]+=(1/b)*111*10^3/(60*60*24); #[m/s]
                                count[i]+=1;
                            end
                        else # the wave is propagating
                            a, b = (MM'*MM) \(MM'*dd); # linear regression
                            speeds[i]+=(1/b)*111*10^3/(60*60*24); #[m/s]
                            count[i]+=1;
                        end
                    end
                end
            end   
        end
    end
    
    return speeds./count;
end




end