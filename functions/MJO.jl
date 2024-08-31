function searchclosest(x,v)
    # x must be a sorted vector!!!
    index = [];
    from = 1;
    to = length(x);
    
    #some dummy values
    index = -10000;
    value = -10000;
    
    # start with a binary search for exact correspondance
    while from <= to 
        mid = floor(Int, ((from + to))/2)
        diff = x[mid]-v;
        if diff == 0 
            index = mid;
            value = x[mid];
            break
        elseif diff < 0 
            from = mid +1;
        elseif diff >0
            to = mid-1
        end
        
    end
    # look for closest value
    if index == value == -10000
        y = x[to:from]
        min_i = argmin(abs.(y.-v))
        value = y[min_i];
        index=to+min_i-1;
    end
        
    return index, value
end

function filtering(X, Nx, nt, L, T, kmin, kmax, wmin, wmax)
    nt = floor(Int,nt)
    κ = (-Nx/2:Nx/2-1);
    ω = (-nt/2:nt/2-1)/T;
    # indices of kmin, kmax, wmin and wmax
    kmini, kmin_val = searchclosest(κ, kmin);
    kmaxi, kmax_val = searchclosest(κ, kmax);
    wmini, wmin_val = searchclosest(ω, wmin);
    wmaxi, wmax_val = searchclosest(ω, wmax);
    kmini1, kmin1_val = searchclosest(κ, -kmin);
    kmaxi1, kmax1_val = searchclosest(κ, -kmax);
    wmini1, wmin1_val = searchclosest(ω, -wmin);
    wmaxi1, wmax1_val = searchclosest(ω, -wmax);
    if wmini1>=wmaxi1
        term=wmaxi1; wmaxi1=wmini1; wmini1=term; 
    end
    if kmini1>=kmaxi1
        term=kmaxi1; kmaxi1=kmini1; kmini1=term; 
    end
    
    
    # spatial transform
    X0 = Complex.(X);
    for kts=1:nt
        X0[:,kts]=fftshift(fft(X0[:,kts]))
    end
    Xfs=X0*0; 
    Xfs[kmini:kmaxi,:]=X0[kmini:kmaxi,:]; 
    for kts=1:nt
        Xfs[:,kts]=ifft(ifftshift(Xfs[:,kts])) 
    end
    ####
    # temporal transform
    for i=1:Nx
        Xfs[i,:]=fftshift(fft(Xfs[i,:]))
    end 
    Xfst=Xfs*0; Xfst[:,wmini:wmaxi]=Xfs[:,wmini:wmaxi]; 
    for i=1:Nx
        Xfst[i,:]=ifft(ifftshift(Xfst[i,:])); 
    end
    ###
    
    # for negative values
    X1=Complex.(X);
    # spatial
    for kts=1:nt
        X1[:,kts]=fftshift(fft(X1[:,kts]));
    end
    X1fs=X1*0; X1fs[kmini1:kmaxi1,:]=X1[kmini1:kmaxi1,:]; 
    for kts=1:nt
        X1fs[:,kts]=ifft(ifftshift(X1fs[:,kts]))
    end
    ###
    # temporal
    for i=1:Nx
        X1fs[i,:]=fftshift(fft(X1fs[i,:])); 
    end
    X1fst=X1*0; X1fst[:,wmini1:wmaxi1]=X1fs[:,wmini1:wmaxi1];
    for i=1:Nx
        X1fst[i,:]=ifft(ifftshift(X1fst[i,:])); 
    end
    ###
    
    Xfst=Xfst+X1fst;
    
end