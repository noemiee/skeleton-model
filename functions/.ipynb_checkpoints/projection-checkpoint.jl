function compute_projection(KK, RR, AA, ZZ, Nx, nt, L, T, dofilter, kmin=1, kmax=3, wmin=1/90, wmax=1/30)
    nt = floor(Int, nt); 
    # X = [K, R, A, Z]
    # fomat x x Nt
    Kf = KK; Rf = RR; Af = AA; Zf = ZZ;
    # first filter all signals if needed 
    if dofilter==1
        Kf = filtering(KK, Nx, nt, L, T, kmin, kmax, wmin, wmax);
        Rf = filtering(RR, Nx, nt, L, T, kmin, kmax, wmin, wmax);
        Af = filtering(AA, Nx, nt, L, T, kmin, kmax, wmin, wmax);
        Zf = filtering(ZZ, Nx, nt, L, T, kmin, kmax, wmin, wmax);
    end
    ####################################
    
    # second compute the spatial (zonal) FFT of all signal as a function of time
    X = Kf;
    for ts in 1:nt
        X[:, ts] = fftshift(fft(X[:,ts]));
    end
    Kf = X;
    X = Rf;
    for ts in 1:nt
        X[:, ts] = fftshift(fft(X[:,ts]));
    end
    Rf = X;
    X = Af;
    for ts in 1:nt
        X[:, ts] = fftshift(fft(X[:,ts]));
    end
    Af = X;
    X = Zf;
    for ts in 1:nt
        X[:, ts] = fftshift(fft(X[:,ts]));
    end
    Zf = X;
    ###################################################
    
    # third compute eigenmodes and projection
    κ = Int.(-Nx/2:Nx/2-1);
    Xproj = Complex.(zeros(Nx, nt));
    for kval in κ
        ωᵣₑ, ωᵢₘₐ, V = compute_eigenmodes(kval, Hbar, Qbar, Γ, γ);
        vec = V[:,2];
        # by convension we choose K with phase 0
        phaseK = angle(vec[1]);
        vec = vec.*exp(-im*phaseK);
        for ts in 1:nt
            Xsim = [Kf[kval+33,ts], Rf[kval+33,ts], Af[kval+33,ts], Zf[kval+33,ts]]; # +33 is just to obtain a positive index because kval are in [-32,32]
            Xproj[kval+33,ts] = dot(Xsim, vec)
        end    
    end
    for ts in 1:nt
        Xproj[:,ts] = ifft(ifftshift(Xproj[:,ts])) #ifft((Xproj[:,ts]))
    end
    return Xproj
end