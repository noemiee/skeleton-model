module ev_functions

export evolve_anom
export evolve_stocha_anom
export evolve
export evolve_stocha

using FFTW
using Statistics
using Random, Distributions
using Measures
using StatsBase




function evolve_anom(K0, R0, Z0, A0, Hbar, Qbar, Γ, Δt, κ, dx, As)
    γ = 1.0/π^(1.0/4.0);
    # first evolve K and R in Fourier space
    
    # Kelvin wave
    K0hat = fft(K0);
    K0hat = fftshift(K0hat); # rearanges the data so that negative frequencies are on the left
    α = κ;
    rhs = -1.0/√2 .* (Hbar .*(A0.-As))*γ/γ^2;
    rhs_hat = fft(rhs);
    rhs_hat = fftshift(rhs_hat);
    
    K1hat = K0hat .* exp.(-α .*Δt) .+ rhs_hat .*(1 .- exp.(-α .*Δt)) ./α;
    index0 = findall(k -> k==0, κ)
    K1hat[index0] = K0hat[index0] + rhs_hat[index0]*Δt;
    
    K1hat = ifftshift(K1hat);
    K1 = real(ifft(K1hat));
    ##OK
 
    # Rossby wave
    R0hat = fft(R0);
    R0hat = fftshift(R0hat); # rearanges the data so that negative frequencies are on the left
    α = -1.0/3.0*κ;
    rhs = -2.0*√2/3.0 .*(Hbar .*(A0.-As))*γ/γ^2;
    rhs_hat = fft(rhs);
    rhs_hat = fftshift(rhs_hat);
    R1hat = R0hat .*exp.(-α .*Δt) .+ rhs_hat .*(1 .- exp.(-α .*Δt)) ./α;
    index0 = findall(k -> k==0, κ)
    R1hat[index0] = R0hat[index0] + rhs_hat[index0]*Δt;
   
    
    R1hat = ifftshift(R1hat);
    R1 = real(ifft(R1hat));    
    ## OK
    
    # use the new K and R (K1 and R1) to evolve Z and a
    Z1 = Z0 .+ Δt .*(Qbar - 1.0)*(Hbar .*(A0.-As));
    
    theta = -1/√2 .*(K1 .+ R1 ./ 2.0)*γ; #OK #theta = -1/√2 .*(K0 .+ R0 ./ 2.0)*γ; #OK
    Q1 = Z1 .- Qbar .*theta; #Q1 = Z0 .- Qbar .*theta;
    
    A1 = A0.*exp.(Γ*Q1*Δt); 
    
    return K1, R1, A1, Z1
    
end


function evolve_stocha_anom(K0, R0, Z0, A0, Hbar, Qbar, Γ, Δt, κ, Δa, Nx, As)
    ## K, R, Z, A are anomalies ##
    γ = 1.0/π^(1.0/4.0);
    # first evolve K and R in Fourier space
    
    # Kelvin wave
    K0hat = fft(K0);
    K0hat = fftshift(K0hat); # rearanges the data so that negative frequencies are on the left
    α = κ;
    rhs = -1.0/√2 .* (Hbar .*A0)*γ/γ^2;
    rhs_hat = fft(rhs);
    rhs_hat = fftshift(rhs_hat);
    
    K1hat = K0hat .* exp.(-α .*Δt) .+ rhs_hat .*(1 .- exp.(-α .*Δt)) ./α;
    index0 = findall(k -> k==0, κ)
    K1hat[index0] = K0hat[index0] + rhs_hat[index0]*Δt;
    
    K1hat = ifftshift(K1hat);
    K1 = real(ifft(K1hat));
    ##OK
 
    # Rossby wave
    R0hat = fft(R0);
    R0hat = fftshift(R0hat); # rearanges the data so that negative frequencies are on the left
    α = -1.0/3.0*κ;
    rhs = -2.0*√2/3.0 .*(Hbar .*A0)*γ/γ^2;
    rhs_hat = fft(rhs);
    rhs_hat = fftshift(rhs_hat);
    R1hat = R0hat .*exp.(-α .*Δt) .+ rhs_hat .*(1 .- exp.(-α .*Δt)) ./α;
    index0 = findall(k -> k==0, κ)
    R1hat[index0] = R0hat[index0] + rhs_hat[index0]*Δt;
   
    
    R1hat = ifftshift(R1hat);
    R1 = real(ifft(R1hat));    
    ## OK
    
    # use the new K and R (K1 and R1) to evolve Z and a
    # stochastic update with inner stochastic time
    
    # stochastic loop initial contitions
    θ = -1/√2 .*(K1 .+ R1 ./ 2.0)*γ;
    Z1 = Z0;
    η = floor.(Int, A0/Δa); 
    
    for point in 1:Nx
        ts = 0; 
        while ts<Δt
            # compute the transition rates λ (up) and μ (down)
            q = Z1[point] .- Qbar .*θ[point];
            if η[point]>0 
                if q>=0
                    λ = Γ*q*η[point]; # up rate
                    μ = 0.0; # down rate
                else
                    λ = 0.0;
                    μ = -Γ*q*η[point];
                end
            else
                λ = 1.0;
                μ = 0.0; 
            end
            # compute time to the next transisiton
            rr = rand();
            while rr==0.0
                rr = rand();
            end
            τ = -1/(λ+μ)*log(rr);
            # if not the last jump, 
            if ts+τ<=Δt
                # update Z and then η according to transition probabilities
                Z1[point] = Z1[point] + τ*(Qbar-1)*Hbar*(η[point]*Δa);
                rrr = rand();
                while rrr == 0.0 
                    rrr=rand();
                end
                if rrr <= λ/(λ+μ)
                    η[point]+=1;
                else
                    η[point]-=1;
                end
                if η[point]<0
                    η[point]=0;
                end
            else
                # only update Z until the end of the time step and do not change η
                Z1[point] = Z1[point] + (Δt-ts)*(Qbar-1)*Hbar*(η[point]*Δa);
            end
            ts = ts+τ;           
        end
    end
    
    q = Z1 .- Qbar.*θ;
    A1 = η.*Δa .+ Γ*q.*As ;  # Exp[η.*Δa] =  Γ*Q_anom*A_anom  
    return K1, R1, A1, Z1
    
end





function evolve(K0, R0, Z0, A0, Hbar, Qbar, Stheta, Sq, Γ, Δt, κ, dx)
    γ = 1.0/π^(1.0/4.0);
    # first evolve K and R in Fourier space
    
    # Kelvin wave
    K0hat = fft(K0);
    K0hat = fftshift(K0hat); # rearanges the data so that negative frequencies are on the left
    α = κ;
    rhs = -1.0/√2 .* (Hbar .*A0 .- Stheta)*γ/γ^2;
    rhs_hat = fft(rhs);
    rhs_hat = fftshift(rhs_hat);
    
    K1hat = K0hat .* exp.(-α .*Δt) .+ rhs_hat .*(1 .- exp.(-α .*Δt)) ./α;
    index0 = findall(k -> k==0, κ)
    K1hat[index0] = K0hat[index0] + rhs_hat[index0]*Δt;
    
    K1hat = ifftshift(K1hat);
    K1 = real(ifft(K1hat));
    ##OK
 
    # Rossby wave
    R0hat = fft(R0);
    R0hat = fftshift(R0hat); # rearanges the data so that negative frequencies are on the left
    α = -1.0/3.0*κ;
    rhs = -2.0*√2/3.0 .*(Hbar .*A0 .- Stheta)*γ/γ^2;
    rhs_hat = fft(rhs);
    rhs_hat = fftshift(rhs_hat);
    R1hat = R0hat .*exp.(-α .*Δt) .+ rhs_hat .*(1 .- exp.(-α .*Δt)) ./α;
    index0 = findall(k -> k==0, κ)
    R1hat[index0] = R0hat[index0] + rhs_hat[index0]*Δt;
   
    
    R1hat = ifftshift(R1hat);
    R1 = real(ifft(R1hat));    
    ## OK
    
    
    # use the new K and R (K1 and R1) to evolve Z and a
    Z1 = Z0 .+ Δt .*((Qbar-1.0)*(Hbar.*A0) .+ (Sq .- Qbar*Stheta) );
    
    theta = -1/√2 .*(K1 .+ R1 ./ 2.0)*γ; #OK #theta = -1/√2 .*(K0 .+ R0 ./ 2.0)*γ; #OK
    Q1 = Z1 .- Qbar .*theta; #Q1 = Z0 .- Qbar .*theta;
        
    A1 = A0.*exp.(Γ*Q1*Δt); 
    
    return K1, R1, A1, Z1
end


function evolve_stocha(K0, R0, Z0, A0, Hbar, Qbar, Stheta, Sq, Γ, Δt, κ, Δa, Nx)
    γ = 1.0/π^(1.0/4.0);
    # first evolve K and R in Fourier space
    
    # Kelvin wave
    K0hat = fft(K0);
    K0hat = fftshift(K0hat); # rearanges the data so that negative frequencies are on the left
    α = κ;
    rhs = -1.0/√2 .* (Hbar .*A0 .- Stheta)/γ;
    rhs_hat = fft(rhs);
    rhs_hat = fftshift(rhs_hat);
    
    K1hat = K0hat .* exp.(-α .*Δt) .+ rhs_hat .*(1 .- exp.(-α .*Δt)) ./α;
    index0 = findall(k -> k==0, κ)
    
    K1hat[index0] = K0hat[index0] + rhs_hat[index0]*Δt;
    
    K1hat = ifftshift(K1hat);
    K1 = real(ifft(K1hat));
    ##OK
 
    # Rossby wave
    R0hat = fft(R0);
    R0hat = fftshift(R0hat); # rearanges the data so that negative frequencies are on the left
    α = -1.0/3.0*κ;
    rhs = -2.0*√2/3.0 .*(Hbar .*A0 .- Stheta)/γ;
    rhs_hat = fft(rhs);
    rhs_hat = fftshift(rhs_hat);
    R1hat = R0hat .*exp.(-α .*Δt) .+ rhs_hat .*(1 .- exp.(-α .*Δt)) ./α;
    index0 = findall(k -> k==0, κ)
    R1hat[index0] = R0hat[index0] + rhs_hat[index0]*Δt;
   
    
    R1hat = ifftshift(R1hat);
    R1 = real(ifft(R1hat));    
    ## OK
    
    # use the new K and R (K1 and R1) to evolve Z and a
    # stochastic update with inner stochastic time
    
    # stochastic loop initial contitions
    θ = -1/√2 .*(K1 .+ R1 ./ 2.0)*γ;
    Z1 = Z0;
    η = floor.(Int, A0/Δa); 
    
    for point in 1:Nx
        ts = 0; 
        while ts<Δt
            # compute the transition rates λ (up) and μ (down)
            q = Z1[point] .- Qbar .*θ[point];
            if η[point]<0
                println("oups")
            end
            if η[point]>0 
                if q>=0
                    λ = Γ*q*η[point]; # up rate
                    μ = 0.0; # down rate
                else
                    λ = 0.0;
                    μ = -Γ*q*η[point];
                end
            else
                λ = 1.0;
                μ = 0.0; 
            end
            # compute time to the next transisiton
            rr = rand();
            while rr==0.0
                rr = rand();
            end
            τ = -1/(λ+μ)*log(rr);
            # if not the last jump, 
            if ts+τ<=Δt
                # update Z and then η according to transition probabilities
                Z1[point] = Z1[point] + τ*(Qbar*(Hbar*η[point]*Δa - Stheta[point]) -(Hbar*η[point]*Δa - Sq[point]));
                #Z1[point] = Z1[point] + τ*(1-Qbar)*(Sq[point] -Hbar*η[point]*Δa);
                rrr = rand();
                while rrr == 0.0 
                    rrr=rand();
                end
                if rrr <= λ/(λ+μ)
                    η[point]+=1;
                else
                    η[point]-=1;
                end
                if η[point]<0
                    η[point]=0;
                end
            else
                # only update Z until the end of the time step and do not change η
                Z1[point] = Z1[point] + (Δt-ts)*(Qbar*(Hbar*η[point]*Δa - Stheta[point]) -(Hbar*η[point]*Δa - Sq[point]));
            end
            ts = ts+τ;           
        end
    end
    
    A1 = η.*Δa;
    return K1, R1, A1, Z1
    
end







function distance(point)
    return (point.x^2+point.y^2)^0.5
end

end