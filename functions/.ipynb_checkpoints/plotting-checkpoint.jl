function circleShape(h,k,r)
    θ = LinRange(0, 2*π, 500);
    h .+ r*sin.(θ), k .+ r*cos.(θ)
end