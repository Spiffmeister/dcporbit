




#=
    MAGNETIC FIELDS
=#
function Bfield(x)
    α = pi/6
    if x[3] > 0
        B = [cos(α),-sin(α),0]
    else
        B = [cos(α),sin(α),0]
    end
    return B
end


p = particle(x₀,v₀,t₀,Δt,MagneticForce,Bfield,Bfield(x₀),1)


#=
    EVENT FUNCTION
=#
function event(xv,t)
    chk = xv[3,2]
    return chk
end








