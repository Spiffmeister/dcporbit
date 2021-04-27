#=
    Performing convergence study
=#

# Is this a guiding center intialisation?
guidingcenter = false

gc₀ = [0.,0.,0.]
v₀ = [1.,0.,0.]

t₀ = 0.
t_f = 1.
Δt = 1.e-1
tol = 1.e-14

# function Bfield(x)
#     α = pi/6
#     if x[3] > 0
#         B = [cos(α),-sin(α),0.]
#     else
#         B = [cos(α),sin(α),0.]
#     end
#     return B
# end

α = pi/6
B₁(x) = [cos(α),-sin(α),0.]
B₂(x) = [cos(α),sin(α),0.]

Bfield = [B₁,B₂]

function event(xv,t)
    chk = xv[3,2]
    return chk
end


p = particle(x₀,v₀,t₀,Δt,MagneticForce,Bfield,Bfield(x₀),1)


integrate(p,t_f,eventfn=eventfn)



