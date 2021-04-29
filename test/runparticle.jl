#=
    Performing convergence study
=#

# Is this a guiding center intialisation?
guidingcentre = false

x₀ = [0.,0.,0.]
v₀ = [1.,0.,0.]

t₀ = 0.
t_f = 1.
Δt = 1.e-1
tol = 1.e-14


Bfield(x) = [0.,0.,1.]

p = particle(x₀,v₀,guidingcentre,t₀,Δt,MagneticForce,Bfield,1)

integrate!(p,t_f)



exactsolve(x₀,v₀,p.t,Bfield,crossing=false)