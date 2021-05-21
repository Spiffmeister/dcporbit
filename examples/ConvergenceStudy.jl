#=
    Performing convergence study
=#

# Is this a guiding center intialisation?


tol = 1.e-14


# Simulation parameters
nparts = 8

guidingcenter = false
gc₀ = [0., 0., 0.]
v₀ = [1.,0.,0.]
Δt = [2^-i for i in 1:nparts]
t_f = 1.

α = π/6
B₁(x) = [cos(α),-sin(α),0.]
B₂(x) = [cos(α),sin(α),0.]

Bfield = [B₁,B₂]

function event(xv,t)
    chk = xv[3,2]
    return chk
end



f = sim(nparts,gc₀,v₀,guidingcenter,Δt,dvdt,Bfield,1)

run_sim!(f,dvdt,t_f)



pe = analytic_solve(x₀,v₀,p.t,Bfield,crossing=false)

# pe = analytic_solve(x₀,v₀,p.t,Bfield,crossing=false)


