using Plots: length
#=
    Performing convergence study
=#

# Is this a guiding center intialisation?
using Pkg
Pkg.activate(".")
using orbit

using LinearAlgebra
using LaTeXStrings
using Plots
using LsqFit
pyplot()




tol = 1.e-14


# Simulation parameters
nparts = 6

guidingcenter = :fullorbit
gc₀ = [0., 0., -4.]
v₀ = [1.,0.,0.]
# Δt = [2.0^-i for i in 1:nparts]
Δt = [1.e-1/2^(i-1) for i in 1:nparts]
t_f = 200.

α = π/6
B₁(x) = [cos(α),sin(α),0.]
B₂(x) = [cos(α),-sin(α),0.]

Bfield = [B₁,B₂]

function event(xv)
    chk = sign(xv[3,1])*sign(xv[3,2])
    return chk
end


ODE = forces(Bfield,event=event)

f = sim(nparts,gc₀,v₀,guidingcenter,Δt,Bfield,1)

run_sim!(f,ODE,t_f)

g = analytic_solve(f,Bfield,crossing=true,eventfn=event)











#========== DEF PLOTS ==========#


function diag_momentum(part)
    p = zeros(2,length(part.t))
    p[1,:] = part.v[1,:] - abs.(part.x[3,:]) * sin(π/6)
    p[2,:] = part.v[2,:] - part.x[3,:] * cos(π/6)
    return p
end

function plt_momentum_conservation(sim,n::Vector{Int64})
    # Plot momentum conservation
    p_x = plot(dpi=600,legend=:topleft,guidefontsize=16,framestyle=:classic)
    p_y = plot(dpi=600,legend=:topleft,guidefontsize=16,framestyle=:classic)

    # sim.sp[i].t / 2π

    for i = n
        p = diag_momentum(sim.sp[i])

        plot!(p_x,sim.sp[i].t/(2π),log10.(abs.(p[1,:]/p[1,1] .- 1)),label=string(L"\Delta t = 1/",10*2^(i-1)))
        plot!(p_y,sim.sp[i].t/(2π),log10.(abs.(p[2,:]/p[2,1] .- 1)),label=string(L"\Delta t = 1/",10*2^(i-1)))

    end

    xlabel!(p_x,"t (s)")
    xlabel!(p_y,"t (s)")

    ylabel!(p_x,L"log_{10}(\Delta p_x/p_{x,0} - 1)")
    ylabel!(p_y,L"log_{10}(\Delta p_y/p_{y,0} - 1)")

    savefig(p_x,"Figures/momentumcons_x.pdf")
    savefig(p_y,"Figures/momentumcons_y.pdf")
end




function plt_error(sim,exact)

    x_err = plot(dpi=600,legend=true,guidefontsize=16,framestyle=:classic)
    v_err = plot(dpi=600,legend=true,guidefontsize=16,framestyle=:classic)
    
    Δt = [x.Δt for x in f.sp]
    
    Δx = zeros(length(sim.sp))
    Δv = zeros(length(sim.sp))

    for i = 1:length(sim.sp)
        # Loop over simulations
        for j = 1:length(sim.sp[i].t)
            # Loop over points and take norms
            x_error = norm(sim.sp[i].x[:,j] - exact.sp[i].x[:,j])#/norm(exact.sp[i].x[:,j])
            x_error > Δx[i] ? Δx[i] = x_error : nothing
            
            v_error = norm(sim.sp[i].v[:,j] - exact.sp[i].v[:,j])
            v_error > Δv[i] ? Δv[i] = v_error : nothing
        end

    end

    # 4th order error
    # p_er = Δt.^4

    @. x_model(x,p) = p[1]*x + p[2]
    fit_x = curve_fit(x_model,log10.(Δt[2:end-1]),log10.(Δx[2:end-1]),[4.,0.])
    plot!(x_err,log10.(Δt),log10.(Δt)*4.0 .+ fit_x.param[2],color="black",linestyle=:dash,label=L"\mathcal{O}([\Delta t]^4)")
    
    @. v_model(v,p) = p[1]*v + p[2]
    fit_v = curve_fit(v_model,log10.(Δt[2:end-1]),log10.(Δv[2:end-1]),[4.,0.])
    plot!(v_err,log10.(Δt),log10.(Δt)*4.0 .+ fit_v.param[2],color="black",linestyle=:dash,label=L"\mathcal{O}([\Delta t]^4)")
    
    # Print the params to check
    println("x error order: ",fit_x.param[1])
    println("v error order: ",fit_v.param[1])
    
    # Plot the curves
    plot!(x_err,log10.(Δt),log10.(Δx),marker=true,label="Computed error")
    plot!(v_err,log10.(Δt),log10.(Δv),marker=true,label="Computed error")
    
    # Labels
    xlabel!(x_err,L"log_{10}(\Delta t)")
    xlabel!(v_err,L"log_{10}(\Delta t)")
    
    ylabel!(x_err,L"log_{10}(max(x_{e} - x))")
    ylabel!(v_err,L"log_{10}(max(v_{e} - v))")

    # Save
    savefig(x_err,"Figures/x_error.pdf")
    savefig(v_err,"Figures/v_error.pdf")
end



function plt_full_error(sim,exact,npart)

    x = zeros(Float64,size(sim.sp[npart].x)[2])
    for i = 1:size(sim.sp[npart].x)[2]
        x[i] = norm(sim.sp[npart].x[:,i] .- exact.sp[npart].x[:,i])/norm(exact.sp[npart].x[:,i])
    end

    plta = plot(sim.sp[npart].t,x)
    savefig(plta,"Figures/x_total_error.pdf")

end



#========== PLOTTING ==========#
plt_momentum_conservation(f,[1,4,6])
plt_error(f,g)

plt_full_error(f,g,6)






