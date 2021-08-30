#=
    Performing convergence study
=#
using Pkg
Pkg.activate(".")
using orbit


#=
Guiding center and full orbits moving through Z=0 plane
=#
guidingcenter = :fullorbit
gc₀ = [0., 0., 0.]
v₀ = [1.,0.,0.]
Δt = 1.e-3
t_f = 30.
function event(xv)
    chk = sign(xv[3,1])*sign(xv[3,2])
    return chk
end

x₀ = [0.,0.,-0.5]

#=
Run simulation and solve analytically
=#
α = π/6
B₁(x) = [cos(α),sin(α),0.]
B₂(x) = 2.0*[cos(α),-sin(α),0.]
Bfield = [B₁,B₂]

ODE1 = forces(Bfield,event=event)

p1 = particle(x₀,v₀,guidingcenter,Δt,Bfield,1,gc_initial=false)
solve_orbit!(p1,ODE1,t_f)
pe1 = analytic_solve(p1,Bfield,crossing=true,eventfn=event)


#=
Other field
=#
α = π/6
B₁(x) = [cos(α),sin(α),0.]
B₂(x) = [cos(α),-sin(α),0.]
Bfield = [B₁,B₂]

ODE2 = forces(Bfield,event=event)

p2 = particle(x₀,v₀,guidingcenter,Δt,Bfield,1,gc_initial=false)
solve_orbit!(p2,ODE2,t_f)
pe2 = analytic_solve(p2,Bfield,crossing=true,eventfn=event)




#=
    Compute the GC positions to add to plot
=#
ODE_GC₁ = forces(B₁,force=MagneticForce_GC)
gcsim₁ = particle(gc₀,v₀,:guidingcentre,Δt,Bfield[1],1)
ODE_GC₂ = forces(B₂,force=MagneticForce_GC)
gcsim₂ = particle(gc₀,v₀,:guidingcentre,Δt,Bfield[2],1)

solve_orbit!(gcsim₁,ODE_GC₁,t_f)
solve_orbit!(gcsim₂,ODE_GC₂,t_f)



using Plots
using Printf
pyplot()

function plt_gcprojection(f,gc₁,gc₂)
    # GC PROJECTION
    cp = palette(:tab20)

    plt = plot(xlabel="x",ylabel="y",legend=true)
    for k = 1:f.nparts
        gcp = zeros(2,length(f.sp[k].lvol))
        gcp .= NaN
        gcm = zeros(2,length(f.sp[k].lvol))
        gcm .= NaN

        z_p = findall(x->x == 1,f.sp[k].lvol)
        z_m = findall(x->x == 2,f.sp[k].lvol)

        gcp[:,z_p] = f.sp[k].gc[1:2,z_p]
        gcm[:,z_m] = f.sp[k].gc[1:2,z_m]

        plot!(plt,gcp[1,:],gcp[2,:],color=cp[k],label=string("Δt=",f.sp[k].Δt[1],"z<0"))
        plot!(plt,gcm[1,:],gcm[2,:],color=cp[k],label=string("Δt=",f.sp[k].Δt[1],"z>0"))
    end

    plot!(plt,gc₁.x[1,:],gc₁.x[2,:],linestyle=:dash)
    plot!(plt,gc₂.x[1,:],gc₂.x[2,:],linestyle=:dash)

    return plt
end


function plt_avprojection(fe,gc₁,gc₂)
    # AVERAGE ORBIT
    # cp = palette(:tab20)

    plt = plot(xlabel="x",ylabel="y",legend=true,dpi=600,framestyle=:classic)
    for k = 1:f.nparts
        t = range(fe.sp[k].t_boundary[2],stop=fe.sp[k].t[end],length=100)
        ave = zeros(3,length(t))
        for i = 1:length(t)
            ave[:,i] = fe.sp[k].avetraj[:,2] * t[i]
        end
        # plot!(plt,ave[1,:],ave[2,:],label=string("gc₀=",f.sp[k].gc[3,1]))
        plot!(plt,ave[1,:],ave[2,:],label="gc_z=$(@sprintf("%.1e",f.sp[k].gc[3,1]))")
    end

    plot!(plt,gc₁.x[1,:],gc₁.x[2,:],linestyle=:dash,label="z<0 field")
    plot!(plt,gc₂.x[1,:],gc₂.x[2,:],linestyle=:dash,label="z>0 field")

    return plt
end


function plt_td(p1,p2)
    plt = plot3d(p1.x[1,:],p1.x[2,:],p1.x[3,:])
    plt = plot3d!(p2.x[1,:],p2.x[2,:],p2.x[3,:])
    return plt
end


# Plot in 3D
pav = plt_td(p1,p2)


# Plot gcs
pgc = plot(xlabel="x",ylabel="y",legend=true)
plot!(pgc,p1.gc[1,:],p1.gc[2,:])
plot!(pgc,p2.gc[1,:],p2.gc[2,:])

