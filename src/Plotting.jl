


#=
    INTERESTING PLOTS
=#





"""
    plt_orbit(p::Union{particle,exact_particle};lab="Computed")
"""
function plt_orbit(p::Union{particle,exact_particle};lab="Computed")
    # Plot the solution in 3D
    plt = plot3d(p.x[1,:],p.x[2,:],p.x[3,:],label=lab,guidefontsize=16)
    xlabel!("x"); ylabel!("y");
    return plt
end

"""
    plt_gyroradius(v::Array,B::Function,t::Array)
"""
function plt_gyroradius(v::Array,B::Function,t::Array)
    # Gyoradius
    n = size(v)[2]
    gr = zeros(n)
    for i = 1:n
        gr[i] = m*norm(cross(v[:,i],B[:,i]),2)/(q*norm(B[:,i],2))
    end
    plt = plot(t,gr,guidefontsize=16,framestyle=:classic)
    xlabel!("x")
    xlabel!("y")
    return plt
end



#=

#=
    CONVERGENCE TESTS
=#

function plt_error(f::sim,g::analytic_sim)
    n = f.nparts

end



function plt_error(p::Dict,err;ylab="xerr",type="t",tdc=[],savfig=false,n=[0],format="png",modname="")
    if n[1] == 0
        n = collect(range(1,length(p[err]),step=1))
    end
    
    if in(err,["mexerr","meverr"])
        # weird case
        plt = plot(dpi=600,legend=true,guidefontsize=16,framestyle=:classic)
        plot!(p["fo"][n[1]].t,log10.(p[err][n[1]]))
        xlabel!(L"t (s)")
    elseif type=="t"
        # If plotting vs time
        plt = plot(dpi=600,legend=true,guidefontsize=16,framestyle=:classic)
        for i = n
            plot!(p["fo"][i].t,(p[err][i]),label=string(p["Δt"][i][1]))
        end
        if !isempty(tdc)
            vline!(tdc,label="Discontinuity")
        end
        xlabel!("t (s)")
    elseif type=="dt"
        # If plotting vs dt
        plt = plot(dpi=600,legend=true,guidefontsize=16,framestyle=:classic)
        Δx = zeros(n[end])
        Δt = zeros(n[end])
        for i = n
            Δx[i] = maximum(p[err][i])
            # Δx[i] = log10(Δx[i])
            Δt[i] = p["Δt"][i][1]
        end
        plot!(log10.(Δt),log10.(Δx),marker=true,label="Computed error")
        xlabel!(L"log_{10}(\Delta t)")

        # plot 4th order error
        p_er = Δt.^4

        @. model(x,p) = p[1]*x + p[2]
        fit = curve_fit(model,log10.(Δt[2:end-1]),log10.(Δx[2:end-1]),[4.,0.])
        plot!(log10.(Δt),log10.(Δt)*fit.param[1].+fit.param[2],color="black",linestyle=:dash,label=L"\mathcal{O}([\Delta t]^4)")

        # plot!(log10.(Δt),log10.(p_er),color="black",linestyle=:dash,label=L"\mathcal{O}([\Delta t]^4)")

    end

    ylabel!(ylab)

    if savfig
        savefig(plt,string("../figures/",err,"_",type,"_error",modname,".",format))
    end
    return plt
end

# Plot the energy conservation
function plt_energycons(p::Dict;type="t",savfig=false,n=[0],format="pdf",tdc=[])
    # Plot lines of energy consveration vs time
    if n[1] == 0
        n = collect(range(1,length(p["E"]),step=1))
    end
    # plt = plot(dpi=600,legend=true,yformatter=x->string(Int(x/1.e10),"*pow10"))
    
    if type == "t"
        # plot vs t
        plt = plot(dpi=600, legend=true,guidefontsize=16,framestyle=:classic)
        for i = n
            ΔE = abs.(p["E"][i]/p["E"][i][1] .- 1)
            ΔE = log10.(ΔE)
            ΔE[findall(x->x==-Inf,ΔE)] .= -16.
            plot!(p["fo"][i].t,ΔE,label=string(p["fo"][i].Δt))
        end
        xlabel!("t (s)")
        ylabel!(L"log_{10}(E/E_0 - 1)")
        
        if (!isempty(tdc))
            vline!(tdc,label="Crossing")
        end
        
    elseif type == "dt"
        # plot vs dt
        plt = plot(dpi=600, legend=true,guidefontsize=16,framestyle=:classic)
        ΔE=zeros(n[end])
        Δt=zeros(n[end])
        for i = n
            ΔE[i] = maximum(abs.(p["E"][i]/p["E"][i][1] .- 1))
            ΔE[i] = log10(ΔE[i])
            Δt[i] = p["Δt"][i][1]
        end
        plot!(log10.(Δt),ΔE,marker=true,label="Energy loss")
        xlabel!(L"log_{10}(\Delta t)")
        ylabel!(L"log_{10} (max(E/E_0 - 1))")

        # plot 4th order error
        p_er = Δt.^4
        # plot!(log10.(Δt),log10.(p_er),marker=true,line=(:dash),color="black",label=L"\mathcal{O}([\Delta t]^4)")
        
        @. model(x,p) = p[1]*x + p[2]
        fit = curve_fit(model,log10.(Δt[1:end-3]),ΔE[1:end-3],[4.,0.])
        plot!(log10.(Δt),log10.(Δt)*fit.param[1].+fit.param[2],color="black",linestyle=:dash,label=L"\mathcal{O}([\Delta t]^5)")#label=string(fit.param[1]))

    end

    if savfig
        savefig(plt,string("../figures/dE","_",type,".",format))
    end
    return plt
end

# Plot conservation of other properties i.e. velocity
function plt_cons(p::Dict,val,tdc;ylab="",savfig=false,n=[0],format="png")
    if n[1] == 0
        n = collect(range(1,length(p["fo"]),step=1))
    end
    plt = plot(dpi=600,legend=true,guidefontsize=16,framestyle=:classic)

    for i = n
        x = getfield(p["fo"][i],val)
        Δx = zeros(size(x)[2])
        for j = 1:size(x)[2]
            Δx[j] = abs(norm(x[:,j],2) .- norm(x[:,1],2))
        end
        
        plot!(p["fo"][i].t,log10.(Δx),label=string(p["fo"][i].Δt))
    end

    vline!(tdc,label="Discontinuity")
    xlabel!("t (s)")

    if ylab != ""
        ylabel!(ylab)
    end
    if savfig
        savefig(plt,string("../figures/",val,"_cons.",format))
    end
    return plt
end


function plt_moderr(err,moderr,tdc,t;lab="xerr",savfig=false,format="png")
    # Plot error of exact soln with computational boundary
    plt = plot(dpi=600,legend=true,guidefontsize=16,framestyle=:classic)
    plot!(plt,t,log10.(err),label="Exact interface")
    plot!(plt,t,log10.(moderr),label="Computation interface")
    vline!(plt,tdc,label="Discontinuity")

    ylabel!(L"log_{10}(x_{e} - x)")
    xlabel!("t (s)")

    if savfig
        savefig(plt,string("../figures/moderr_",lab,".",format))
    end
    return plt
end

function plt_momcons(p,xy;n=[],tdc=[],lab="p_x",savfig=false,format="png")
    # Plot momentum conservation
    if isempty(n)
        n = collect(range(1,length(p["fo"]),step=1))
    end
    plt = plot(dpi=600,legend=:topleft,guidefontsize=16,framestyle=:classic)

        # HARD CODE WARNING
        # leg = ["1/10","1/80","1/1280"]

    for i = n
        plot!(plt,p["fo"][i].t,log10.(abs.(p["p"][i][xy,:]/p["p"][i][xy,1] .- 1)),label=string(L"\Delta t = 1/",10*2^(i-1)))
    end
    if !isempty(tdc)
        vline!(tdc,label="Discontinuity")
    end

    xlabel!("t (s)")
    if xy == 1
        ylabel!(L"log_{10}(\Delta p_x/p_0 - 1)")
        fname = string("../figures/momentumcons_x.",format)
    else
        ylabel!(L"log_{10}(\Delta p_y/p_0 - 1)")
        fname = string("../figures/momentumcons_y.",format)
    end

    if savfig
        savefig(plt,fname)
    end
    return plt
end


=#