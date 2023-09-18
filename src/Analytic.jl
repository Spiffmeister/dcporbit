#===
    Full orbit and gyrocenter equations
===#



"""
    exact_particle{T<:Real}
"""
struct exact_particle{T<:Real}
    # Structure for hold the exact solutions to particles
    x           :: Array{T}
    v           :: Array{T}
    x_boundary  :: Array{T}
    v_boundary  :: Array{T}
    t           :: Vector{T}
    t_boundary  :: Vector{T}
    avetraj     :: Array{T}
end


"""
    analytic_sim
"""
mutable struct analytic_sim
    # Container for multiple analytic particles
    nparts  :: Int64
    sp      :: Array{exact_particle}

    function analytic_sim(nparts::Int64)
        parts = Array{exact_particle}(undef,nparts)
        new(nparts,parts)
    end
end





#===
    METHODS FOR SOLVING PARTICLES
===#
"""
    analytic_solve
"""
function analytic_solve end
"""
    analytic_solve(p::particle,Bfield;crossing=true,eventfn=nothing)
"""
function analytic_solve(p::particle,Bfield;crossing=true,eventfn=nothing)
    # Method for taking in a particle and returning the analytic solve
    x₀ = p.gc[:,1]
    v₀ = p.v[:,1]
    t = p.t
    lvol = p.lvol[1]
    crossing != diff(p.lvol .== 0)
    # Compute the analytic solution based on the simulation
    pe = analytic_solve(x₀,v₀,t,Bfield,crossing=crossing,eventfn=eventfn,lvol=lvol)
    return pe
end
"""
    analytic_solve(simulation::sim,Bfield;crossing=true,eventfn=nothing::Union{Nothing,Function})
"""
function analytic_solve(simulation::sim,Bfield;crossing=true,eventfn=nothing::Union{Nothing,Function})
    # Method for taking in a simulation
    asym = analytic_sim(simulation.nparts)
    # Loop over simulations and compute analytic simulations
    for i = 1:simulation.nparts
        x₀ = simulation.sp[i].gc[:,1]
        v₀ = simulation.sp[i].v[:,1]
        t = simulation.sp[i].t
        lvol = simulation.sp[i].lvol[1]
        # Check if there are crossings
        crossing = any(diff(simulation.sp[i].lvol) .!= 0)
        
        if (crossing) & (typeof(eventfn) <: Nothing)
            error("Cannot have a crossing case with no event function")
        end

        asym.sp[i] = analytic_solve(x₀,v₀,t,Bfield,crossing=crossing,eventfn=eventfn,lvol=lvol)
    end
    return asym
end

#
# ANALYTIC SOLUTIONS FOR ORBITS
#
"""
    analytic_solve(x₀::Vector{Float64},v₀::Vector{Float64},t::Vector{Float64},Bfield::Union{Function,Array{Function}};crossing=true,eventfn=nothing,lvol=1::Int64)
"""
function analytic_solve(x₀::Vector{Float64},v₀::Vector{Float64},t::Vector{Float64},Bfield::Union{Function,Array{Function}};crossing=true,eventfn=nothing,lvol=1::Int64)
    # Works for static fields
    maxvol = 2

    if typeof(Bfield) <: Function
        Bf = Bfield
    elseif typeof(Bfield) <: Array
        Bf = Bfield[lvol]
    end

    ω = [abs(q)/m * norm(Bf(x₀,t[1]))]
    x = x_b = x₀ - guiding_center(v₀,Bf(x₀,t[1]))
    v = v_b = v₀
    
    pB = zeros(Float64,(3,1))
    pB .= Bf(x₀,t[1])
    x̄ = Array{Float64,2}(undef,3,0)
    b₁ = Array{Float64,2}(undef,3,0)
    # Crossing time
    τ_b = Vector{Float64}(undef,0)

    if !crossing
        _analytic_solve_nocrossing(x₀,v₀,t,Bf)
    end
    k = 0
    while crossing
        # If there are crossings solve per region
        B = Bfield[lvol](x_b[:,end],t)
        b = magcoords(v₀,B)
        b₁ = hcat(b₁,b[1]) #Store the main vector so we can compute average orbits at the end
        # ω = abs(q)/m * norm(B)
        append!(ω,abs(q)/m * norm(B))
        τ = bound_time(ω[end],b) #Compute crossing time
        if τ != 0 & (τ < t[end])
            # If crossing time found store it
            append!(τ_b,τ)
            #Get all values of t in this section
            tᵢ = findall(x->sum(τ_b[1:end-1])<=x<=sum(τ_b),t)
            t_f = t[tᵢ[end]] #For exit condition
            # zero the crossing time
            tᵢ = t[tᵢ] .- sum(τ_b[1:end-1])
            # Eqns of motion
            xᵢ = exact_x(v₀,x₀,b,tᵢ,ω[end])
            vᵢ = exact_v(v₀,b,tᵢ,ω[end])
            # Positions at boundary
            x_b = hcat(x_b,exact_x(v₀,x₀,b,τ,ω[end]))
            v_b = hcat(v_b,exact_v(v₀,b,τ,ω[end]))
            # Store positions for returning
            x = hcat(x,xᵢ)
            v = hcat(v,vᵢ)
            # Update the field to use
            lvol = mod1(lvol+1,maxvol)
        else
            # If crossing time cannot be computed then iterate
            if length(τ_b) == 0
                tᵢ = t
            else
                # If happens after a crossing zero the time
                tᵢ = t - τ_b
                tᵢ = tᵢ[tᵢ .<= 0]
            end
            i = 2
            xₜₘₚ = x[:,end]
            vₜₘₚ = v[:,end]
            while (t[i] < t[end])
                # Loop over timesteps until boundary crossing
                xₜₘₚ = hcat(xₜₘₚ,exact_x(v₀,x₀,b,tᵢ[i],ω[end]))
                vₜₘₚ = hcat(vₜₘₚ,exact_v(v₀,b,tᵢ[i],ω[end]))
                # x = hcat(x,exact_x(v₀,x₀,b,tᵢ[i],ω[end]))
                # v = hcat(v,exact_v(v₀,b,tᵢ[i],ω[end]))
                chk = eventfn(xₜₘₚ[:,end-1],xₜₘₚ[:,end])
                if chk < 0.
                    # If crossing detected comput the exact position
                    ex(t) = exact_x(v₀,x₀,b,t,ω[end])[3]
                    τ = find_zero(ex,(tᵢ[i-1],tᵢ[i+1]),Bisection(),atol=1.e-15)
                    x_b = hcat(x_b,exact_x(v₀,x₀,b,τ,ω[end]))
                    v_b = hcat(v_b,exact_v(v₀,b,τ,ω[end]))
                    x = hcat(x,xₜₘₚ[:,2:end-1])
                    v = hcat(v,vₜₘₚ[:,2:end-1])
                    append!(τ_b,τ)
                    lvol = mod1(lvol+1,maxvol)
                    t_f = τ #For exit condition
                    break
                end
                i += 1
            end
        end
        k += 1
        if mod(k,2) == 0
            x_bp = x_bar(v_b[:,k-1],b₁[:,k-1],τ_b[k-1],ω[end-1])
            x_bm = x_bar(v_b[:,k],b₁[:,k],τ_b[k],ω[end])
            x̄ = hcat(x̄,average_orbit(x_bp,x_bm,τ_b[k-1],τ_b[k]))
        end
        # Set all params for next phase
        v₀ = v_b[:,end]
        x₀ = gc(x_b[:,end],v₀,Bfield[lvol](x_b[:,end],t))
        if t_f >= t[end] # If not finished yet
            τ_b[end] > t[end] ? τ_b=τ_b[1:end-1] : τ_b
            # If finished
            break
        end
    end
    return exact_particle(x,v,x_b,v_b,t,τ_b,x̄)
end


function _analytic_solve_nocrossing(x₀::Vector{T},v₀::Vector{T},t::Vector{T},Bf::Function) where T
    ω = [abs(q)/m * norm(Bf(x₀,t[1]))]
    x = x_b = x₀ - guiding_center(v₀,Bf(x₀,t[1]))
    v = v_b = v₀

    b₁ = Array{T,2}(undef,3,0)
    
    # If there are no crossings then solve the entire system at once
    B = Bf(x₀,t[1])
    b = magcoords(v₀,B)
    x = exact_x(v₀,x₀,b,t,ω[1])
    v = exact_v(v₀,b,t,ω[1])

    x̄ = Vector{T}(undef,0)
    τ_b = Vector{T}(undef,0)


    return exact_particle(x,v,x_b,v_b,t,τ_b,x̄)
end
#===
    Supporting functions
===#
"""
    exact_x(v₀,X₀,b,tᵢ,ω)
"""
function exact_x(v₀,X₀,b,tᵢ,ω) #position
    x = v_para(v₀,b)*tᵢ' + 1/ω * norm(v₀)*(-b[3] * sin.(ω*tᵢ)' + b[2] * cos.(ω*tᵢ)') .+ X₀
    return x
end
"""
    exact_v(v₀,b,tᵢ,ω)
"""
function exact_v(v₀,b,tᵢ,ω)
    v = v_para(v₀,b) .- norm(v₀)*(b[3] * cos.(ω*tᵢ)' + b[2] * sin.(ω*tᵢ)')
    return v
end
"""
    v_para(v₀,b)
"""
function v_para(v₀,b) #parallel velocity
    return dot(v₀,b[1])*b[1]
end
"""
    gc(x₀,v₀,B)
"""
function gc(x₀,v₀,B) #GC position
    return x₀ + cross(v₀,B)/norm(B,2)^2
end

### Functions for computing the average drift ###
"""
    x_bar(v₀,b,τ,ω)
"""
function x_bar(v₀,b,τ,ω)
    # ONLY WORKS SINCE ω=1
    x_bar = dot(v₀,b)*b*τ + 2/ω*dot(v₀,[0,0,1])*cross([0,0,1],b)
    return x_bar
end
"""
    average_orbit(x_bp,x_bm,τ_p,τ_m)
"""
function average_orbit(x_bp,x_bm,τ_p,τ_m)
    # Adds the average drift to the original position for adding to plots
    x_bar = (x_bp + x_bm)/(τ_p + τ_m)
    return x_bar
end
"""
    magcoords(v,B)
"""
function magcoords(v,B)
    # Compute the field aligned coordinates
    B_0     = norm(B,2)
    b1      = B/B_0
    b2      = cross(b1,v)/norm(v,2)
    b3      = cross(b1,b2)
    return [b1,b2,b3]
end
"""
    bound_time(ω,b)
"""
function bound_time(ω,b)
    # Compute the boundary crossing time
    τ = 2/ω * acot(-b[2][3]/b[3][3])
    if τ < 0
        τ += 2*pi
    end
    return τ
end


