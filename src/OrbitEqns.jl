#=
    Full orbit and gyrocenter equations
=#

function analytic_solve(sim::analytic_sim,Bfield)
    for i = 1:sim.nparts
        analytic_solve(sim.sp[i].x,sim.sp[i].v,sim.sp[i].t,Bfield)
    end
end

function analytic_solve(x₀::Vector{Float64},v₀::Vector{Float64},t::Vector{Float64},Bfield::Union{Function,Array{Function}};crossing=true,eventfn=Nothing)
    # Works for static fields
    if typeof(Bfield) <: Function
        B = Bfield
    elseif typeof(Bfield) <: Array
        B = Bfield[1]
    end

    ω = abs(q)/m * norm(B(x₀)) #Specific to |B|=1

    x = x_b = x₀
    v = v_b = v₀

    n = length(t)
    τ_b = []

    if !crossing
        # If there are no crossings then solve the entire system at once
        B = Bfield(x₀)
        b = magcoords(v₀,B)
        x = exact_x(v₀,x₀,b,t,ω)
        v = exact_v(v₀,b,t,ω)
    end

    while crossing
        # If there are crossings solve per region
        B = Bfield(xᵢ)
        b = magcoords(v₀,B)
        ω = abs(q)/m * norm(B(x₀)) #Specific to |B|=1
        τ = bound_time(ω,b)
        if τ == 0
            # If crossing time cannot be computed the iterate
            i = 1
            while true
                x_tmp = exact_x(v₀,X₀,b,t[i],ω)
                v_tmp = exact_v(v₀,b,t[i],ω)
                chk = eventfn(x_tmp)
                if chk
                    ex(t) = exact_x(v₀,X₀,b,t,ω)
                    τ = bisection(ex,a=t[i-1],b=t[i+1])
                    break
                end
                i += 1
            end
        else
            # Store crossing time
            append!(τ_b,τ)
            #Get all values of t in this section
            tᵢ = findall(x->sum(τ_b[1:end-1])<=x<=sum(τ_b),t) 
            # zero the crossing time
            t_sub = t(tᵢ) .- sum(τ_b[1:end])
            # Eqns of motion
            xᵢ = exact_x(v₀,x₀,b,tᵢ,ω)
            vᵢ = exact_v(v₀,b,tᵢ,ω)
            # Store positions for returning
            x = hcat(x_b,xᵢ)
            v = hcat(v_b,vᵢ)
            # Positions at boundary
            x_b = hcat(x_b,exact_x(v₀,x₀,b,τ,ω))
            v_b = hcat(v_b,exact_v(v₀,b,τ,ω))

        end

        # Set all params for next phase
        v₀ = v_b
        x₀ = gc(x_b[:,end],v₀,B)

        if tᵢ[end] == t[end]
            # If finished
            break
        end

    end
    return exact_particle(x,v,x_b,v_b,t,τ_b)
end

#=
    Supporting functions
=#
function exact_x(v₀,X₀,b,tᵢ,ω) #position
    x = v_para(v₀,b)*tᵢ' + 1/ω * norm(v₀)*(-b[3] * sin.(ω*tᵢ)' + b[2] * cos.(ω*tᵢ)') .+ X₀
    return x
end
function exact_v(v₀,b,tᵢ,ω)
    v = v_para(v₀,b) .- norm(v₀)*(b[3] * cos.(ω*tᵢ)' + b[2] * sin.(ω*tᵢ)')
    return v
end
function v_para(v₀,b) #parallel velocity
    return dot(v₀,b[1])*b[1]
end
function gc(x₀,v₀,B) #GC position
    return x₀ + cross(v₀,B)/norm(B,2)^2
end

#=
    Functions for computing the average drift
=#

function x_bar(v₀,b,τ)
    # ONLY WORKS SINCE ω=1
    x_bar = dot(v₀,b[1])*b[1]*τ + 2*dot(v₀,[0,0,1])*cross([0,0,1],b[1])
    return x_bar
end

function average_orbit(x_bp,x_bm,τ_p,τ_m)
    # Adds the average drift to the original position for adding to plots
    x_bar = (x_bp + x_bm)/(τ_p + τ_m)
    return x_bar
end

function magcoords(v,B)
    # Compute the field aligned coordinates
    B_0     = norm(B,2)
    b1      = B/B_0
    b2      = cross(b1,v)/norm(v,2)
    b3      = cross(b1,b2)
    return [b1,b2,b3]
end

function bound_time(ω,b)
    # Compute the boundary crossing time
    τ = 2/ω * acot(-b[2][3]/b[3][3])
    if τ < 0
        τ += 2*pi
    end
    return τ
end

