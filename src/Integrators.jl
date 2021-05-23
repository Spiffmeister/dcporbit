# Integrators

#=
    ORBIT SOLVER
=#
function solve_orbit!(p::particle,ODE,t_f;method="RK4")

    # Used as an iterator this will update the particle object provided

    # Set the integrator
    if method=="RK4"
        integrator = RK4
    end
    # Check the magnetic field type
    Btype = typeof(ODE.MagneticField)

    crossing = false

    # if !p.gc_or_fo
    #     EOM(x,t,B) = ODE.EOM(x,t,B)
    # else
    #     EOM(x,t,B) = ODE.EOM_GC(x,t,B)
    # end

    
    # Main time loop
    while p.t[end] < t_f
        xv = vcat(p.x[:,end], p.v[:,end])
        tᵢ = [p.t[end], p.t[end]+p.Δt]
        if Btype <: Vector{Function}
            B = ODE.MagneticField[p.lvol[end]]
        else
            B = ODE.MagneticField
        end
        
        fₓ(x,t) = ODE.EOM(x,t,B(p.x[end,:]))
        
        if !crossing
            xv, t, Δt_tmp, crossing = integrator(fₓ,xv,tᵢ,eventfn=ODE.event)
        else
            xv, t, Δt_tmp, crossing = integrator(fₓ,xv,tᵢ)
        end
        # Storage
        p.x = hcat(p.x,xv[1:3,2])
        p.v = hcat(p.v,xv[4:6,2])
        p.t = vcat(p.t,t[2])
        p.B = hcat(p.B,B(xv[1:3,2]))
        p.gc = hcat(p.gc,xv[1:3,2]+guiding_center(xv[1:3,2],xv[4:6,2],B))
        # Updates the volume based on the direction of crossing
        if !crossing
            append!(p.lvol,p.lvol[end])
        else
            append!(p.lvol,p.lvol[end]+sign(xv[6,2]))
        end
        
    end

end



function run_sim!(f::sim,ODE,t_f;method="RK4")
    # INTERFACE FOR solving simulations
    for i = 1:f.nparts
        solve_orbit!(f.sp[i],ODE,t_f,method=method)
    end
    # return f
end


function integrate!(ODE::Function,xv::Vector{Float64},t;eventfn=nothing,method="RK4")

    if method=="RK4"
        integrator = RK4
    end

    # Turn event location on/off
    typeof(eventfn) <: Nothing ? dchk = false : dchk = true

    while i < m
        x, t, h, crossing = integrator(ODE,xv,t)

        if dchk
            event = eventfn(x[:,i:i+1])
            if event < 0.
                t, h, xn = event_loc(fₓ,x[:,i:i+1],k,h,t)
                t[i+1] = t[i] + h
                x[:,i+1] = xn
                crossing = true
            elseif event == 0
                crossing = true
            end
        end

    end

    return x, t, Δt, crossing
end


#== INTEGRATORS ==#
function RK4(fₓ::Function,x₀::Vector{Float64},t::Vector{Float64};eventfn=nothing,tol=1.e-15)
    # Standard RK4 with discontinuity detection

    s = 4 #Stages

    n = length(x₀)
    m = length(t)
    x = zeros(n,m)
    k = zeros(n,s)
    h = 0.
    x[:,1] = x₀

    aᵢⱼ = [[1/2] , [0 , 1/2] , [0 , 0 , 1]]
    bᵢ = [1/6, 1/3, 1/3, 1/6]
    cᵢ = [1/2, 1/2, 1]

    i = 1

    crossing = false
    dchk = false
    if typeof(eventfn)!=Nothing
        dchk = true
    end

    while i < m
        # Integrate
        h = abs(t[i] - t[i+1])
        k[:,1] = fₓ(x[:,i],t[i])
        for j = 1:(s-1)
            k[:,j+1] = fₓ(x[:,i]+h*k[:,1:j]*aᵢⱼ[j],t[i]+h*cᵢ[j])
        end
        x[:,i+1] = x[:,i] + h*k*bᵢ

        if dchk
            event = eventfn(x[:,i:i+1])
            if event < 0.
                t, h, xn = event_loc(fₓ,x[:,i:i+1],k,h,t)
                t[i+1] = t[i] + h
                x[:,i+1] = xn
                crossing = true
            elseif event == 0
                crossing = true
            end
        end
        i += 1
    end
    return x, t, h, crossing
end

#=
Tools for integrators
=#

function event_loc(dvdt::Function,x::Array{Float64},k::Array{Float64},h::Float64,t::Array{Float64};bound=0.,tol=1.e-14)
    # Endpoint derivative evaluation
    f₁ = dvdt(x[:,2],t[2])
    # Cubic hermite for bootstrapping
    # Evaluate the derivative at the h/10 point
    fₘ = dvdt(cubic_hermite(1/10,h,x,k,f₁),t[1] + h/10.)
    
    # Initial guess based on linear interpolation
    τₙ = (bound - x[3,1])/(x[3,2] - x[3,1])
    
    xtilde(τ) = boostrapped(τ,h,x,k[:,1],f₁,fₘ)[3] #- bound
    xtilde_dt(τ) = dt_bootstrapped(τ,h,x,k[:,1],f₁,fₘ)[3] #- bound
    
    τₙ = f_root!(xtilde,xtilde_dt,τₙ)
    xₙ = boostrapped(τₙ,h,x,k[:,1],f₁,fₘ)
    h = τₙ*h
    t[2] = t[1] + h
    
    return t, h, xₙ
end

function cubic_hermite(τ::Float64,h::Float64,x::Array{Float64},k::Array{Float64},f₁::Array{Float64})
    bar_x = (1-τ)*x[:,1] + τ*x[:,2] + τ*(τ-1)*((1-2*τ)*(x[:,2]-x[:,1]) + (τ-1)*h*k[:,1] + τ*h*f₁)
    return bar_x
end

function boostrapped(τ::Float64,h::Float64,x::Array{Float64},k₁::Array{Float64},f₁::Array{Float64},fₘ::Array{Float64})
    # Bootstrapped to quartic
    d₀ = (τ-1)^2 * (15/4 * τ^2 + 2*τ + 1)
    d₁ = τ*(τ-1)^2 * (1-τ*35/8)
    d₂ = 1-d₀
    d₃ = τ^2 * (τ-1)*(85*τ-13)/72
    d₄ = τ^2 * (τ-1)^2 * 125/18
    tilde_x = d₀*x[:,1] + d₁*h*k₁ + d₂*x[:,2] + d₃*h*f₁ + d₄*h*fₘ
    return tilde_x
end

function boostrapped(τ::Float64,h::Float64,x::Array{Float64},k₁::Float64,f₁::Float64,fₘ::Float64)
    # Bootstrapped to quartic
    d₀ = (τ-1)^2 * (15/4 * τ^2 + 2*τ + 1)
    d₁ = τ*(τ-1)^2 * (1-τ*35/8)
    d₂ = 1-d₀
    d₃ = τ^2 * (τ-1)*(85*τ-13)/72
    d₄ = τ^2 * (τ-1)^2 * 125/18
    tilde_x = d₀*x[1] + d₁*h*k₁ + d₂*x[2] + d₃*h*f₁ + d₄*h*fₘ
    return tilde_x
end

function dt_bootstrapped(τ::Float64,h::Float64,x::Array{Float64},k₁::Array{Float64},f₁::Array{Float64},fₘ::Array{Float64})

    d₀ = 2*(τ-1)*(15/4*τ^2 + 2*τ + 1) + (τ-1)^2 * (15*τ/2 + 2)
    d₁ = (τ-1)^2 * (1-35*τ/8) + 2*τ*(τ-1) * (1-35*τ/8) - 35*τ*(τ-1)^2 / 8
    d₂ = -1*d₀
    d₃ = τ*(τ-1)*(85*τ-13)/36 + τ^2 * (85*τ-13)/72 + 85*τ^2 * (τ-1)/72
    d₄ = 125*τ*(τ-1)^2 / 9 + 125*τ^2 * (τ-1)/9
    dtilde_xdt = d₀*x[:,1] + d₁*h*k₁ + d₂*x[:,2] + d₃*h*f₁ + d₄*h*fₘ
    return dtilde_xdt
end

function dt_bootstrapped(τ::Float64,h::Float64,x::Array{Float64},k₁::Float64,f₁::Float64,fₘ::Float64)

    d₀ = 2*(τ-1)*(15/4*τ^2 + 2*τ + 1) + (τ-1)^2 * (15*τ/2 + 2)
    d₁ = (τ-1)^2 * (1-35*τ/8) + 2*τ*(τ-1) * (1-35*τ/8) - 35*τ*(τ-1)^2 / 8
    d₂ = -1*d₀
    d₃ = τ*(τ-1)*(85*τ-13)/36 + τ^2 * (85*τ-13)/72 + 85*τ^2 * (τ-1)/72
    d₄ = 125*τ*(τ-1)^2 / 9 + 125*τ^2 * (τ-1)/9
    dtilde_xdt = d₀*x[1] + d₁*h*k₁ + d₂*x[2] + d₃*h*f₁ + d₄*h*fₘ
    return dtilde_xdt
end

#=
=#

function f_root!(f,df,x₀;maxits=20,atol=1.e-14)
    x₀,converged = newtons!(f,df,x₀,maxits=maxits,atol=atol)
    if !converged
        # If newtons oscillates revert to bisection
        x₀,converged = bisection(f,a=x₀-2*dx,b=x₀+2*dx)
    end
    if !(0. < x₀ < 1.)
        # If newtons fails call bisection
        x₀,converged = bisection(f)
    end
    return x₀
end

function newtons!(f,df,xₙ;maxits=20,atol=1.e-14)
    # xₙ = x₀
    # if df == Nothing
    converged = false
    for n = 1:maxits
        xₒ = xₙ
        xₙ = xₒ - f(xₒ)/df(xₒ)
        if (abs(xₒ - xₙ) < atol)
            converged = true
            break
        end
    end
    return xₙ, converged
end

function bisection(f::Function;maxits=50,atol=1.e-10,a=0.,b=1.)
    
    converged = false
    c = 0.
    if sign(f(a))==sign(f(b))
        return Nothing, false
    end

    for i = 1:maxits
        c = (b-a)/2.
        if sign(f(c)) == sign(f(a))
            a = c
        else
            b = c
        end
        if (b-a) < atol
            converged = true
            break
        end
    end
    return c, converged
end
