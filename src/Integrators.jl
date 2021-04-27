# Integrators




@everywhere function integrate(p::Union{particle,Array{particle}},t_f;eventfn=Nothing,method="RK4")

    if method=="RK4"
        integrator = RK4
    end

    while p.t[end] < t_f
        xv = vcat(p.x[:,end], p.v[:,end])
        tᵢ = [p.t[end], p.t[end]+Δt]
        B = p.Bfield[lvol]
        
        fₓ(x,t) = p.dvdt(x,t,B)
        
        x_tmp, t_tmp, Δt_tmp, crossing = integrator(fₓ,xv,tᵢ,eventfn=eventfn)
        
        # Storage
        p.x = hcat(p.x,x_tmp[1:3,2])
        p.v = hcat(p.v,x_tmp[4:6,2])
        p.t = vcat(p.t,t_tmp[2])
        p.B = hcat(p.B,p.Bfield(x_tmp[1:3,2]))
        # Updates the volume based on the direction of crossing

        if !crossing
            append!(p.lvol,p.lvol[end])
        else
            append!(p.lvol,p.lvol[end]+sign(v_tmp[3,2]))
        end
        
    end

    return p
end


#== INTEGRATORS ==#
function RK4(fₓ::Function,x₀::Vector{Float64},t::Vector{Float64};eventfn=Nothing,tol=1.e-15)
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
            if event == 0
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
    
    xtilde(τ) = boostrapped(τ,h,x[3,:],k[1,3],f₁,fₘ) - bound
    xtilde_dt(τ) = dt_bootstrapped(τ,h,x[3,:],k[1,1],f₁,fₘ) - bound

    f_root!(xtilde,xtilde_dt,τₙ)
    
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

function f_root!(f,df,x₀;maxits=20,atol=1.e-10)
    x₀,dx,converged = Newtons(f,df,x₀)
    if !converged
        # If Newtons oscillates revert to bisection
        x₀ = bisection(f,a=x₀-2*dx,b=x₀+2*dx)
    end
    if x₀ ∉ [0,1]
        # If Newtons fails call bisection
        x₀ = bisection(f)
    end
    return x₀
end

function Newtons(f,df,x₀;maxits=20,atol=1.e-10)
    xₒ = x₀
    # if df == Nothing
    converged = true
    for n = 1:maxits
        xₒ = xₙ
        xₙ = xₒ - f(xₒ)/df(xₒ)
        dx = abs(xₒ - xₙ)
        if (dx < atol)
            break
        end
    end
    if (n == maxits) & (dx > atol)
        converged = false
    end
    return xₙ, dx, converged
end

function bisection(f::Function;maxits=30,atol=1.e-10,a=0.,b=1.)
    
    converged = true
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
            break
        end
    end
    if i == maxits
        converged = false
    end
    return c, converged
end
