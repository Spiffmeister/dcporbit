"""
    event_loc(dvdt::Function,x::Array{Float64},k::Array{Float64},h::Float64,t::Array{Float64};bound=0.,tol=1.e-14)
Event location function using bootstrapping for higher order interpolation
"""
function event_loc(dvdt::Function,x::Array{Float64},k::Array{Float64},h::Float64,t::Array{Float64};bound=0.,tol=1.e-14)
    # Endpoint derivative evaluation
    f₁ = dvdt(x[:,2],t[2])
    # Evaluate the derivative at h/10 point
    fₘ = dvdt(cubic_hermite(1/10,h,x,k,f₁),t[1] + h/10.)
    
    # Newtons initial guess based on linear interpolation
    τₙ = (bound - x[3,1])/(x[3,2] - x[3,1])
    
    xtilde(τ) = boostrapped(τ,h,x,k[:,1],f₁,fₘ)[3] #- bound
    xtilde_dt(τ) = dt_bootstrapped(τ,h,x,k[:,1],f₁,fₘ)[3] #- bound
    
    τₙ = f_root!(xtilde,xtilde_dt,τₙ)
    xₙ = boostrapped(τₙ,h,x,k[:,1],f₁,fₘ)
    h = τₙ*h
    t[2] = t[1] + h
    
    return t, h, xₙ
end


#=
    Rootfinding methods
=#
"""
    f_root!(f,df,x₀;maxits=20,atol=1.e-14)
"""
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

"""
    newtons!(f,df,xₙ;maxits=20,atol=1.e-14)
"""
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

"""
    bisection(f::Function;maxits=50,atol=1.e-10,a=0.,b=1.)
"""
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


#= INTERPOLATING FUNCTIONS =#
"""
    cubic_hermite(τ::Float64,h::Float64,x::Array{Float64},k::Array{Float64},f₁::Array{Float64})
"""
function cubic_hermite(τ::Float64,h::Float64,x::Array{Float64},k::Array{Float64},f₁::Array{Float64})
    bar_x = (1-τ)*x[:,1] + τ*x[:,2] + τ*(τ-1)*((1-2*τ)*(x[:,2]-x[:,1]) + (τ-1)*h*k[:,1] + τ*h*f₁)
    return bar_x
end
"""
    boostrapped

Methods:
```julia
    boostrapped(τ::Float64,h::Float64,x::Array{Float64},k₁::Array{Float64},f₁::Array{Float64},fₘ::Array{Float64})
    boostrapped(τ::Float64,h::Float64,x::Array{Float64},k₁::Float64,f₁::Float64,fₘ::Float64)
```
"""
function boostrapped end
function boostrapped(τ::Float64,h::Float64,x::Array{Float64},k₁::Array{Float64},f₁::Array{Float64},fₘ::Array{Float64})
    # Vector function Bootstrapped to quartic
    d₀ = (τ-1)^2 * (15/4 * τ^2 + 2*τ + 1)
    d₁ = τ*(τ-1)^2 * (1-τ*35/8)
    d₂ = 1-d₀
    d₃ = τ^2 * (τ-1)*(85*τ-13)/72
    d₄ = τ^2 * (τ-1)^2 * 125/18
    tilde_x = d₀*x[:,1] + d₁*h*k₁ + d₂*x[:,2] + d₃*h*f₁ + d₄*h*fₘ
    return tilde_x
end
function boostrapped(τ::Float64,h::Float64,x::Array{Float64},k₁::Float64,f₁::Float64,fₘ::Float64)
    # Scalar function bootstrapped to quartic
    d₀ = (τ-1)^2 * (15/4 * τ^2 + 2*τ + 1)
    d₁ = τ*(τ-1)^2 * (1-τ*35/8)
    d₂ = 1-d₀
    d₃ = τ^2 * (τ-1)*(85*τ-13)/72
    d₄ = τ^2 * (τ-1)^2 * 125/18
    tilde_x = d₀*x[1] + d₁*h*k₁ + d₂*x[2] + d₃*h*f₁ + d₄*h*fₘ
    return tilde_x
end

"""
    dt_bootstrapped

Methods:
```julia
dt_bootstrapped(τ::Float64,h::Float64,x::Array{Float64},k₁::Array{Float64},f₁::Array{Float64},fₘ::Array{Float64})
dt_bootstrapped(τ::Float64,h::Float64,x::Array{Float64},k₁::Float64,f₁::Float64,fₘ::Float64)
```
"""
function dt_bootstrapped end
function dt_bootstrapped(τ::Float64,h::Float64,x::Array{Float64},k₁::Array{Float64},f₁::Array{Float64},fₘ::Array{Float64})
    # Derivative of vector function
    d₀ = 2*(τ-1)*(15/4*τ^2 + 2*τ + 1) + (τ-1)^2 * (15*τ/2 + 2)
    d₁ = (τ-1)^2 * (1-35*τ/8) + 2*τ*(τ-1) * (1-35*τ/8) - 35*τ*(τ-1)^2 / 8
    d₂ = -1*d₀
    d₃ = τ*(τ-1)*(85*τ-13)/36 + τ^2 * (85*τ-13)/72 + 85*τ^2 * (τ-1)/72
    d₄ = 125*τ*(τ-1)^2 / 9 + 125*τ^2 * (τ-1)/9
    dtilde_xdt = d₀*x[:,1] + d₁*h*k₁ + d₂*x[:,2] + d₃*h*f₁ + d₄*h*fₘ
    return dtilde_xdt
end

function dt_bootstrapped(τ::Float64,h::Float64,x::Array{Float64},k₁::Float64,f₁::Float64,fₘ::Float64)
    # Derivative of scalar function
    d₀ = 2*(τ-1)*(15/4*τ^2 + 2*τ + 1) + (τ-1)^2 * (15*τ/2 + 2)
    d₁ = (τ-1)^2 * (1-35*τ/8) + 2*τ*(τ-1) * (1-35*τ/8) - 35*τ*(τ-1)^2 / 8
    d₂ = -1*d₀
    d₃ = τ*(τ-1)*(85*τ-13)/36 + τ^2 * (85*τ-13)/72 + 85*τ^2 * (τ-1)/72
    d₄ = 125*τ*(τ-1)^2 / 9 + 125*τ^2 * (τ-1)/9
    dtilde_xdt = d₀*x[1] + d₁*h*k₁ + d₂*x[2] + d₃*h*f₁ + d₄*h*fₘ
    return dtilde_xdt
end
