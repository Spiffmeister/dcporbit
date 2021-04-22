# Integrators




function integrate(p::Union{particle,Array{particle}},t_f;eventfn=Nothing,method="RK4")

    if method=="RK4"
        integrator = RK4
    end

    while p.t[end] < t_f
        xv = vcat(p.x[:,end], p.v[:,end])
        tᵢ = [p.t[end], p.t[end]+Δt]
        B = p.Bfield[lvol]
        
        f_x(x,t) = p.dvdt(x,t,B)
        
        x_tmp, t_tmp, Δt_tmp = integrator(f_x,xv,tᵢ,eventfn=eventfn)
        
        # Storage
        p.x = hcat(p.x,x_tmp[1:3,2])
        p.v = hcat(p.v,x_tmp[4:6,2])
        p.t = vcat(p.t,t_tmp[2])
        p.B = hcat(p.B,p.Bfield(x_tmp[1:3,2]))
        # Updates the volume based on the direction of crossing
        append!(p.lvol,p.lvol[end]+sign(v_tmp[3,2]))
        
    end

    return p
end


#== INTEGRATORS ==#
function RK4(f_x::Function,x_0::Array{Float64},t::Vector{Float64};eventfn=Nothing,tol=1.e-15)
    # Standard RK4 with discontinuity detection

    s = 4 #Stages

    n = length(x_0)
    m = length(t)
    x = zeros(n,m)
    k = zeros(n,s)
    h = 0
    x[:,1] = x_0

    a_ij = [[1/2] , [0 , 1/2] , [0 , 0 , 1]]
    b_i = [1/6, 1/3, 1/3, 1/6]
    c_i = [1/2, 1/2, 1]

    i = 1

    crossing = 0
    dchk = false
    if typeof(eventfn)!=Nothing
        dchk = true
    end

    while i < m
        # Integrate
        h = abs(t[i] - t[i+1])
        k[:,1] = f_x(x[:,i],t[i])
        for j = 1:(s-1)
            k[:,j+1] = f_x(x[:,i]+h*k[:,1:j]*a_ij[j],t[i]+h*c_i[j])
        end
        x[:,i+1] = x[:,i] + h*k*b_i

        if dchk
            event = eventfn(x[:,i:i+1])
            if event < 0.
                t, h, xn = event_loc(f_x,x[:,i:i+1],k,h,t)
                t[i+1] = t[i] + h
                x[:,i+1] = xn
                crossing = 1
            if event == 0
                crossing = 1
            end
        end
        i += 1
    end
    return x, t, h, crossing
end

#=
Tools for integrators
=#

function event_loc(dvdt,x::Array{Float64},k::Array{Float64},h::Float64,t::Array{Float64};bound=0,tol=1.e-14)
    # Endpoint derivative evaluation
    f_1 = dvdt(x[:,2],t[2])
    # Cubic hermite for bootstrapping
    cubic_hermite(τ)    = (1-τ)*x[:,1] + τ*x[:,2] + τ*(τ-1)*((1-2*τ)*(x[:,2]-x[:,1]) + (τ-1)*h*k[:,1] + τ*h*f_1)
    # Evaluate the derivative at the h/10 point
    f_m = dvdt(cubic_hermite(1/10),t[1] + h/10)

    # Bootstrapped quartic polynomials
    function boostrapped(τ::Float64,h::Float64,x::Array{Float64},k1::Array{Float64},f_1::Array{Float64},f_m::Array{Float64})

        d_0 = (τ-1)^2 * (15/4 * τ^2 + 2*τ + 1)
        d_1 = τ*(τ-1)^2 * (1-τ*35/8)
        d_2 = 1-d_0
        d_3 = τ^2 * (τ-1)*(85*τ-13)/72
        d_4 = τ^2 * (τ-1)^2 * 125/18
        tilde_x = d_0*x[:,1] + d_1*h*k1 + d_2*x[:,2] + d_3*h*f_1 + d_4*h*f_m
        return tilde_x
    end
    function dt_bootstrapped(τ::Float64,h::Float64,x::Array{Float64},k1::Array{Float64},f_1::Array{Float64},f_m::Array{Float64})

        d_0 = 2*(τ-1)*(15/4*τ^2 + 2*τ + 1) + (τ-1)^2 * (15*τ/2 + 2)
        d_1 = (τ-1)^2 * (1-35*τ/8) + 2*τ*(τ-1) * (1-35*τ/8) - 35*τ*(τ-1)^2 / 8
        d_2 = -1*d_0
        d_3 = τ*(τ-1)*(85*τ-13)/36 + τ^2 * (85*τ-13)/72 + 85*τ^2 * (τ-1)/72
        d_4 = 125*τ*(τ-1)^2 / 9 + 125*τ^2 * (τ-1)/9
        dtilde_xdt = d_0*x[:,1] + d_1*h*k1 + d_2*x[:,2] + d_3*h*f_1 + d_4*h*f_m
        return dtilde_xdt
    end
    
    # NEWTONS METHOD
    xtilde(τ)       = boostrapped(τ,h,x,k[:,1],f_1,f_m)[3]
    xtilde_dt(τ)    = dt_bootstrapped(τ,h,x,k[:,1],f_1,f_m)[3]
    
    # Initial guess based on linear interpolation
    τₙ = (bound - x[3,1])/(x[3,2] - x[3,1])

    x_n = boostrapped(τₙ,h,x,k[:,1],f_1,f_m)

    τₙ = f_root(xtilde,xtilde_dt)

    h = τₙ*h
    t[2] = t[1] + h

    return t, h, x_n
end



#=
=#

function f_root(f,df,x₀;maxits=20,atol=1.e-10)
    xₙ,converged = Newtons(f,df,x₀)
    if !converged
        # If Newtons oscillates revert to bisection
        xₙ = bisection(f,a=xₙ-2*dx,b=xₙ+2*dx)
    end
    if xₙ ∉ [0,1]
        # If Newtons fails call bisection
        xₙ = bisection(f)
    end
    return xₙ
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
    if n == maxits && dx > atol
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
