# Integrators


#=
"""
    RungeKutta{TT,STAGES}
"""
struct RungeKutta{TT,STAGES}
    k   :: Vector{Vector{TT}}
    s   :: Int64

    aᵢⱼ :: Vector{Vector{TT}}
    bᵢ  :: Vector{TT}
    cᵢ  :: Vector{TT}
end
"""
    RungeKutta{TT}(order,len)
"""
function RungeKutta{TT}(order::Int64,len::Int64)
    if order == 4
        aᵢⱼ = [[1/2.] , [0. , 1/2.] , [0. , 0. , 1.]]
        bᵢ = [1/6., 1/3., 1/3., 1/6.]
        cᵢ = [1/2., 1/2., 1.]

        k = [zeros(TT,len) for i = 1:s]
        s = order

        return RungeKutta{TT,order}(k,s,aᵢⱼ,bᵢ,cᵢ)
    end
end

"""
    RK4(p,f,t)
"""
function (RK4::RungeKutta{TT,4})(p::particle,f::forces,t:TT,h::TT) where {TT,S}
    x[:,1] = x₀

    i = 1

    # Integrate
    RK.k[:,1] = fₓ(x[:,i],t[1])
    for j = 1:(s-1)
        RK.k[:,j+1] = fₓ(x[:,i]+h*RK.k[:,1:j]*RK.aᵢⱼ[j],t[i]+h*RK.cᵢ[j])
    end
    x[:,i+1] = x[:,i] + h*RK.k*RK.bᵢ
    h = t[i+1] - t[i]
    i += 1

    p
end
=#

