module OrbitEqns

    using LinearAlgebra
    using Roots

    #=
        Full orbit and gyrocenter equations
    =#

    function exactsolve(x₀,v₀,t,B;q=1,m=1)

        function Bfield(x)
            α = pi/6
            if x > 0
                Bf = [cos(α),-sin(α),0]
            else
                Bf = [cos(α),sin(α),0]
            end
            return Bf
        end

        ω = q/m * norm(B(x₀)) #Specific to |B|=1

        dc = -1 #Starting field

        x = x₀
        v = v₀
        x_b = x₀
        v_b = x₀

        tn = length(t)
        τᵦ = []

        while true

            B = Bfield(xᵢ)
            b = magcoords(v₀,B)

            τ = bound_time(ω,b)
            if τ == 0
                # Iterate to boundary
            else
                # Store crossing time
                append!(τᵦ,τ)

                #Get all values of t in this section
                tᵢ = findall(x->sum(τᵦ[1:end-1])<=x<=sum(τᵦ),t) 
                # zero the crossing time
                t_sub = t(tᵢ) .- sum(τᵦ[1:end])
                
                # Eqns of motion
                xᵢ = p_x(v₀,x₀,b,tᵢ)
                vᵢ = p_v(v₀,b,tᵢ)
                # Store positions for returning
                x = hcat(x_b,xᵢ)
                v = hcat(v_b,vᵢ)
                # Positions at boundary
                x_b = hcat(x_b,p_x(v₀,x₀,b,τ))
                v_b = hcat(v_b,p_v(v₀,b,τ))

            end

            # Set all params for next phase
            v₀ = v_b
            x₀ = gc(x_b[:,end],v₀,B)
            dc = -dc #switch field

            if tᵢ[end] == t[end]
                break
            end

        end

    end

    function p_x(v₀,X₀,b,tᵢ) #position
        x = v_para(v₀,b).*tᵢ' + 1/ω * norm(v₀)*(-b[3].*sin(ω*tᵢ)' + b[2].*cos(ω*tᵢ)') + X₀
        return x
    end
    function p_v(v₀,b,tᵢ)
        v = v_para(v₀,b) - norm(v₀)*(b[3].*cos(ω*tᵢ)' + b[2].*sin(ω*tᵢ)')
        return v
    end
    function v_para(v₀,b) #parallel velocity
        return dot(v₀,b[1])*b[1]
    end
    function gc(x₀,v₀,B) #GC position
        return x₀ + cross(v₀,B)/norm(B,2)^2
    end


    function position(Bfield,x₀,v₀,t,Δt;usetdc=[],q=1,m=1)
        #=
        B       - Magnetic field (function)
        x₀     - initial position (guiding center)
        v₀     - initial velocity
        t       - array of times
        =#
        tn = length(t)
        t_out = zeros(tn)
        x = zeros(3,tn)
        v = zeros(3,tn)
        gr = zeros(tn)
        xdc = zeros(3)
        vdc = v₀
        Δx_dc = zeros(3)
        tdc = []
        t_t = 0
        Δx_o = 0 #For average drift
        τ_o = 0 #For crossing time
        avetraj = zeros(3)
        τᵦ = [] #For crossing time
        
        # Move particle from GC to initial position
        dc = sign(x₀[3])*ones(3) #Set current field to use
        if dc[3] == 0
            dc = -1*ones(3)
        end

        # Initial mag field
        B       = Bfield(dc)
        B_0     = norm(B,2)
        b       = magcoords(v₀,B)

        crossing = 1

        i = 1
        k = 1

        # Shift from GC to FO position
        # x₀ = x₀ - cross(B,v₀)/norm(B,2)^2
        
        # Gyrofreqency
        # ω   = q/(m*B_0)
        ω = 1.


        while k <= size(x)[2]

            tᵢ = t[i] - t_t
            
            # Position
            x[:,k] = v_par(v₀,b)*tᵢ + x_per(v₀,ω,b,tᵢ) + x₀
            t_out[k] = t[i]
            v[:,k] = vel(v_par(v₀,b),v₀,ω,b,tᵢ)

            gr[k] = norm(v[:,k] - v_par(v₀,b),2)

            if isempty(usetdc)    
                # If the crossing times are supplied
                if (k != 1) && (crossing == 1) && (x[3,k-1]*x[3,k] < 0) && (crossing != 2)
                    crossing = 0
                    append!(τᵦ,bound_time(ω,b))
                end
            elseif !isempty(usetdc) & (i > 1)
                # If the crossing times are not supplied
                if !isempty(findall(x->x==t[i],usetdc))
                    crossing = 0
                    if bound_time(ω,b) != 0
                        append!(τᵦ,bound_time(ω,b))
                    end
                end
                # if length(tdc)==1
                #     crossing = 0
                    # append!(τᵦ,tdc[1])
                # end
            end

            if crossing == 2
                crossing = 1
            end

            if (crossing == 0) && isempty(usetdc)
                # Find the intercept
                bifn(t) = (dot(v₀,b[1])*b[1]*t + 1/ω * norm(v₀)* (-sin(ω*t)*b[3] + cos(ω*t)*b[2]) + x₀)[3]
                tmpt = find_zero(bifn, (tᵢ-abs(t[i]-t[i-1]), tᵢ), Bisection(),atol=1.e-15,rtol=0)
                # Compute new velocity at intercept time
                v₀ = vel(v_par(v₀,b),v₀,ω,b,tmpt)
                # Compute position at intercept time
                tmpx = v_par(v₀,b)*tmpt + x_per(v₀,ω,b,tmpt) + x₀
                # Compute new magnetic field
                dc = -dc
                B = Bfield(dc)
                B_0 = norm(B,2)
                b = magcoords(v₀,B)
                # Add an extra column for intercept, move forward and add to arrays
                x[:,k] = tmpx
                v[:,k] = v₀
                t_out[k] = tmpt + t_t
                # Update guiding center
                x₀ = tmpx + cross(v₀,B)/norm(B,2)^2
                t_t = t_t + tmpt
                # Record the exact boundary
                xdc = hcat(xdc,tmpx)
                append!(tdc,tmpt)
                # Screw with counter
                crossing = 2
            elseif (crossing == 0) && !isempty(usetdc)
                dc = -dc
                B = Bfield(dc)
                B_0 = norm(B,2)

                v₀ = v[:,k]
                x₀ = x[:,k] + cross(v₀,B)/norm(B,2)^2

                b=magcoords(v₀,B)
                
                t_t = t[i]
                xdc = hcat(xdc,x[:,k])
                append!(tdc,tᵢ)

                crossing = 2
            end
            if (crossing == 2)
                if length(tdc)==1
                    append!(τᵦ,tdc[1])
                end
                vdc = hcat(vdc,v₀)
                Δx_dc = hcat(Δx_dc,x_bar(vdc[:,end-1],magcoords(vdc[:,end-1],Bfield(-dc)),τᵦ[end]))

                # Update the average drift
                # τ_n = τ_o
                # τ_o = bound_time(ω,b)
                # Δx_n = Δx_o
                # Δx_o = x_bar(v₀,b,τᵦ[end])
                if (mod(length(tdc),2) == 0) & (length(tdc)>0)
                    # Δx_n = x_bar(vdc[:,end],magcoords(vdc[:,end],Bfield(dc)),τᵦ[end])
                    # Δx_o = x_bar(vdc[:,end-1],magcoords(vdc[:,end-1],Bfield(-dc)),τᵦ[end-1])
                    # println(Δt,' ',τᵦ[end-1],' ',τᵦ[end],' ',Δx_dc[:,end-1],' ',Δx_dc[:,end])
                    avetraj = hcat(avetraj,average_orbit(Δx_dc[:,end-1],Δx_dc[:,end],τᵦ[end],τᵦ[end-1]))
                end
            end
            k = k + 1
            i = i + 1
        end
        p = Dict("x"=>x, "v"=>v, "t"=>t_out, "tdc"=>tdc, "xdc"=>xdc[:,2:end], "gr"=>gr,"avetraj"=>avetraj[:,2:end], "tau_b"=>τᵦ[2:end],"dx_dc"=>Δx_dc[:,2:end])
        return p
    end
    
    #=
    These are used to compute and plot the average drift of a particle
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
    
    function v_par(v₀,b)
        # parallel velocity
        return dot(v₀,b[1])*b[1]
    end

    function x_per(v₀,ω,b,t)
        # perpendicular position
        return 1/ω * norm(v₀)*(-sin(ω*t)*b[3] + cos(ω*t)*b[2])
    end

    function vel(v_par,v₀,ω,b,t)
        # velocity
        return v_par - norm(v₀)*(cos(ω*t)*b[3] + sin(ω*t)*b[2]) #Velocity
    end


    function magcoords(v,B)
        B_0     = norm(B,2)
        b1      = B/B_0
        b2      = cross(b1,v)/norm(v,2)
        b3      = cross(b1,b2)
        
        return [b1,b2,b3]
    end

    function bound_time(ω,b)
        τ = 2/ω * acot(-b[2][3]/b[3][3])
        if τ < 0 
            τ += 2*pi
        end
        # println(τ)
        return τ
    end


end