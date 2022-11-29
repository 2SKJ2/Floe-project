# Module for the simulation of circular sea-ice floes with inelastic collisions and coupling to ocean and atmosphere forcings
# Dynamic coupling is based on the rheology of Herman et al. (2026)
# Thermodynamic coupling is based on

module FloeModule

    import Base.@kwdef
    using StructArrays
    using Statistics

    @kwdef mutable struct Floe
        # Floe struct containing positional and dimensional information of a given floe
        x::Float64 = 0. # x-position [m]
        y::Float64 = 0.  # y-position [m]
        θ::Float64 = 0. # angle wrt x-axis (+ve anticlockwise) [radians]
        u::Float64 = 0. # Mean x-velocity [m/s]
        v::Float64 = 0. # Mean y-velocity [m/s]
        ζ::Float64 = 0. # rotation rate [s-1]
        r::Float64 = 0. # radius [m]
        h::Float64 = 0. # thickness [m]
        ID::Int64 = 0 # Permanent floe ID
        cID::Int64 = 0 # ID of mirrored floe (for periodic boundary conditions)
        collision::Bool = false # True if this floe is under collision
    end

    @kwdef mutable struct FloeGrid
        # FloeGrid struct containing geometrical information about the surface of the domain
        size::Tuple{Int64,Int64} # Domain size
        extent::Tuple{Float64,Float64} # Domain extent
        xc::Vector{Float64} = Vector{Float64}(undef,size[1]) # Grid center coordinate in x
        yc::Vector{Float64} = Vector{Float64}(undef,size[2]) # Grid center coordinate in y
        xg::Vector{Float64} = Vector{Float64}(undef,size[1]+1) # Grid edge coordinate in x
        yg::Vector{Float64} = Vector{Float64}(undef,size[2]+1) # Grid edge coordinate in y
        function FloeGrid(size::Tuple{Int64,Int64},extent::Tuple{Float64,Float64},xc::Vector{Float64},yc::Vector{Float64},xg::Vector{Float64},yg::Vector{Float64})
            Δx = extent[1] / size[1] # Assumes regular grid (Δx=Δy)
            xc = (Δx/2):Δx:(extent[1]-Δx/2)
            yc = (Δx/2):Δx:(extent[2]-Δx/2)
            xg = 0:Δx:extent[1]
            yg = 0:Δx:extent[2]
            return new(size, extent, xc, yc, xg, yg)
        end
    end

    @kwdef mutable struct FloeModel
        # FloeModel struct containing all information pertaining to floes and their coupling
        floes::StructArray{Floe} # Floes
        mc::NamedTuple # Floe model constants
        grid::FloeGrid # Floe model constants
        Δt::Float64 # Floe model timestep (assumed constant)
        T_O::Array{Float64,2} = zeros(grid.size[1],grid.size[2]) # Surface ocean temperature [°C]
        S_O::Array{Float64,2} = zeros(grid.size[1],grid.size[2]) # Surface ocean salinity [psu]
        u_O::Array{Float64,2} = zeros(grid.size[1],grid.size[2]) # Surface ocean u-velocity [m/s]
        v_O::Array{Float64,2} = zeros(grid.size[1],grid.size[2]) # Surface ocean v-velocity [m/s]
        T_A::Array{Float64,2} = zeros(grid.size[1],grid.size[2]) # Surface atmopshere temperature [°C]
        u_A::Array{Float64,2} = zeros(grid.size[1],grid.size[2]) # Surface atmopshere u-velocity [m/s]
        v_A::Array{Float64,2} = zeros(grid.size[1],grid.size[2]) # Surface atmopshere v-velocity [m/s]
        taux::Array{Float64,2} = zeros(grid.size[1],grid.size[2]) # Net x-stress on ocean [Nm⁻²]
        tauy::Array{Float64,2} = zeros(grid.size[1],grid.size[2]) # Net y-stress on ocean [Nm⁻²]
        melt_rate::Array{Float64,2} = zeros(grid.size[1],grid.size[2]) # Melt rate [ms⁻¹]
        SIfract::Array{Float64,2} = zeros(grid.size[1],grid.size[2]) # Areal sea-ice fraction
        hflx::Array{Float64,2} = zeros(grid.size[1],grid.size[2]) # Net surface heat flux [Wm⁻²]
    end

    @inline function prepareCollisionFactors(E::Float64,ν::Float64,er::Float64)
        # Preparing constants for collisions according to Herman et al. (2016)
        #
        Ec = 0.5 * E / (1 - ν^2)
        Gc = 0.25 * E / (2 + ν) / (1 - ν)
        βc = log(er) / sqrt(log(er)^2 + π^2)
        return Ec, Gc, βc
    end

    @inline function getCollisionFactors(Ec::Float64,Gc::Float64,βc::Float64,rho_I::Float64,r1::Float64,r2::Float64,h1::Float64,h2::Float64)
        # Preparing constants for collisions according to Herman et al. (2016)
        if r2 == -1 # Wall collisions
            hm = h1
            mc = 1
        else
            hm = min(h1,h2)
            mc = π*rho_I*r1^2*r2^2 / (r1^2 + r2^2)
        end
        kn = π/4*hm*Ec
        kt = 6*Gc/Ec*kn
        γn = -βc*sqrt(5*kn*mc)
        γt = -2*βc*sqrt(5*Gc*kn*mc/Ec)
        return kn,kt,γn,γt
    end

    @inline function FloeShape(h::Float64,k::Float64,r::Float64)
        θ = LinRange(0,2*π,200)
        h .+ r*sin.(θ), k .+ r*cos.(θ)
    end

    @inline function FloeOrient(x::Float64,y::Float64,r::Float64,θ::Float64)
        l = LinRange(0,r,100)
        l*cos(θ) .+ x, l*sin(θ) .+ y
    end

    @inline function initFloes!(ngrid::Float64,Lx::Float64,Ly::Float64,Δx::Float64,rmin::Float64,rmax::Float64,hmin::Float64,hmax::Float64,floes::Float64)
        # Initialize floe in a regular lattice
        xxc = Lx/ngrid/2
        yyc = Ly/2 + Ly/ngrid/2
        nfloes = length(floes)
        for i in 1:nfloes
            floes.x[i] = xxc
            floes.y[i] = yyc
            floes.r[i] = rmin + (rmax-rmin)*rand()
            floes.h[i] = hmin + (hmax-hmin)*rand()
            floes.ID[i] = i
            xxc += Lx/ngrid
            if xxc>Lx
                xxc = Lx/ngrid/2
                yyc += Ly/ngrid
            end
        end
        return floes
    end

    @inline function atan2(y::Float64,x::Float64)
        return sign(x)^2 * atan(y/x) + (1 - sign(x))/2*(1 + sign(y) - sign(y)^2)*π
    end

    @inline function collideTwoFloes!(Fx::Vector{Float64}, Fy::Vector{Float64}, Trq::Vector{Float64}, idx1, f::NamedTuple, mc::NamedTuple)
        # Collision routine modifying Fx, Fy, and Trq
        # idx1 is the index of the primary floe being collided
        # f contains relevant information about the two floes being collided

        # Getting collision constants (not using kt)
        kn, kt, γn, γt = getCollisionFactors(mc.Ec,mc.Gc,mc.βc,mc.rho_I,f.r1,f.r2,f.h1,f.h2)

        # Rotating coordinate system
        Q = [cos(f.theta) sin(f.theta) ; -sin(f.theta) cos(f.theta)] # Rotation matrix
        vel1 = Q*[f.u1 f.v1]' # Initial velocity of floe 1 in new coord system. First component normal, second tangential
        vel2 = Q*[f.u2 f.v2]' # Initial velocity of floe 2 in new coord system. First component normal, second tangential

        # Normal force (based on Herman et al. 2016 - Supplementary Information)
        δndt = vel1[1] - vel2[1] # Rate of change of collision distance
        Fn = kn*f.δn - γn*δndt # Normal force

        # Tangential force based on tangential velocity, chord length and friction factor
        tmp = (-f.d + f.r1 - f.r2)*(-f.d - f.r1 + f.r2)*(-f.d + f.r1 + f.r2)*(f.d + f.r1 + f.r2)
        if tmp > 0
            chord_length = 1/f.d*sqrt( tmp )
        else
            chord_length = 2*min(r1,r2)
        end
        δtdt = -(f.ζ1*f.r1 + f.ζ2*f.r2) + vel2[2] + vel1[2] # Velocity difference in tangential direction
        Ft = mc.frict_fac*chord_length*δtdt # Tangential force. No static coefficient

        # Updating Trq and rotating forces back to update Fx and Fy
        Trq[idx1] = Ft*f.r1
        Fx[idx1],Fy[idx1] = Q'*[Fn Ft]'
    end

    @inline function getLowerbound(g::Vector{Float64},p::Float64)
        # Returns the index of the lowest grid point that is closest to the point p
        idx = findmin(abs.(g .- p))[2]
        if (p < g[idx]) && (p > g[1])
            idx -= 1
        end
        return idx
    end

    @inline function getUpperbound(g::Vector{Float64},p::Float64)
        # Returns the index of the highest grid point that is closest to the point p
        idx = findmin(abs.(g .- p))[2]
        if (p > g[idx]) && (p < g[end])
            idx += 1
        end
        return idx
    end

    @inline function CellAreaRatio(xg1::Float64, yg1::Float64, ref::Int64 , Δ::Float64, xfloe::Float64, yfloe::Float64, rfloe::Float64)
        # Δ is the refined spacing
        xc_r = range(xg1+Δ/2,step=Δ,length=ref)
        yc_r = transpose(range(yg1+Δ/2,step=Δ,length=ref))
        return sum( (xfloe .- xc_r).^2 .+ (yfloe .- yc_r).^2 .<= rfloe^2  ) / ref^2
    end

    @inline function FloeAreaRatio(xfloe::Float64,yfloe::Float64,rfloe::Float64,xg::Vector{Float64},yg::Vector{Float64})
        # Returns indeces of centered grid cells that encompass the floe along with the area fraction covered by ice

        Δx = xg[2] - xg[1] # Grid spacing

        # Indeces of the box that fully encompasses the circular floe
        xmin_idx = getLowerbound(xg, xfloe - rfloe)
        ymin_idx = getLowerbound(yg, yfloe - rfloe)
        xmax_idx = getUpperbound(xg, xfloe + rfloe)
        ymax_idx = getUpperbound(yg, yfloe + rfloe)

        # Index ranges for the box that encompasses the circle
        nx = xmax_idx - xmin_idx # number of centered cells in the x-direction
        ny = ymax_idx - ymin_idx # number of centered cells in the y-direction

        # Getting indeces and area ratios
        xidx = zeros(Int64,nx*ny)
        yidx = zeros(Int64,nx*ny)
        area_ratios = zeros(Float64,nx*ny)
        ref = 10 # grid refinement factor
        cnt = 1
        for j=1:ny
            for i=1:nx
                xidx[cnt] = i
                yidx[cnt] = j
                area_ratios[cnt] = CellAreaRatio(xg[xmin_idx]+(i-1)*Δx, yg[ymin_idx]+(j-1)*Δx, ref, Δx/ref, xfloe, yfloe, rfloe)
                cnt += 1
            end
        end
        xidx .+= xmin_idx .- 1
        yidx .+= ymin_idx .- 1
        idx = CartesianIndex.(Tuple.(eachrow(hcat(xidx,yidx)))) # indices on ocean grid corresponding to floe location
        return xidx,yidx,idx,area_ratios
    end

    @inline function deleteFloeID!(floes,ID)
        # Deletes a floe via its ID number (avoids loop confusion)
        ifloe = findfirst(isequal(ID),floes.ID)
        deleteat!(floes.x, ifloe)
        deleteat!(floes.y, ifloe)
        deleteat!(floes.θ, ifloe)
        deleteat!(floes.u, ifloe)
        deleteat!(floes.v, ifloe)
        deleteat!(floes.ζ, ifloe)
        deleteat!(floes.r, ifloe)
        deleteat!(floes.h, ifloe)
        deleteat!(floes.ID, ifloe)
        deleteat!(floes.cID, ifloe)
        deleteat!(floes.collision, ifloe)
    end

    @inline function appendFloe!(floes,x,y,θ,u,v,ζ,r,h,ID,cID,collision)
        # Appends a floe to floes
        push!(floes.x,x)
        push!(floes.y,y)
        push!(floes.θ,θ)
        push!(floes.u,u)
        push!(floes.v,v)
        push!(floes.ζ,ζ)
        push!(floes.r,r)
        push!(floes.h,h)
        push!(floes.ID,ID)
        push!(floes.cID,cID)
        push!(floes.collision,collision)
    end

    @inline function makeChannelPeriodic!(m::FloeModel)
        # Channel model configuration (periodic in x, walls in y)
        # Local referencing
        xleft = 0
        xright = m.grid.extent[1]
        floes = m.floes

        delete_IDs = [] # IDs of duplicate floes to delete
        for ifloe in 1:length(floes)
            x1 = floes.x[ifloe]
            r1 = floes.r[ifloe]
            xdr = xright - x1
            xdl = x1 - xleft
            ID = floes.ID[ifloe]
            cID = floes.cID[ifloe]
            collision = floes.collision[ifloe]

            # right side of domain
            if (xdr <= r1)

                # Create duplicate (enter right)
                if cID == 0
                    appendFloe!(floes,-xdr,floes.y[ifloe],floes.θ[ifloe],floes.u[ifloe],floes.v[ifloe],floes.ζ[ifloe],r1,floes.h[ifloe],-ID,ID,collision)
                    floes.cID[ifloe] = -ID

                # Original has exited domain from the right
                elseif (xdr < 0) & (-xdr < r1) & (ID > 0)
                    my_idx = findfirst(isequal(-ID),floes.ID) # index position of the corresponding duplicate
                    floes.ID[my_idx] = ID # duplicate becomes original
                    floes.cID[my_idx] = -ID # duplicate becomes original
                    floes.ID[ifloe] = -ID # original becomes duplicate
                    floes.cID[ifloe] = ID # original becomes duplicate

                # Delete duplicate (exit right)
                elseif (xdr < 0) & (-xdr > r1)
                    push!(delete_IDs,ID)
                    my_idx = findfirst(isequal(-ID),floes.ID) # index position of the corresponding original
                    floes.cID[my_idx] = 0
                end

            # left side of domain
            elseif xdl <= r1

                # Create duplicate (enter left)
                if cID == 0
                    appendFloe!(floes,xright+xdl,floes.y[ifloe],floes.θ[ifloe],floes.u[ifloe],floes.v[ifloe],floes.ζ[ifloe],r1,floes.h[ifloe],-ID,ID,collision)
                    floes.cID[ifloe] = -ID

                # Original has exited domain from the left
                elseif (xdl < 0) & (-xdl < r1) & (ID > 0)
                    my_idx = findfirst(isequal(-ID),floes.ID) # index position of the corresponding duplicate
                    floes.ID[my_idx] = ID # duplicate becomes original
                    floes.cID[my_idx] = -ID # duplicate becomes original
                    floes.ID[ifloe] = -ID # original becomes duplicate
                    floes.cID[ifloe] = ID # original becomes duplicate

                # Delete duplicate (exit left)
                elseif (xdl < 0) & (-xdl > r1)
                    push!(delete_IDs,ID)
                    my_idx = findfirst(isequal(-ID),floes.ID) # index position of the corresponding original
                    floes.cID[my_idx] = 0
                end
            end

        end

        # Delete duplicate floes that have exited
        for i=unique(delete_IDs)
            deleteFloeID!(floes,i)
        end
    end

    @inline function makeDoublyPeriodic!(m::FloeModel)
        # Doubly periodic model configuration (periodic in x, periodic in y)
        # Local referencing
        xleft = 0
        ybot = 0
        xright = m.grid.extent[1]
        ytop = m.grid.extent[2]
        floes = m.floes

        delete_IDs = [] # IDs of duplicate floes to delete
        for ifloe in 1:length(floes)
            x1 = floes.x[ifloe]
            y1 = floes.y[ifloe]
            r1 = floes.r[ifloe]
            xdr = xright - x1
            xdl = x1 - xleft
            ydt = ytop - y1
            ydb = y1 - ybot
            ID = floes.ID[ifloe]
            cID = floes.cID[ifloe]
            collision = floes.collision[ifloe]
            my_idx = findfirst(isequal(-ID),floes.ID) # index position of the corresponding mirrored floe

            # right side of domain
            if (xdr <= r1)

                # Create duplicate (enter right)
                if cID == 0
                    appendFloe!(floes,-xdr,floes.y[ifloe],floes.θ[ifloe],floes.u[ifloe],floes.v[ifloe],floes.ζ[ifloe],r1,floes.h[ifloe],-ID,ID,collision)
                    floes.cID[ifloe] = -ID

                # Exit right corner
                elseif (cID < 0) && (abs(floes.x[my_idx] - x1) < r1/1e3) && ((ydt <= r1) || (ydb <= r1))
                    floes.x[my_idx] = -xdr

                # Original has exited domain from the right
                elseif (xdr < 0) && (-xdr < r1) & (ID > 0)
                    my_idx = findfirst(isequal(-ID),floes.ID) # index position of the corresponding duplicate
                    floes.ID[my_idx] = ID # duplicate becomes original
                    floes.cID[my_idx] = -ID # duplicate becomes original
                    floes.ID[ifloe] = -ID # original becomes duplicate
                    floes.cID[ifloe] = ID # original becomes duplicate

                # Delete duplicate (exit right)
                elseif (xdr < 0) && (-xdr > r1)
                    push!(delete_IDs,ID)
                    my_idx = findfirst(isequal(-ID),floes.ID) # index position of the corresponding original
                    floes.cID[my_idx] = 0
                end

            # left side of domain
            elseif xdl <= r1
                # Create duplicate (enter left)
                if cID == 0
                    appendFloe!(floes,xright+xdl,floes.y[ifloe],floes.θ[ifloe],floes.u[ifloe],floes.v[ifloe],floes.ζ[ifloe],r1,floes.h[ifloe],-ID,ID,collision)
                    floes.cID[ifloe] = -ID

                #  Exit left corner
                elseif (cID < 0) && (abs(floes.x[my_idx] - x1) < r1/1e3) && ((ydb <= r1) || (ydt <= r1))
                    floes.x[my_idx] = xright+xdl

                # Original has exited domain from the left
                elseif (xdl < 0) && (-xdl < r1) && (ID > 0)
                    my_idx = findfirst(isequal(-ID),floes.ID) # index position of the corresponding duplicate
                    floes.ID[my_idx] = ID # duplicate becomes original
                    floes.cID[my_idx] = -ID # duplicate becomes original
                    floes.ID[ifloe] = -ID # original becomes duplicate
                    floes.cID[ifloe] = ID # original becomes duplicate

                # Delete duplicate (exit left)
                elseif (xdl < 0) && (-xdl > r1)
                    push!(delete_IDs,ID)
                    my_idx = findfirst(isequal(-ID),floes.ID) # index position of the corresponding original
                    floes.cID[my_idx] = 0
                end
            end

            # top side of domain
            if (ydt <= r1)
                # Enter top
                if cID == 0
                    appendFloe!(floes,floes.x[ifloe],-ydt,floes.θ[ifloe],floes.u[ifloe],floes.v[ifloe],floes.ζ[ifloe],r1,floes.h[ifloe],-ID,ID,collision)
                    floes.cID[ifloe] = -ID

                # Exit top corner
                elseif (cID < 0) && (abs(floes.y[my_idx] - y1) < r1/1e3) && ((xdr <= r1) || (xdl <= r1))
                    floes.y[my_idx] = -ydt

                # Original has exited domain from the top
                elseif (ydt < 0) && (-ydt < r1) && (ID > 0)
                    my_idx = findfirst(isequal(-ID),floes.ID) # index position of the corresponding duplicate
                    floes.ID[my_idx] = ID # duplicate becomes original
                    floes.cID[my_idx] = -ID # duplicate becomes original
                    floes.ID[ifloe] = -ID # original becomes duplicate
                    floes.cID[ifloe] = ID # original becomes duplicate

                # Exit top
                elseif (ydt < 0) && (-ydt > r1)
                    push!(delete_IDs,ID)
                    my_idx = findfirst(isequal(-ID),floes.ID) # index position of the corresponding original
                    floes.cID[my_idx] = 0
                end

            # bottom side of domain
            elseif ydb <= r1
                # Create duplicate (enter bottom)
                if cID == 0
                    appendFloe!(floes,floes.x[ifloe],ytop+ydb,floes.θ[ifloe],floes.u[ifloe],floes.v[ifloe],floes.ζ[ifloe],r1,floes.h[ifloe],-ID,ID,collision)
                    floes.cID[ifloe] = -ID

                #  Exit bottom corner
            elseif (cID < 0) && (abs(floes.y[my_idx] - y1) < r1/1e3) && ((xdl <= r1) || (xdr <= r1))
                    floes.y[my_idx] = ytop+ydb

                # Original has exited domain from the bottom
                elseif (ydb < 0) && (-ydb < r1) & (ID > 0)
                    my_idx = findfirst(isequal(-ID),floes.ID) # index position of the corresponding duplicate
                    floes.ID[my_idx] = ID # duplicate becomes original
                    floes.cID[my_idx] = -ID # duplicate becomes original
                    floes.ID[ifloe] = -ID # original becomes duplicate
                    floes.cID[ifloe] = ID # original becomes duplicate

                # Delete duplicate (exit bottom)
                elseif (ydb < 0) && (-ydb > r1)
                    push!(delete_IDs,ID)
                    my_idx = findfirst(isequal(-ID),floes.ID) # index position of the corresponding original
                    floes.cID[my_idx] = 0
                end
            end

        end

        # Delete duplicate floes that have exited
        for i=unique(delete_IDs)
            deleteFloeID!(floes,i)
        end

    end

    @inline function getColForces!(Fx::Vector{Float64},Fy::Vector{Float64},Trq::Vector{Float64},floes::StructArray{Floe},extent::Tuple{Float64,Float64},mc::NamedTuple,Δt::Float64)

        # Local referencing
        xleft = 0
        ybot = 0
        xright = extent[1]
        ytop = extent[2]

        nfloes = length(floes)
        for idx1 in 1:nfloes

            x1 = floes.x[idx1]
            y1 = floes.y[idx1]
            r1 = floes.r[idx1]
            h1 = floes.h[idx1]
            m1 = mc.rho_I*h1*π*r1^2
            u1 = floes.u[idx1]
            v1 = floes.v[idx1]
            ζ1 = floes.ζ[idx1]
            floes.collision[idx1] = false # Setting collision flag

            # Checking for collisions
            lst = collect(1:nfloes)
            deleteat!(lst,idx1)

            for idx2 in lst
                x2 = floes.x[idx2]
                y2 = floes.y[idx2]
                r2 = floes.r[idx2]
                h2 = floes.h[idx2]
                m2 = mc.rho_I*h2*π*r2^2
                δn = sqrt( (x2 - x1)^2 + (y2 - y1)^2 ) - (r1+r2) # collision distance
                if δn < 0 # collision occurs
                    u2 = floes.u[idx2]
                    v2 = floes.v[idx2]
                    ζ2 = floes.ζ[idx2]
                    theta = atan2(y2-y1,x2-x1)
                    d = sqrt( (x2 - x1)^2 + (y2 - y1)^2 ) # Distance between centers

                    # Updating Fx, Fy and Trq
                    collideTwoFloes!(Fx, Fy, Trq, idx1, (r1=r1, r2=r2 ,u1 =u1 ,v1=v1 ,u2=u2 ,v2=v2,ζ1=ζ1 ,ζ2=ζ2 ,h1=h1 ,h2=h2 ,d=d, theta=theta, δn=δn), mc)

                    # Setting collision flag
                    floes.collision[idx1] = true

                end
            end

            if mc.domain_config == "Channel" || mc.domain_config == "4Walls"
                # Check for collisions with top and bottom walls
                δn = -(y1 + r1 - ytop) # collision distance
                if δn < 0 # top wall
                    kn, kt, γn, γt = getCollisionFactors(mc.Ec,mc.Gc,mc.βc,mc.rho_I,r1,r1,h1,h1)
                    δndt = v1 # Rate of change of collision distance
                    δtdt = -ζ1*r1 + u1  # Rate of change of shear distance
                    Fy[idx1] += kn*δn - γn*δndt
                    Fx[idx1] += -γt*δtdt
                    Trq[idx1] += -Fx[idx1]*r1
                end

                δn = y1 - r1 - ybot # collision distance
                if δn < 0 # bottom wall
                    kn, kt, γn, γt = getCollisionFactors(mc.Ec,mc.Gc,mc.βc,mc.rho_I,r1,r1,h1,h1)
                    δndt = -v1 # Rate of change of collision distance
                    δtdt = ζ1*r1 + u1  # Rate of change of shear distance
                    Fy[idx1] += -kn*δn + γn*δndt
                    Fx[idx1] += -γt*δtdt
                    Trq[idx1] += Fx[idx1]*r1
                end
            end

            if mc.domain_config == "4Walls"
                # Collisions with left and right wals

                δn = -x1 - r1 + xright # collision distance
                if δn < 0 # right wall
                    kn, kt, γn, γt = getCollisionFactors(mc.Ec,mc.Gc,mc.βc,mc.rho_I,r1,r1,h1,h1)
                    δndt = u1 # Rate of change of collision distance
                    δtdt = ζ1*r1 + v1  # Rate of change of shear distance
                    Fx[idx1] += kn*δn - γn*δndt
                    Fy[idx1] += -γt*δtdt
                    Trq[idx1] += Fy[idx1]*r1
                end

                δn = x1 - r1 - xleft # collision distance
                if (x1 - xleft) <=  r1 # left wall
                    kn, kt, γn, γt = getCollisionFactors(mc.Ec,mc.Gc,mc.βc,mc.rho_I,r1,r1,h1,h1)
                    δndt = -u1 # Rate of change of collision distance
                    δtdt = -ζ1*r1 + v1  # Rate of change of shear distance
                    Fx[idx1] += -kn*δn + γn*δndt
                    Fy[idx1] += -γt*δtdt
                    Trq[idx1] += -Fy[idx1]*r1
                end
            end
        end
    end


    @inline function getOceanCoupling!(m::FloeModel,Fx::Vector{Float64},Fy::Vector{Float64},Trq::Vector{Float64},dVol::Vector{Float64})
        # Local referencing
        floes = m.floes
        nfloes = length(floes)
        mc = m.mc
        Δx = m.grid.extent[1] / m.grid.size[1]

        # Clearing SIfract
        m.SIfract = zeros(Float64,m.grid.size[1], m.grid.size[2])

        for ifloe in 1:nfloes
            xidx_floe,yidx_floe,idx_floe,area_ratios = FloeAreaRatio(floes.x[ifloe],floes.y[ifloe],floes.r[ifloe],m.grid.xg,m.grid.yg)

            if length(area_ratios) > 0
                # Ice stress on ocean
                lx = m.grid.xc[xidx_floe] .- floes.x[ifloe] # x-moment arm [m]
                ly = m.grid.yc[yidx_floe] .- floes.y[ifloe] # y-moment arm [m]
                ui = floes.u[ifloe] .- ly.*floes.ζ[ifloe] # ice velocity vector
                vi = floes.v[ifloe] .+ lx.*floes.ζ[ifloe] # ice velocity vector
                taux_IO = mc.rho_O*mc.Cd_IO*(ui .- m.u_O[idx_floe]).*abs.(ui .- m.u_O[idx_floe])
                tauy_IO = mc.rho_O*mc.Cd_IO*(vi .- m.v_O[idx_floe]).*abs.(vi .- m.v_O[idx_floe])
                m.taux[idx_floe] .= taux_IO.*area_ratios .+ m.taux[idx_floe].*(1 .- area_ratios)
                m.tauy[idx_floe] .= tauy_IO.*area_ratios .+ m.tauy[idx_floe].*(1 .- area_ratios)

                # Forces and torque on floes
                Fxx = (-taux_IO + mc.rho_A*mc.Cd_IA*(m.u_A[idx_floe] .- ui).*abs.(m.u_A[idx_floe] .- ui)).*area_ratios*(Δx^2) # Force at each grid point
                Fyy = (-tauy_IO + mc.rho_A*mc.Cd_IA*(m.v_A[idx_floe] .- vi).*abs.(m.v_A[idx_floe] .- vi)).*area_ratios*(Δx^2) # Force at each grid point
                Fx[ifloe] = sum(Fxx)
                Fy[ifloe] = sum(Fyy)
                Trq[ifloe] = sum(lx.*Fyy .- ly.*Fxx)

                # Ice/ocean heat flux on ocean
                hlfx_IO = mc.q_IO*(mc.Tf .- m.T_O[idx_floe])
                m.hflx[idx_floe] .= hlfx_IO.*area_ratios + m.hflx[idx_floe].*(1 .- area_ratios)

                # Floe volume change due to heat flux and melt rate imposed on ocean & associated melt rate
                dh = m.Δt/(mc.rho_I*mc.Lf)*( mc.Q*(1 - mc.alpha_I) - (mc.A + mc.B*mc.Tf) .- hlfx_IO).*area_ratios
                m.melt_rate[idx_floe] = dh./m.Δt
                dVol[ifloe] = -sum(dh*Δx^2)

                # Mask update
                m.SIfract[idx_floe] .+= area_ratios
            end
        end
    end
    @inline function dudt(m::FloeModel,x::Vector{Float64},y::Vector{Float64},Uf::Vector{Float64},Vf::Vector{Float64},ζ::Vector{Float64},Kx::Vector{Float64},Ky::Vector{Float64},Ktrq::Vector{Float64},dVol::Vector{Float64},stepi::Float64)
        floes = m.floes
        nfloes = length(floes)
        mc = m.mc
        Δx = m.grid.extent[1] / m.grid.size[1]
        taux=copy(m.taux)
        tauy=copy(m.tauy)
        # Clearing SIfract
        m.SIfract = zeros(Float64,m.grid.size[1], m.grid.size[2])

        for ifloe in 1:nfloes
            xidx_floe,yidx_floe,idx_floe,area_ratios = FloeAreaRatio(floes.x[ifloe],floes.y[ifloe],floes.r[ifloe],m.grid.xg,m.grid.yg)    #######what is this 

            if length(area_ratios) > 0
                # Ice stress on ocean
                lx = m.grid.xc[xidx_floe] .- x[ifloe] # x-moment arm [m]
                ly = m.grid.yc[yidx_floe] .- y[ifloe] # y-moment arm [m]
                ui = Uf[ifloe] .- ly.*ζ[ifloe] # ice velocity vector
                vi = Vf[ifloe] .+ lx.*ζ[ifloe] # ice velocity vector
                taux_IO = mc.rho_O*mc.Cd_IO*(ui .- m.u_O[idx_floe]).*abs.(ui .- m.u_O[idx_floe])
                tauy_IO = mc.rho_O*mc.Cd_IO*(vi .- m.v_O[idx_floe]).*abs.(vi .- m.v_O[idx_floe])
                if stepi==2
                    m.taux[idx_floe] .= taux_IO.*area_ratios .+ m.taux[idx_floe].*(1 .- area_ratios)
                    m.tauy[idx_floe] .= tauy_IO.*area_ratios .+ m.tauy[idx_floe].*(1 .- area_ratios)
                end
                taux[idx_floe] .= taux_IO.*area_ratios .+ taux[idx_floe].*(1 .- area_ratios)
                tauy[idx_floe] .= tauy_IO.*area_ratios .+ tauy[idx_floe].*(1 .- area_ratios)

                # Forces and torque on floes
                Fxx = (-taux_IO + mc.rho_A*mc.Cd_IA*(m.u_A[idx_floe] .- ui).*abs.(m.u_A[idx_floe] .- ui)).*area_ratios*(Δx^2) # Force at each grid point
                Fyy = (-tauy_IO + mc.rho_A*mc.Cd_IA*(m.v_A[idx_floe] .- vi).*abs.(m.v_A[idx_floe] .- vi)).*area_ratios*(Δx^2) # Force at each grid point
                Kx[ifloe] = sum(Fxx)
                Ky[ifloe] = sum(Fyy)
                Ktrq[ifloe] = sum(lx.*Fyy .- ly.*Fxx)

                # Ice/ocean heat flux on ocean
                hlfx_IO = mc.q_IO*(mc.Tf .- m.T_O[idx_floe])
                m.hflx[idx_floe] .= hlfx_IO.*area_ratios + m.hflx[idx_floe].*(1 .- area_ratios)

                # Floe volume change due to heat flux and melt rate imposed on ocean & associated melt rate
                dh = m.Δt/(mc.rho_I*mc.Lf)*( mc.Q*(1 - mc.alpha_I) - (mc.A + mc.B*mc.Tf) .- hlfx_IO).*area_ratios
                m.melt_rate[idx_floe] = dh./m.Δt
                dVol[ifloe] = -sum(dh*Δx^2)

                # Mask update
                m.SIfract[idx_floe] .+= area_ratios
            end
        end
    end
    @inline function updateFloes!(m::FloeModel,Fx::Vector{Float64},Fy::Vector{Float64},Trq::Vector{Float64},dVol::Vector{Float64})

        # Local referencing
        floes = m.floes
        nfloes = length(floes)
        mc = m.mc
        Δt = m.Δt

        # Acccounting for duplicates in the forces/torque on floes
        lst = findall(var->var<0, floes.cID) # Indeces of all originals that have duplicates
        for ifloe in lst
            ifloeC = findfirst(isequal(-floes.ID[ifloe]),floes.ID) # index position of the corresponding duplicate
            Fx[ifloe] += Fx[ifloeC]
            Fy[ifloe] += Fy[ifloeC]
            Trq[ifloe] += Trq[ifloeC]
            dVol[ifloe] += dVol[ifloeC]
            Fx[ifloeC] = copy(Fx[ifloe])
            Fy[ifloeC] = copy(Fy[ifloe])
            Trq[ifloeC] = copy(Trq[ifloe])
            dVol[ifloeC] = copy(dVol[ifloe])
            # Updating the collision Boolean
            floes.collision[ifloe] = floes.collision[ifloe] || floes.collision[ifloeC] # Collision True if either original or duplication experiences a collision (or both)
        end

        # Updating floe velocities, rotation rate
        for ifloe=1:nfloes
            mass = (mc.rho_I*floes.h[ifloe]*π*floes.r[ifloe]^2)
            floes.u[ifloe] += Fx[ifloe]/mass*Δt
            floes.v[ifloe] += Fy[ifloe]/mass*Δt
            floes.ζ[ifloe] += Trq[ifloe]*Δt/(0.5*mass*floes.r[ifloe]^2)
        end

        # Updating floe positions, angle and thickness
        for ifloe=1:nfloes
            floes.x[ifloe] += floes.u[ifloe]*Δt
            floes.y[ifloe] += floes.v[ifloe]*Δt
            floes.θ[ifloe] += 2*π*floes.ζ[ifloe]*Δt
            floes.h[ifloe] += dVol[ifloe]/(π*floes.r[ifloe]^2)
            if floes.θ[ifloe]>2*π
                floes.θ[ifloe] -= 2*π
            end
        end

        # Deleting if floe has melted
        delete_IDs = [] # IDs of floes to delete
        for ifloe=1:nfloes
            if floes.h[ifloe] < 0.05
                push!(delete_IDs,floes.ID[ifloe])
            end
        end

        for i=unique(delete_IDs)
            deleteFloeID!(floes,i)
        end

    end

    @inline function timestepFloes!(m::FloeModel)
        floes = m.floes
        nfloes = length(floes)
        Δt = m.Δt
        mc = m.mc
        # Applying periodicity
        if m.mc.domain_config == "Channel"
            makeChannelPeriodic!(m)
        elseif m.mc.domain_config == "DoublyPeriodic"
            makeDoublyPeriodic!(m)
        end

        # Getting coupled forces and fluxes
        Fx = zeros(Float64,length(m.floes))
        Fy = zeros(Float64,length(m.floes))
        Trq = zeros(Float64,length(m.floes))
        Kx = zeros(Float64,length(m.floes))
        Ky = zeros(Float64,length(m.floes))
        Ktrq = zeros(Float64,length(m.floes))
        Kx1 = zeros(Float64,length(m.floes))
        Ky1 = zeros(Float64,length(m.floes))
        Ktrq1 = zeros(Float64,length(m.floes))
        Kx2 = zeros(Float64,length(m.floes))
        Ky2 = zeros(Float64,length(m.floes))
        Ktrq2 = zeros(Float64,length(m.floes))
        Kx3 = zeros(Float64,length(m.floes))
        Ky3 = zeros(Float64,length(m.floes))
        Ktrq3 = zeros(Float64,length(m.floes))
        dVol = zeros(Float64,length(m.floes))
        Uf1=copy(floes.u) 
        Vf1=copy(floes.v) 
        x1=copy(floes.x) 
        y1=copy(floes.y) 
        ζ1=copy(floes.ζ)

        Uf2=copy(floes.u) 
        Vf2=copy(floes.v) 
        x2=copy(floes.x) 
        y2=copy(floes.y) 
        ζ2=copy(floes.ζ)
        
        Uf3=copy(floes.u) 
        Vf3=copy(floes.v) 
        x3=copy(floes.x) 
        y3=copy(floes.y) 
        ζ3=copy(floes.ζ)

        # getOceanCoupling!(m,Fx,Fy,Trq,dVol);
        dudt(m,floes.x,floes.y,floes.u,floes.v,floes.ζ,Kx,Ky,Ktrq,dVol,1.0)
        # getColForces!(Kx,Ky,Ktrq,m.floes,m.grid.extent,mc,Δt/2)
        nfloes = length(floes) 
        for ifloe=1:nfloes
            mass = (mc.rho_I*floes.h[ifloe]*π*floes.r[ifloe]^2)
        # mass = (mc.rho_I*floes.h.*floes.r.^2*π)
    
            Uf1[ifloe]+= Kx[ifloe]/mass*Δt/2;
            Vf1[ifloe]+= Ky[ifloe]/mass*Δt/2;
            x1[ifloe]+= Uf1[ifloe]*Δt/2;
            y1[ifloe]+=Vf1[ifloe]*Δt/2;
            ζ1[ifloe]+= Ktrq[ifloe]*Δt*0.5/(1*mass*floes.r[ifloe]^2)
        end
        dudt(m,x1,y1,Uf1,Vf1,ζ1,Kx1,Ky1,Ktrq1,dVol,1.0)
        # getColForces!(Kx1,Ky1,Ktrq1,m.floes,m.grid.extent,mc,Δt/2)
        nfloes = length(floes)
        for ifloe=1:nfloes
            mass = (mc.rho_I*floes.h[ifloe]*π*floes.r[ifloe]^2)
        # mass = (mc.rho_I*floes.h.*floes.r.^2*π)
    
            Uf2[ifloe]+= Kx1[ifloe]/mass*Δt/2;
            Vf2[ifloe]+= Ky1[ifloe]/mass*Δt/2;
            x2[ifloe]+= Uf2[ifloe]*Δt/2;
            y2[ifloe]+=Vf2[ifloe]*Δt/2;
            ζ2[ifloe]+= Ktrq1[ifloe]*Δt*0.5/(1*mass*floes.r[ifloe]^2)
        end
        dudt(m,x2,y2,Uf2,Vf2,ζ2,Kx2,Ky2,Ktrq2,dVol,1.0)
        # getColForces!(Kx2,Ky2,Ktrq2,m.floes,m.grid.extent,mc,Δt/2)
        nfloes = length(floes)
        for ifloe=1:nfloes
            mass = (mc.rho_I*floes.h[ifloe]*π*floes.r[ifloe]^2)
        # mass = (mc.rho_I*floes.h.*floes.r.^2*π)
    
            Uf3[ifloe]+= Kx2[ifloe]/mass*Δt;
            Vf3[ifloe]+= Ky2[ifloe]/mass*Δt;
            x3[ifloe]+= Uf3[ifloe]*Δt;
            y3[ifloe]+=Vf3[ifloe]*Δt;
            ζ3[ifloe]+= Ktrq2[ifloe]*Δt/(1*mass*floes.r[ifloe]^2)
        end
        dudt(m,x3,y3,Uf3,Vf3,ζ3,Kx3,Ky3,Ktrq3,dVol,1.0)
        # getColForces!(Kx3,Ky3,Ktrq3,m.floes,m.grid.extent,mc,Δt)
        
        Fx=Kx./6.0+Kx1./3.0+Kx2./3.0+Kx3./6.0
        Fy=Ky./6.0+Ky1./3.0+Ky2./3.0+Ky3./6.0
        Trq=Ktrq./6.0+Ktrq1./3.0+Ktrq2./3.0+Ktrq3./6.0
        dudt(m,floes.x,floes.y,floes.u,floes.v,floes.ζ,Kx,Ky,Ktrq,dVol,2.0)
        getColForces!(Fx,Fy,Trq,m.floes,m.grid.extent,m.mc,m.Δt)
        updateFloes!(m,Fx,Fy,Trq,dVol)
        print("[ ",Fx[2]," ]  [ ",Fy[2]," ] ")

    end

    @inline function timestepAtm!(m::FloeModel)
        mc = m.mc
        m.taux .= mc.rho_A*mc.Cd_AO*(m.u_A .- m.u_O).*abs.(m.u_A .- m.u_O)
        m.tauy .= mc.rho_A*mc.Cd_AO*(m.v_A .- m.v_O).*abs.(m.v_A .- m.v_O)
        m.hflx .= mc.Q*(1 - mc.alpha_O) .- (mc.A .+ mc.B*m.T_O )
    end

end # end module
