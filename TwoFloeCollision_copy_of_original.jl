using StructArrays
using BenchmarkTools
using Printf
using NetCDF
using Plots, Statistics

# include("FloeModule_current.jl")
include("FloeModule.jl")
# include("FloeModule - RK4.jl")
include("C:/Users/shiva/Downloads/distribute_rand_floes_shivam.jl")

############################################################################

function plot_floes(xfloes,yfloes,rfloes,hfloes,theta_floes,zeta_floes,Lx,Ly,i)
    nfloes = length(rfloes)
    for ifloe in 1:nfloes
        if hfloes[i,ifloe] != 0.0
            fillalpha = 0.8*hfloes[i,ifloe]/max_floe_h
            plot!(FloeModule.FloeShape(xfloes[i,ifloe]/1000,yfloes[i,ifloe]/1000,rfloes[ifloe]/1000),seriestype = [:shape,], lw = 0.5, c=:blue,linecolor=:black,fillalpha=fillalpha,xlims=(0,Lx/1000),ylims=(0,Ly/1000),aspect_ratio=1,xlabel="x [km]",ylabel="y [km]")
            clr = :red
            if zeta_floes[i,ifloe]<0
                clr=:green
            end
            plot!(FloeModule.FloeOrient(xfloes[i,ifloe]/1000,yfloes[i,ifloe]/1000,rfloes[ifloe]/1000,theta_floes[i,ifloe]),seriestype=[:shape,],lw = 6,linecolor=clr,fillalpha=fillalpha,xlims=(0,Lx/1000),ylims=(0,Ly/1000),aspect_ratio=1,xlabel="x [km]",ylabel="y [km]")
        end
    end
end

function random_distribution(conc,Lx,Ly,rmin,rmax,hmin,hmax)
    # Generates random normal distribution of floes with concentration close to: conc
    eps = sqrt(Lx^2 + Ly^2)/100 # minimum distance between floes and between boundaries and floes
    h = hmin + (hmax - hmin)*rand()
    r = rmin + (rmax - rmin)*rand()
    all_h = [h]
    all_r = [r]
    all_x = [r + eps + (Lx -2*eps -2*r)*rand()]
    all_y = [r + eps + (Ly -2*eps -2*r)*rand()]
    current_conc = 0
    tot_cnt = 0
    max_cnt = 1e5
    while current_conc < conc && tot_cnt < max_cnt
        h = hmin + (hmax - hmin)*rand()
        r = rmin + (rmax - rmin)*rand()
        inters = true
        cnt = 0
        # Finding the x,y coordinates
        while inters == true && tot_cnt < max_cnt
            if cnt > 100
                rmax = max(rmax/2,rmin)
                cnt = 0
            end
            x = r + eps + (Lx -2*eps -2*r)*rand()
            y = r + eps + (Ly -2*eps -2*r)*rand()
            all_dist = sqrt.((x .- all_x).^2 + (y .- all_y).^2)
            all_thresh = eps + r .+ all_r
            if !any(x->x<0, all_dist - all_thresh)
                inters = false
                all_x = [all_x x]
                all_y = [all_y y]
                all_h = [all_h h]
                all_r = [all_r r]
            end
            cnt += 1
            tot_cnt += 1
        end
        # Updating concentration
        current_conc = sum(π*(all_r).^2) / (Lx*Ly)
    end
    # Generating floes
    nfloes = length(all_r)
    floes = StructArray(FloeModule.Floe(ID=i) for i in 1:nfloes)
    for i=1:nfloes
        floes.x[i] = all_x[i]
        floes.y[i] = all_y[i]
        floes.r[i] = all_r[i]
        floes.h[i] = all_h[i]
        floes.ID[i] = i
    end
    (tot_cnt >= max_cnt) && @printf("Exceeed maximum tries to generate floes: count = %d \n", tot_cnt);
    @printf("Generated floes with concentration: %.3f \n",current_conc)
    return floes
end
er_all = 0.1:0.1:0.9;
nsim = length(er_all);
# delta_KE_all=Array{Float64,1}(undef,nsim);
############################################################################
#for irti=1:nsim
## All model constants
# Ocean
α = 4.9467e-05           # Thermal expansion coefficient (ref. at 0C & 34 psu)
β = 7.8137e-04           # Haline contraction coefficient (ref. at 0C & 34 psu)
Ω = 2*π/(3600*24)
f = 2*Ω*sind(70) # Coriolis parameter
Tref = 0.
Sref = 34.
rho_O = 1027.1719 # Mean ocean density [kg/m3]. Calcuated at Tref and Sref
cp_O = 3991. # Sensible heat capacity of ocean (J/kg/K)
alpha_O = 0.4 # Ocean albedo
Tf = -1.8 # Ocean freezing point ( C)
# Atmosphere
Q = 300. # Solar constant for Arctic summer (W/m2)
A = 90. # Upward flux constant A (W/m2)
B = 15. # Upward flux constant B (W/m2/K)
rho_A = 1.25 # Atmosphere density [kg/m3]
Cd_AO = 0.00125 # Atm/ocean momentum drag coefficient. Default: 0.00125 #Change it to 0 for plot of  change in kinetic energy 
# Ice floes
E = 5.0*1e5 # Young's modulus (elastic modulus). Default: 6.0*1e9 (Pa)
er = 0.5 # coefficient of restitution
ν = 0.33 # Poisson's ratio
Ec,Gc,βc = FloeModule.prepareCollisionFactors(E,ν,er)
# Cd_IO =0.0055 # Ice/ocean momentum drag coefficient. Default: 0.0055 #Change it to 0 for plot of  change in kinetic energy 
Cd_IO =0.0
# Cd_IA =0.00125 # Ice/atm momentum drag coefficient # TODO: Find appropriate value for this. I'm using 0.00125 #Change it to 0 for plot of  change in kinetic energy 
Cd_IA =0.0
rho_I = 1000. # Ice floe density [kg/m3]
Lf = 1000*334.  # Latent heat of freezing (J/kg)
alpha_I = 0.7 # Ice albedo
q_IO = 100. # Vertical mixing coeff at ice base (W/m2/K). Default:100
frict_fac = 1000. # Friction factor [Ns/m^2] . Default:100
domain_config = "Channel" # Choose between "4Walls","Channel","DoublyPeriodic"
mc = (f=f,rho_O=rho_O,cp_O=cp_O,alpha_O=alpha_O,Tf=Tf,
      Q=Q,A=A,B=B,rho_A=rho_A,Cd_AO=Cd_AO,
      Ec=Ec,Gc=Gc,βc=βc,Cd_IO=Cd_IO,Cd_IA=Cd_IA,rho_I=rho_I,Lf=Lf,alpha_I=alpha_I,q_IO=q_IO,frict_fac=frict_fac,domain_config=domain_config) # All model constants

## Creating floe model
Nx = 512
Ny = 512
Δy = 250.0
Δx = 250.0
Lx = Δx*Nx
Ly = Δy*Ny
floe_grid = FloeModule.FloeGrid(size=(Nx, Ny), extent=(Δx*Nx, Δx*Nx))
xc = floe_grid.xc
yc = floe_grid.yc

# Timestep
Δt = 60.
n_iter = 20000
n_iter_save = 200
nsave = n_iter_save+1
it_interval = Int(n_iter/n_iter_save)
time_interval = Δt*it_interval
tsave = collect(time_interval*range(1,stop=nsave,length=nsave))

# Floes
# nfloes_orig = 2
# floes = StructArray(FloeModule.Floe(ID=i) for i in 1:nfloes_orig)
# floes.r[1] = 7000
# floes.x[1] = Lx/4 
# floes.y[1] = Ly/2
# floes.h[1] = 1;
# floes.u[1] = itt;

# floes.r[2] = 7000
# floes.x[2] = 3*Lx/4 
# floes.y[2] = Ly/2
# floes.h[2] = 1;
# floes.u[2] = -itt;

# floes.r[3] = 5000
# floes.x[3] = Lx/2
# floes.y[3] = 3*Ly/4
# floes.h[3] = 1;
# floes.v[3] = -0.5;

# floes.r[4] = 5000
# floes.x[4] = Lx/2 
# floes.y[4] = Ly/4
# floes.h[4] = 1;
# floes.v[4] = 0.5;

# floes.r[5] = 5000
# floes.x[5] = Lx/2 
# floes.y[5] = Ly/2
# floes.h[5] = 1;
# floes.u[5] = 0;

rfloes_bin = [2, 3, 5, 10]*1e3
nfloes_bin = [100, 30, 10, 2]
hfloes_bin = [1. 1. 1. 1.]
conc = π*sum(nfloes_bin.*rfloes_bin.^2) / (Lx*Ly)
@printf("concentration is: %f \n",conc)
floes = distributeFloes(rfloes_bin,nfloes_bin,hfloes_bin,Lx,Ly)
nfloes_orig = length(floes)

# conc = 0.4
# rmax = 10000
# rmin = 700
# hmin = 2.5
# hmax = 2.5
umag = 0.1
# floes = random_distribution(conc,Lx,Ly,rmin,rmax,hmin,hmax)
# nfloes_orig = length(floes.x)
floes.u .= 2*umag*(rand(nfloes_orig) .- 0.5)
floes.v .= 2*umag*(rand(nfloes_orig) .- 0.5)


# mass_1=(mc.rho_I*floes.h[1]*π*floes.r[1]^2);
# mass_2=(mc.rho_I*floes.h[2]*π*floes.r[2]^2);
# KE_1=1/2*mass_1*floes.u[1]*floes.u[1]+1/2*mass_2*floes.u[2]*floes.u[2];

# Making floe model
floe_model = FloeModule.FloeModel(mc=mc,floes=floes,grid=floe_grid,Δt=Δt)

## Writer for surface boundary conditions
fprefix = split(split(@__FILE__, string(pwd(), "\\"))[2], ".")[1]
if ~isdir("data/"*fprefix)
mkdir("data/"*fprefix)
end
# fname_bcs = string("data/",fprefix,"/surface_bcs.nc")
# if isfile(fname_bcs)
# rm(fname_bcs)
# end
# nccreate(fname_bcs,"taux","xc",xc,"yc",yc,"t",tsave,atts=Dict("units"=>"N/m2"))
# nccreate(fname_bcs,"tauy","xc",xc,"yc",yc,"t",tsave,atts=Dict("units"=>"N/m2"))
# nccreate(fname_bcs,"SIfract","xc",xc,"yc",yc,"t",tsave,atts=Dict("units"=>""))
# nccreate(fname_bcs,"hflx","xc",xc,"yc",yc,"t",tsave,atts=Dict("units"=>"W/m2"))
# nccreate(fname_bcs,"melt_rate","xc",xc,"yc",yc,"t",tsave,atts=Dict("units"=>"m/s"))
# taux_all = zeros(Nx,Ny,nsave)
# tauy_all = zeros(Nx,Ny,nsave)
# SIfract_all = zeros(Nx,Ny,nsave)
# hflx_all = zeros(Nx,Ny,nsave)
# melt_rate_all = zeros(Nx,Ny,nsave)
# Writer for floes
fname_floes = string("data/",fprefix,"/floes.nc")
if isfile(fname_floes)
rm(fname_floes)
end
nccreate(fname_floes,"x","t",tsave,"floe_ID",nfloes_orig,atts=Dict("units"=>"m"))
nccreate(fname_floes,"y","t",tsave,"floe_ID",nfloes_orig,atts=Dict("units"=>"m"))
nccreate(fname_floes,"u","t",tsave,"floe_ID",nfloes_orig,atts=Dict("units"=>"m/s"))
nccreate(fname_floes,"v","t",tsave,"floe_ID",nfloes_orig,atts=Dict("units"=>"m/s"))
nccreate(fname_floes,"theta","t",tsave,"floe_ID",nfloes_orig,atts=Dict("units"=>"rad"))
nccreate(fname_floes,"zeta","t",tsave,"floe_ID",nfloes_orig,atts=Dict("units"=>"1/s"))
nccreate(fname_floes,"h","t",tsave,"floe_ID",nfloes_orig,atts=Dict("units"=>"m"))
nccreate(fname_floes,"r","floe_ID",nfloes_orig,atts=Dict("units"=>"m"))
nccreate(fname_floes,"collision_time","t",tsave,"floe_ID",nfloes_orig,atts=Dict("units"=>"s"))
xfloes_all = zeros(nsave,nfloes_orig)
yfloes_all = zeros(nsave,nfloes_orig)
ufloes_all = zeros(nsave,nfloes_orig)
vfloes_all = zeros(nsave,nfloes_orig)
θfloes_all = zeros(nsave,nfloes_orig)
ζfloes_all = zeros(nsave,nfloes_orig)
hfloes_all = zeros(nsave,nfloes_orig)
rfloes_all = copy(floes.r)
collision_time_all = zeros(nsave,nfloes_orig)
for ifloe in 1:nfloes_orig
    xfloes_all[1,ifloe] = floes.x[ifloe]
    yfloes_all[1,ifloe] = floes.y[ifloe]
    ufloes_all[1,ifloe] = floes.u[ifloe]
    vfloes_all[1,ifloe] = floes.v[ifloe]
    θfloes_all[1,ifloe] = floes.θ[ifloe]
    ζfloes_all[1,ifloe] = floes.ζ[ifloe]
    hfloes_all[1,ifloe] = floes.h[ifloe]
    rfloes_all[ifloe] = floes.r[ifloe]
end

#Initializing ocean field
mid_idx = Int(Nx/2)
floe_model.u_O[1:mid_idx,:] .=1	#change it to 0 for plots
floe_model.u_O[mid_idx+1:end,:] .= -1	#change it to 0 for plots we are working with ocean at rest 

## Running
cnt_save = 2
t_tot = 0
for i=1:n_iter
    global t_tot, cnt_save
    t_tot += Δt

    ## Running floe model
    nfloes = length(floes)
    if nfloes>0
        FloeModule.timestepFloes!(floe_model)
        nfloes = length(floes)
    end

    ## Updating collision freqency (do every time step)
    for ID in 1:nfloes_orig
        ifloe = findfirst(isequal(ID),floes.ID)
        if ~isnothing(ifloe)
            collision_time_all[cnt_save,ID] += Δt*floe_model.floes.collision[ifloe]
        end
    end

    # Print iteration
    println(i)

    # Printing and saving
    if i%it_interval == 0
        # Printing a progress message
        flush(stdout)
        # Floe output
        for ID in 1:nfloes_orig
            ifloe = findfirst(isequal(ID),floes.ID)
            if ~isnothing(ifloe)
                xfloes_all[cnt_save,ID] = floe_model.floes.x[ifloe]
                yfloes_all[cnt_save,ID] = floe_model.floes.y[ifloe]
                ufloes_all[cnt_save,ID] = floe_model.floes.u[ifloe]
                vfloes_all[cnt_save,ID] = floe_model.floes.v[ifloe]
                θfloes_all[cnt_save,ID] = floe_model.floes.θ[ifloe]
                ζfloes_all[cnt_save,ID] = floe_model.floes.ζ[ifloe]
                hfloes_all[cnt_save,ID] = floe_model.floes.h[ifloe]
            end
        end
        # # # Surface BC output
        # taux_all[:,:,cnt_save] = floe_model.taux
        # tauy_all[:,:,cnt_save] = floe_model.tauy
        # SIfract_all[:,:,cnt_save] = floe_model.SIfract
        # hflx_all[:,:,cnt_save] = floe_model.hflx
        # melt_rate_all[:,:,cnt_save] = floe_model.melt_rate

        cnt_save += 1
    end

end
# KE_2=1/2*mass_1*floes.u[1]*floes.u[1]+1/2*mass_2*floes.u[2]*floes.u[2];
# delta_KE_all[1]=KE_2-KE_1;

## Writing output
# Floe
ncwrite(xfloes_all,fname_floes,"x")
ncwrite(yfloes_all,fname_floes,"y")
ncwrite(ufloes_all,fname_floes,"u")
ncwrite(vfloes_all,fname_floes,"v")
ncwrite(θfloes_all,fname_floes,"theta")
ncwrite(ζfloes_all,fname_floes,"zeta")
ncwrite(rfloes_all,fname_floes,"r")
ncwrite(hfloes_all,fname_floes,"h")
ncwrite(collision_time_all,fname_floes,"collision_time")
# # Surface BCs
# ncwrite(taux_all,fname_bcs,"taux")
# ncwrite(tauy_all,fname_bcs,"tauy")
# ncwrite(SIfract_all,fname_bcs,"SIfract")
# ncwrite(hflx_all,fname_bcs,"hflx")
# ncwrite(melt_rate_all,fname_bcs,"melt_rate")

# ################################################################################
## Making movie

# Floe data
xfloes = xfloes_all
yfloes = yfloes_all
ufloes = ufloes_all
vfloes = vfloes_all
θfloes = θfloes_all
ζfloes = ζfloes_all
hfloes = hfloes_all
rfloes = rfloes_all
theta_floes = θfloes_all
zeta_floes = ζfloes_all
max_floe_h = maximum(hfloes)
nfloes = length(rfloes)

# # BCs data
# SIfract = ncread(fname_bcs,"SIfract")

title_plt = plot(title=string(t_tot), grid = false, showaxis = false, bottom_margin = -30Plots.px)
t_tot = 0
anim = @animate for i=1:nsave
    global t_tot
    t_tot += Δt*n_iter/n_iter_save

    floes_plt = plot(xlabel="x [km]", ylabel="y [km]",title=string(" Euler"),legend=false)
    plot_floes(xfloes,yfloes,rfloes,hfloes,theta_floes,zeta_floes,Lx,Ly,i)

    plot(floes_plt,size=(500, 400))

end
mp4(anim,string("videos/"," RK4 ",fprefix,itt,".mp4"), fps = 15) # hide

