
# Monte Carlo for Annual Review; 
# Quantitative Spatial Model;
# Helpman (1998) model;
# Countries and regions;

# July, 2016;

# *********************; 
# **** Choose User ****; 
# *********************;



cd("/Users/Pol/Dropbox/Mac/Desktop/Scolarite/Master/M1 - économie du développement/S2/Urban economics/qsm/ARE_MatlabPrograms/replication_file") 

# ************************; 
# **** Initialization ****;
# ************************;

# Set default random number stream;
using Random, LinearAlgebra, Statistics

s = Random.MersenneTwister(1)
Random.seed!(s)


# *************************; 
# **** Distance matrix ****;
# *************************;

# Locations on a N * N latitude and longitude grid
# Implies N * N locations
# Distance matrix is N * N 

N=30
NN=N*N
# Other latitude-longitude grid
ltd = range(0, 4, length=N)'
lgd = range(0, 4, length=N)


# Transport costs for each point on grid;
tt=1;
tau=zeros(N);
tau .=tt; #### créer un vecteur de 1 ????
tau = ones(N);

# Compute weighted distance using these transport costs;
## Distances
# Start by creating a matrix of coordinates.
coords = zeros(NN,2)
coords[:,1] = floor.(collect((0:(NN-1))/N))
coords[:,2] = (collect(0:(NN-1)) - (coords[:,1])*N)

# Then create the matrix of distances (quelques doutes sur sa création)
dx = coords[:,1] .- coords[:,1]'
dy = coords[:,2] .- coords[:,2]'
dist  = real(sqrt.((coords[:,1] .- coords[:,1]').^2 + (coords[:,2] .- coords[:,2]').^2))

# Own iceberg transport costs are one
for i in 1:length(diag(dist))
    dist[i,i] = 1
end

# Trade costs are a power function of effective distance;
dist=dist.^0.33;

# Define east and west as two countries;
Iwest = zeros(N,N)
Ieast = zeros(N,N)
Iwest[:, 1:div(N, 2)] .= 1
Ieast[:, div(N, 2)+1:N] .= 1
Iwest = reshape(Iwest, NN, 1)
Ieast = reshape(Ieast, NN, 1)


# Border friction between grid points;
bord = ones(NN, NN)
bord .= 2
for i in 1:length(diag(bord))
    bord[i,i] = 1
end


# Border friction between countries
bordcty = ones(NN, NN)

for i in (div(NN, 2)+1):NN
    for j in (1:(div(NN, 2)))
        bordcty[i,j] = 2
        bordcty[j,i] = 2
    end
end

# Counterfactual border friction between grid points;
cbord=ones(NN,NN);

# Counterfactual border friction between countries;
cbordcty=ones(NN,NN);

# **************************;
# **** Parameterization ****;
# **************************;

# Share of goods in consumption expenditure (1-housing share);
alpha=0.75; 
# Elasticity of substitution;
sigma=5;

# ************************************;
# **** Random productivity shocks ****;
# ************************************;
a = randn(NN, 1) .* 1 .+ 0
a = exp.(a)
a[1:div(NN, 2),1] = a[1:div(NN, 2),1] ./ mean(a[1:div(NN, 2),1])
a[(div(NN, 2)+1):NN,1] = a[(div(NN, 2)+1):NN,1] ./ mean(a[(div(NN, 2)+1):NN,1])


display("Summary statistics productivities");
display("mean(a) std(a) max(a) min(a)");
[mean(a) std(a) maximum(a) minimum(a)]

# **************************;
# **** Other Parameters ****;
# **************************:

# Land area;
H=100 .* ones(NN,1);
# Aggregate labor Supply;
LL=153889;          # US civilian labor force 2010 (Statistical Abstract, millions);

LLwest=(sum(Iwest, dims=1)./(sum(Iwest, dims=1)+sum(Ieast, dims=1))).*LL;
LLeast=(sum(Ieast, dims=1)./(sum(Iwest, dims=1)+sum(Ieast, dims=1))).*LL;
# Fixed production cost;
F=1;

# ********************************;
# **** Matrix of Fundamentals ****;
# ********************************;

fund = zeros(NN, 4)
fund[:, 1] = a
fund[:, 2] = H
fund[:, 3] = Iwest
fund[:, 4] = Ieast

# ****************************************************************************;
# **** Open Economy Solve for Endogenous Variables in Initial Equilibrium ****;
# ****************************************************************************;
using SpecialFunctions, Statistics, BenchmarkTools
include("functions/solveHLwCtyOpen_E.jl")


display(">>>> Start Wage and Population Convergence <<<<");
w, L, tradesh, dtradesh, converge, xtic = solveHLwCtyOpen_E(fund,dist,bord,bordcty,NN);
display(">>>> Wage and Population System Converged <<<<");
display(">>>> Check Wage and Population Convergence (Should ==1) <<<<");

[converge]
display(">>>> Elapsed Time in Seconds <<<<");
xtic

# Price index;
include("functions/Hpindex.jl")
P = Hpindex(fund,L,w,dtradesh);

# Land prices;
include("functions/Hlandprice.jl")
r=Hlandprice(fund,L,w);

# Real wage;
include("functions/Hrealw.jl")
realwage=Hrealw(fund,L,w,tradesh);

# ************************************************************************;
# ***** Counterfactual Eliminating Border Frictions Between Countries ****;
# ************************************************************************;

display(">>>> Start Wage and Population Convergence <<<<");
cw,cL,ctradesh,cdtradesh,cconverge,xtic=solveHLwCtyOpen_E(fund,dist,bord,cbordcty,NN);
display(">>>> Wage and Population System Converged <<<<");
display(">>>> Check Wage and Population Convergence (Should ==1) <<<<");
[cconverge]
display(">>>> Elapsed Time in Seconds <<<<");
xtic

# Price index;
cP = Hpindex(fund,cL,cw,cdtradesh);

# Land prices;
cr=Hlandprice(fund,cL,cw);

# Real wage;
crealwage=Hrealw(fund,cL,cw,ctradesh);

# Welfare gains;
include("functions/Hwelfaregains.jl")
welfgain = Hwelfaregains(ctradesh,tradesh,cL,L);
display(">>>> Welfare Gains <<<<");
welfgain=round.(welfgain.*(10 .^4));
welfgain=welfgain./(10 .^4);
W1 = unique(welfgain)
if length(W1)==1
    W1 = repeat(W1, 2)
end

# **************************************************************************;
# ***** Counterfactual Eliminating Border Frictions Between Grid Points ****;
# **************************************************************************;

display(">>>> Start Wage and Population Convergence <<<<");
ccw,ccL,cctradesh,ccdtradesh,ccconverge,xtic=solveHLwCtyOpen_E(fund,dist,cbord,bordcty,NN);
display(">>>> Wage and Population System Converged <<<<");
display(">>>> Check Wage and Population Convergence (Should ==1) <<<<");
[ccconverge]
display(">>>> Elapsed Time in Seconds <<<<");
xtic

# Price index;
ccP = Hpindex(fund,ccL,ccw,ccdtradesh);

# Land prices;
ccr=Hlandprice(fund,ccL,ccw);

# Real wage;
ccrealwage=Hrealw(fund,ccL,ccw,cctradesh);

# Welfare gains;
welfgain=Hwelfaregains(cctradesh,tradesh,ccL,L);
display(">>>> Welfare Gains <<<<");
welfgain=round.(welfgain.*(10 .^4));
welfgain=welfgain./(10 .^4);
W2 = unique(welfgain)
if length(W2)==1
    W2 = repeat(W2, 2)
end

# ***********************************************;
# **** Three-Dimensional Initial Equilibrium ****;
# ***********************************************;

# LOG PRODUCTIVITY;
amat=reshape(log.(a),N,N);

# LOG POPULATION;
Lmat=reshape(log.(L),N,N);

# LOG WAGE;
wmat=reshape(log.(w),N,N);

# LOG RELATIVE LAND PRICE;
rmat=reshape(log.(r),N,N);

# PRICE INDEX;
Pmat=reshape(log.(P),N,N);

# MULTI-PANEL FIGURE 1;
using Plots
color = cgrad(:haline, rev = false, scale = :level);

# Productivity;
H_cnty_initial_prod = heatmap(amat , title = "Log Productivity", colorbar=true, xlabel="Longitude", ylabel="Latitude" ,  c = color)
savefig(H_cnty_initial_prod,"graph/H_cnty_initial_prod.pdf")


# MULTI-PANEL FIGURE 2;
# Population;
p1 = heatmap(Lmat , title = "Panel A : Log Population", colorbar=true, xlabel="Longitude", ylabel="Latitude" ,  c = color, xguidefontsize=8,yguidefontsize=8,titlefontsize=8)
# Wages;
p2 = heatmap(wmat , title = "Panel B : Log Wages", colorbar=true, xlabel="Longitude", ylabel="Latitude" ,  c = color, xguidefontsize=8,yguidefontsize=8,titlefontsize=8)
# Land prices;
p3 = heatmap(rmat , title = "Panel C : Log Land Prices", colorbar=true, xlabel="Longitude", ylabel="Latitude" ,  c = color, xguidefontsize=8,yguidefontsize=8,titlefontsize=8)
# Price Index;
p4 = heatmap(Pmat , title = "Panel D : Log Price Index", colorbar=true, xlabel="Longitude", ylabel="Latitude" ,  c = color, xguidefontsize=8,yguidefontsize=8,titlefontsize=8)

H_cnty_initial = [p1 p2 p3 p4 ];
H_cnty_initial = plot(H_cnty_initial..., layout = (2,2))
savefig(H_cnty_initial,"graph/H_cnty_initial.pdf")

# **************************************************************************;
# **** Three-Dimensional Eliminating Border Frictions Between Countries ****;
# **************************************************************************;

# POPULATION;
dL=cL./L; ldL=log.(dL);
dLmat=reshape(ldL,N,N);

# PRICE INDEX;
dP=cP./P; ldP=log.(dP);
dPmat=reshape(ldP,N,N);

# WAGE;
dw=cw./w; ldw=log.(dw);
dwmat=reshape(ldw,N,N);

# RELATIVE LAND PRICE;
dr=cr./r; ldr=log.(dr);
drmat=reshape(ldr,N,N);

# MULTI-PANEL FIGURE 3: DIFFERENCES;
# Population;
p5 = heatmap(dLmat , title = "Panel A : Log Relative Population", colorbar=true, xlabel="Longitude", ylabel="Latitude" ,  c = color, xguidefontsize=8,yguidefontsize=8,titlefontsize=8)
# Wage;
p6 = heatmap(dwmat , title = "Panel B : Log Relative Wages", colorbar=true, xlabel="Longitude", ylabel="Latitude" ,  c = color, xguidefontsize=8,yguidefontsize=8,titlefontsize=8)
# Land price;
p7 = heatmap(drmat , title = "Panel C : Log Relative Land Rents", colorbar=true, xlabel="Longitude", ylabel="Latitude" ,  c = color, xguidefontsize=8,yguidefontsize=8,titlefontsize=8)
# Price index;
p8 = heatmap(dPmat , title = "Panel D : Log Relative Price Index", colorbar=true, xlabel="Longitude", ylabel="Latitude" ,  c = color, xguidefontsize=8,yguidefontsize=8,titlefontsize=8)

H_cnty_c = [p5 p6 p7 p8 ];
H_cnty_c = plot(H_cnty_c..., layout = (2,2))
savefig(H_cnty_c,"graph/H_cnty_c.pdf")



# MULTI-PANEL FIGURE 4 : LEVELS;
# Population;
p9 = heatmap(Lmat + dLmat , title = "Panel A : Log Population", colorbar=true, xlabel="Longitude", ylabel="Latitude" ,  c = color, xguidefontsize=8,yguidefontsize=8,titlefontsize=8)
# Wage;
p10 = heatmap(wmat + dwmat , title = "Panel B : Log Wages", colorbar=true, xlabel="Longitude", ylabel="Latitude" ,  c = color, xguidefontsize=8,yguidefontsize=8,titlefontsize=8)
# Land price;
p11 = heatmap(rmat + drmat , title = "Panel C : Log Land Rents", colorbar=true, xlabel="Longitude", ylabel="Latitude" ,  c = color, xguidefontsize=8,yguidefontsize=8,titlefontsize=8)
# Price index;
p12 = heatmap(Pmat + dPmat , title = "Panel D : Log Price Index", colorbar=true, xlabel="Longitude", ylabel="Latitude" ,  c = color, xguidefontsize=8,yguidefontsize=8,titlefontsize=8)

H_cnty_c_lev = [p9 p10 p11 p12 ];
H_cnty_c_lev = plot(H_cnty_c_lev..., layout = (2,2))
savefig(H_cnty_c_lev,"graph/H_cnty_c_lev.pdf")


# ****************************************************************************;
# **** Three-Dimensional Eliminating Border Frictions Between Grid Points ****;
# ****************************************************************************;

# POPULATION;
ddL=ccL./L; lddL=log.(ddL);
ddLmat=reshape(lddL,N,N);

# PRICE INDEX;
ddP=ccP./P; lddP=log.(ddP);
ddPmat=reshape(lddP,N,N);

# WAGE;
ddw=ccw./w; lddw=log.(ddw);
ddwmat=reshape(lddw,N,N);

# RELATIVE LAND PRICE;
ddr=ccr./r; lddr=log.(ddr);
ddrmat=reshape(lddr,N,N);

# MULTI-PANEL FIGURE 5 : DIFFERENCES;
# Population;
p13 = heatmap(ddLmat, title = "Panel A : Log Relative Population", colorbar=true, xlabel="Longitude", ylabel="Latitude" ,  c = color, xguidefontsize=8,yguidefontsize=8,titlefontsize=8)
# Wage;
p14 = heatmap(ddwmat, title = "Panel B : Log Relative Wages (Truncated)", colorbar=true, xlabel="Longitude", ylabel="Latitude" ,  c = color, xguidefontsize=8,yguidefontsize=8,titlefontsize=8)
# Land price;
p15 = heatmap(ddrmat, title = "Panel C : Log Relative Land Rents", colorbar=true, xlabel="Longitude", ylabel="Latitude" ,  c = color, xguidefontsize=8,yguidefontsize=8,titlefontsize=8)
# Price index;
p16 = heatmap(ddPmat, title = "Panel D : Log Relative Price Index (Truncated)", colorbar=true, xlabel="Longitude", ylabel="Latitude" ,  c = color, xguidefontsize=8,yguidefontsize=8,titlefontsize=8)

H_cnty_cc = [p13 p14 p15 p16 ];
H_cnty_cc = plot(H_cnty_cc..., layout = (2,2))
savefig(H_cnty_cc,"graph/H_cnty_cc.pdf")



# MULTI-PANEL FIGURE 6 : LEVELS;
# Population;
p17 = heatmap(Lmat + ddLmat , title = "Panel A : Log Population", colorbar=true, xlabel="Longitude", ylabel="Latitude" ,  c = color, xguidefontsize=8,yguidefontsize=8,titlefontsize=8)
# Wage;
p18 = heatmap(wmat + ddwmat , title = "Panel B : Log Wages", colorbar=true, xlabel="Longitude", ylabel="Latitude" ,  c = color, xguidefontsize=8,yguidefontsize=8,titlefontsize=8)
# Land price;
p19 = heatmap(rmat + ddrmat , title = "Panel C : Log Land Rents", colorbar=true, xlabel="Longitude", ylabel="Latitude" ,  c = color, xguidefontsize=8,yguidefontsize=8,titlefontsize=8)
# Price index;
p20 = heatmap(Pmat + ddPmat , title = "Panel D : Log Price Index", colorbar=true, xlabel="Longitude", ylabel="Latitude" ,  c = color, xguidefontsize=8,yguidefontsize=8,titlefontsize=8)

H_cnty_cc_lev = [p17 p18 p19 p20 ];
H_cnty_cc_lev = plot(H_cnty_cc_lev..., layout = (2,2))
savefig(H_cnty_cc_lev,"graph/H_cnty_cc_lev.pdf")

#WELFARE GAINS OF TRADE LIBERALIZATIONS (INTERNAL vs EXTERNAL)
W1 = (W1 .- 1)*100
W1 = round.(W1, digits=3)
Ext = Vector{Any}(W1)
push!(Ext,"External")
reverse!(Ext)
Ext= reshape(Ext, 1, 3)

W2 = (W2 .- 1)*100
W2 = round.(W2, digits=3)
Int = Vector{Any}(W2)
push!(Int,"Internal")
reverse!(Int)
Int= reshape(Int, 1, 3)

using LaTeXTabulars
using LaTeXStrings              
latex_tabular("graph/table.tex",
               Tabular("lcl"),
               [Rule(:top),
               [MultiColumn(3, :c, "Welfare Gains in pct")], # ragged!
                ["Trade Liberalization" ,"West", "East"],
                Rule(:mid),
                Ext,
                Rule(:mid),
                Int,
                Rule(:bottom)])      
                

            
