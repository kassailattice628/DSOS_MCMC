module My_MCMC_funcs
# DS/OS 判定のために
# Free angle stimulus の結果を
# 1) Single Von-Mises (VM)
# 2) Sum of two Von-Mises
# 3) Flat (constant + noise) の 3 モデルを使って，パラメタを MCMC でベイズ推定する．

using DataFrames
using Distributions: Cauchy, Normal, MvNormal, Truncated, Uniform, VonMises
using HDF5
using JLD2
using LogExpFunctions: logsumexp
using NaNMath
using Random: seed!
using SpecialFunctions: besseli
using Statistics
using StatsFuns
using StatsPlots
using Turing
using FileIO

export Setting, Setting_MCMC, Run_main, Load_Data,
    f_1VM, f_2VM


######## Parameters ##########
struct Setting
    mouseid::Int
    date::Int
    n_rec::Int
end

struct Setting_MCMC
    n_warm_up::Int
    n_sampling::Int
    n_threads::Int
    seed::Int
end
##############################

######## File IO ###########################################
#
# Original data is save as HDF5 from Matlab
#
############################################################
function Load_Data(mouseid, date, n_rec)
    # Read HDF5 data which is made by Matlab
    # export
    #   ang: moving bar direction (rad)
    #   peak: peak dF/F (z-scored?)
    # 

    f_h5_name = "./data/M$mouseid.hdf5";
    f_h5 = h5open(f_h5_name, "r");

    data_path_h5 = "/$mouseid/$date/Rec$n_rec/Analysys"
    #HDF5 File
    ang =  read(f_h5, "$data_path_h5/angle");
    peak = read(f_h5, "$data_path_h5/peak");

    if ndims(peak) == 3
        #when stim is "Rand12"
        ang = repeat(ang, size(peak, 1));
        ang = reshape(ang, :, 1);

        p = reshape(peak, :,1);
        p = reshape(p, :, size(peak,3)); #刺激サイズで reshapeしてた．修正
        peak = p;
    end

    close(f_h5)

    if size(ang,1) == 1
        ang = ang';
    end

    return ang, peak
end

function SaveData(f_save_name, data_path, ROI,
            selected_model, chain, waic,
            mean_MCMC, lower, upper)
    # Save posterior distribution of MCMC
    #

    jldopen(f_save_name, "a+") do file
        #If the data is already tested and 
        if haskey(file, "$data_path/$ROI/MCMC/model")
            close(file)
            UpdateData(f_save_name, data_path, ROI,
                selected_model, chain, waic, mean_MCMC, lower, upper)
        else
            if selected_model != "Flat"
                file["$data_path/$ROI/MCMC/chain"] = chain;
                file["$data_path/$ROI/MCMC/waic"] = waic;
            end
            file["$data_path/$ROI/MCMC/model"] = selected_model;
            file["$data_path/$ROI/Plot/MCMC_mean"] = mean_MCMC
            file["$data_path/$ROI/Plot/lower"] = lower;
            file["$data_path/$ROI/Plot/upper"] = upper;
        end
    end
end

function UpdateData(f_save_name, data_path, ROI,
    selected_model, chain, waic, mean_MCMC, lower, upper)
    # To update exsistent data use "load" and "save"
    # Overwirte old data to result from new MCMC sampling, so that jld2 file is not editable
    #

    data = load(f_save_name);
    if selected_model != "Flat"
        data["$data_path/$ROI/MCMC/chain"] = chain;
        data["$data_path/$ROI/MCMC/waic"] = waic;
    end
    data["$data_path/$ROI/MCMC/model"] = selected_model;
    data["$data_path/$ROI/Plot/MCMC_mean"] = mean_MCMC
    data["$data_path/$ROI/Plot/lower"] = lower;
    data["$data_path/$ROI/Plot/upper"] = upper;

    println("Update $f_save_name for ROI:$ROI...")
    save(f_save_name, data);
end

function OpenCheckData(f_save, data_path, ROI)
    # When running MCMC in batch condition.
    # if skip == false,
    # the ROI is going to re-analyze and overwrite exsistent jld2 file.
    #
    println("check: $data_path")
    if isfile(f_save)
        f = jldopen(f_save, "r")
        if haskey(f, "$data_path/$ROI")
            skip = true;
        else
            skip = false;
        end
    else
        skip = false;
    end

    return skip
end


######## Analyze data ###########################################
#
# Models for MCMC, Sampling, Cluculate WAIC
#
#################################################################

function Get_VectorSum(ang, peak, ROI)
    # Remove NaN peak and estimate preferred angle of specific ROI
    # ang, peak is extracted from "LoadData"

    # Get index of the vector after removing NaN
    i = findall(x-> !isnan.(x), peak[:, ROI]);

    # Calculate preferred direction as vector sum.
    z = sum(peak[i, ROI] .* exp.(im .* ang[i]));

    # Extract DSI, Preferred angle
    r = abs(z)/length(ang); #DSI
    θ = rem2pi(angle(z), RoundDown); #Preferred direction (rad)
    return r, θ, i
end

function Get_VectorSum_Orientation(ang, peak, ROI)
    # Remove NaN peak and estimate preferred angle of specific ROI
    # ang, peak is extracted from "LoadData"
    # if the ROI is OS, preferred disction does not work.
    # So, the preferred orientation is used to estimate θ₀.
    
    i = findall(x-> !isnan.(x), peak[:, ROI]);
    z = sum(peak[i, ROI] .* exp.(im .* 2*ang[i]));

    r = abs(z)/length(ang); #OSI
    θ = rem2pi(angle(z), RoundDown); #Preferred orientation (rad)
    θ = θ/2;
    return r, θ, i
end


########## Models ########################################
# Prepare ranged prior Distributions
# "Truncated limit the range of parameters
#
##########################################################

HalfCauchy() = Truncated(Cauchy(), 0.0, Inf) # for σ
truncatedVonMises(θ,κ) = Truncated(VonMises(θ,κ), θ-π/6, θ+π/6); #<-2π/3
truncatedVonMises_2nd(θ,κ) = Truncated(VonMises(θ,κ), θ-π/2, θ+π/2)

f_1VM(ang, R, θ₀, κ₁, base) = 
    R*exp(κ₁ * cos(ang - θ₀))/(2π * besseli(0, κ₁)) + base;

f_2VM(ang, R, p, θ₀, θ₁, κ₁, κ₂, base) = 
    R*(p*exp(κ₁ * cos(ang - θ₀))/(2π * besseli(0, κ₁)) +
        (1-p)*exp(κ₂ * cos(ang - (θ₀+θ₁)))/(2π * besseli(0, κ₂))
    ) + base;

########## Model 1 ##########
@model function Single_Von_Mises(ang, peak, θ)
    # Prior distributions
    θ₀ ~ truncatedVonMises(θ, 3.0)
    κ₁ ~ Uniform(0.1, 30)
    R ~ Uniform(minimum(peak), 1.2*maximum(peak))
    base ~ Uniform(0, mean(peak))
    σ ~ HalfCauchy()

    ######
    peak ~ MvNormal(f_1VM.(ang, R, θ₀, κ₁, base), σ)
end

########## Model 2 ##########
@model function Sum_of_two_Von_Mises_model(ang, peak, θ)
    # Prior distribitions
    #θ₀ ~ truncatedVonMises(θ, 3)
    θ₀ ~ Uniform(θ-π/6, θ+π/6)
    θ₁ ~ truncatedVonMises_2nd(π, 3.0)
    κ₁ ~ Uniform(0.1, 30)
    κ₂ ~ Uniform(0.1, 30)
    p ~ Uniform(0.55, 0.99)

    R ~ Uniform(minimum(peak), 1.2*maximum(peak))
    base ~ Uniform(0, mean(peak))
    σ ~ HalfCauchy()
    
    ######
    peak ~ MvNormal(f_2VM.(ang, R, p, θ₀, θ₁, κ₁, κ₂, base), σ)
end

########## Model 3 ##########
@model function Flat(peak)
    # Prior distributions
    μ ~ Normal(mean(peak), std(peak))
    σ ~ HalfCauchy()

    ######
    peak ~ Normal(μ, σ)
end

######### MCMC and Comapre models using WAIC ##########
function Run_MCMC_3models(ang, peak, ROI,
    n_warm_up=1000, n_sampling=5000, n_threads=4)
    
    println("Calculating ROI:$ROI...")
    
    DSI, θ_dir, ind = Get_VectorSum(ang, peak, ROI);
    OSI, θ_ori, ~ = Get_VectorSum_Orientation(ang, peak, ROI);

    # Extract Nan-removed index:ind of ROI
    ang_select = ang[ind];
    peak_select = peak[ind, ROI];
    
    #model = [];
    chain = Array{Chains}(undef, 3); #chain = [];
    waic = Array{Float64}(undef, 3); #waic = [];

    #for i in 1:3
    Threads.@threads for i in 1:3
        println("Start sampling for ROI:$ROI, Model:$i...")

        if i == 1
            println("Prior θdir = $θ_dir.")
            model_ = Single_Von_Mises(ang_select, peak_select, θ_dir);
        elseif i == 2
            if θ_dir > π
                θ_ori = θ_ori + π;
            end
            println("Prior θori = $θ_ori.")
            model_ = Sum_of_two_Von_Mises_model(ang_select, peak_select, θ_ori);
        elseif i ==  3
            model_ = Flat(peak_select);
        end

        #push!(model, model_);

        chain_ = sample(model_, NUTS(n_warm_up, 0.65), MCMCThreads(), n_sampling, n_threads);
        waic_ = Get_WAIC(model_, chain_);

        #push!(chain, chain_);
        #push!(waic, waic_);
        chain[i] = chain_;
        waic[i] = waic_;
    end
    
    return chain, waic
end

########## WAIC ################################
#
# https://zenn.dev/yng/articles/turing_intro
#
################################################

logmeanexp(x) = logsumexp(x) - log(length(x))

function Get_WAIC(model, chain)
    model_params = chain.name_map[:parameters]
    lppd = pointwise_loglikelihoods(model, chain[model_params])
    lppd = values(lppd)
    pointwise_waic = -2*(logmeanexp.(lppd) - var.(lppd))
    return sum(pointwise_waic)
end

########## Plot data ###########################
#
# Extracte plot data from MCMC chains
#
################################################

function Extract_MCMC_posterior(test_ang, selected_model, chain)

    df = DataFrame(chain)

    CIrange = (0.025, 1-0.025);

    if selected_model == 1
        # Generate array posterior distributions
        arr = [f_1VM.(x, df.R, df.θ₀, df.κ₁, df.base) for x in test_ang];

        #Mean
        mean_MCMC = [mean(v) for v in arr];
        #95% range
        quantiles = [quantile(v, [0.025, 0.975]) for v in arr];
        lower = [m - q[1] for (q, m) in zip(quantiles, mean_MCMC)];
        upper = [q[2] - m for (q, m) in zip(quantiles, mean_MCMC)];

    elseif selected_model == 2
        arr = [f_2VM.(x, df.R, df.p, df.θ₀, df.θ₁, df.κ₁, df.κ₂, df.base) for x in test_ang];

        mean_MCMC = [mean(v) for v in arr];
        quantiles = [quantile(v, [0.025, 0.975]) for v in arr];
        lower = [m - q[1] for (q, m) in zip(quantiles, mean_MCMC)];
        upper = [q[2] - m for (q, m) in zip(quantiles, mean_MCMC)];

    elseif selected_model == 3
        
        mean_MCMC = mean(df.μ) .* ones(length(test_ang),1)
        quantiles = quantile(df.μ, [0.025, 0.975]);
        lower = mean_MCMC[1] - quantiles[1];
        upper = quantiles[2] - mean_MCMC[1];
    end

    return mean_MCMC, lower, upper
end


########## Run batch ###############################
function Run_main(setting_DATA, setting_MCMC, selected_ROIs="All", ForceReDo=false::Bool)
    # setting includs
    #   mouseid, date, n_rec: for
    # setting_MCMC includes sampling parameters
    #   n_warm_up for NUTS sampling
    #   n_sampling, n_threads

    mouseid = setting_DATA.mouseid;
    date = setting_DATA.date;
    n_rec = setting_DATA.n_rec;

    println("Loading data.....")
    ang, peak = Load_Data(mouseid, date, n_rec);
    n_ROI = size(peak, 2);
    
    # NUTS sampling parameters.
    n_warm_up = setting_MCMC.n_warm_up;
    n_sampling = setting_MCMC.n_sampling;
    n_threads = setting_MCMC.n_threads;
    seed!(setting_MCMC.seed)

    # For Output file.
    f_save_name = "./result/M$mouseid"*"_MCMC.jld2";
    data_path = "$mouseid/$date/Rec$n_rec";
    test_ang = range(0, 2π, 300);

    ######### Batch: Runinng MCMC ##########
    #
    #
    if selected_ROIs == "All"
        run_ROIs = 1:n_ROI
    else
        run_ROIs = selected_ROIs
    end

    for ROI in run_ROIs

        if ForceReDo
            skip = false;
        else
            #Check file
            skip = OpenCheckData(f_save_name, data_path, ROI)
        end

        if skip
            println("Data already analyzed.")
            println("Skip ROI#$ROI...")
        else

            ##########
            @time chain, waic = Run_MCMC_3models(
                ang, peak, ROI, 
                n_warm_up, n_sampling, n_threads)

            #Select Best Models (find model with minimum WAIC)
            ~, selected_model = findmin(waic)
            if selected_model == 1
                SelectedModel = "1VM";
            elseif selected_model == 2
                SelectedModel = "2VM"
            elseif selected_model == 3
                SelectedModel = "Flat"
            end

            println("$SelectedModel is selected.")
            # Check rhat for valid posterior distributions
            check_list = Check_rhat(chain[selected_model])
            if sum(check_list) > 0
                println("Selected model is bad convergence of posetiors...")
                SelectedModel = "BAD_$SelectedModel";
            end

            # Extract Plot data
            mean_MCMC, lower, upper = Extract_MCMC_posterior(test_ang, selected_model, chain[selected_model]);

            println("Saving data of ROI:$ROI.....\n")
            SaveData(f_save_name, data_path, ROI,
                SelectedModel, chain[selected_model], waic, mean_MCMC, lower, upper);
        end
    end 
    ####    # 
    println("Finished Run_main_MCMC.\n")
    return nothing
end


function Check_rhat(chain)
    #
    # check whether the result of MCMC sampling is valid or not
    # output: if check_list contains 1, those posterior distributions were not converged..
    #
    a = summarize(chain);
    params = a[:, :parameters];
    rhat = a[:, :rhat];
    check_list = [];
    for i in eachindex(params)
        if rhat[i] < 1.1
            push!(check_list, 0)
        else
            push!(check_list, 1)
        end
    end
    #
    return check_list
end



# END of Module: My_MCMC_funcs
end