#Check MCMC result
push!(LOAD_PATH, "./code")

using My_MCMC_funcs
using My_plot_funcs
using StatsPlots


##########
mouseid = 1441;
date = 20230731;
n_rec = 1;

s_data = Setting(mouseid, date, n_rec);

##########
ang, peak = Load_Data(mouseid, date, n_rec);
test_ang = range(0, 2π, 300);

##########

Plot_6posteriors(s_data, ang, peak, test_ang, 100)