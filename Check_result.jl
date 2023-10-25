#Check MCMC result
push!(LOAD_PATH, "./code")

using My_MCMC_funcs
using My_plot_funcs
using StatsPlots


##########
mouseid = 1441;
date = 20230731;
n_rec = 2;

s_data = Setting(mouseid, date, n_rec);

##########
ang, peak = Load_Data(mouseid, date, n_rec);
test_ang = range(0, 2π, 300);

##########

Plot_6posteriors(s_data, ang, peak, test_ang, 1)

s1 = Setting(mouseid, date, 1)
s2 = Setting(mouseid, date, 2)
Compare_plot(s1, s2, 1)