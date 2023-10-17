# MCMC のデータを読み込んでPlot する

module My_plot_funcs
using JLD2
using FileIO
using HDF5
using StatsPlots

export Load_MCMC, Plot_MCMC, Check_rhat, Plot_6posteriors

#####
function Load_MCMC(mouseid, date, n_rec, ROI)
    data_path = "$mouseid/$date/Rec$n_rec/$ROI";
    jldopen("./result/M$mouseid"*"_MCMC.jld2", "r") do file
        model = file["$data_path/MCMC/model"];
        mean_MCMC = file["$data_path/Plot/MCMC_mean"];
        lower = file["$data_path/Plot/lower"];
        upper = file["$data_path/Plot/upper"];
        return (model, mean_MCMC, lower, upper)
    end
end

function Plot_MCMC(x, y, x_new, y_MCMC, lower, upper, ROI, model)
    # Raw data
    p = scatter(x, y[:,ROI], ms=2);#, label="Data"); 
    plot!(x_new, y_MCMC, ribbon=(lower, upper));#, label="Estimation");
    #plot!(legend=:bottomleft)
    title!("ROI = $ROI, Model: $model");
    xlabel!("Direction of moving bar");
    ylabel!("z-scored dF/F");

    return p
end

function Plot_6posteriors(s_data, ang, peak, test_ang, startROI)
    plot_list = [];

    for ROI in startROI:startROI+5
        model, mean_MCMC, lower, upper = Load_MCMC(s_data.mouseid, s_data.date, s_data.n_rec, ROI);
        p = Plot_MCMC(ang, peak, test_ang, 
                mean_MCMC, lower, upper, ROI, model)
        push!(plot_list, p)
    end
    #plot_list... の意味は
    plot(plot_list..., layout=(3,2), legend=false)
end
###################
function Get_chain(s, ROI)
    #
    # s_data = Setting(mouseid, date, n_rec);
    #
    mouseid = s.mouseid;
    date = s.date;
    n_rec = s.n_rec;

    f_save_name = "./result/M$mouseid"*"_MCMC.jld2";
    data_path = "$mouseid/$date/Rec$n_rec";

    println("Under developping.")
    data = load(f_save_name);
    ch = data["$data_path/$ROI/MCMC/chain"]
end


#End of My_plot_funcs
end