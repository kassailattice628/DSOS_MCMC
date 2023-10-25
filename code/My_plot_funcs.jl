# MCMC のデータを読み込んでPlot する

module My_plot_funcs
using JLD2
using FileIO
using HDF5
using StatsPlots
using My_MCMC_funcs: Load_Data

export Load_MCMC, Plot_MCMC, Check_rhat, Plot_6posteriors, Compare_plot

#####
function Load_MCMC(mouseid, date, n_rec::Int, ROI)    
    
    # Data path in jld2 file.
    data_path = "$mouseid/$date/Rec$n_rec/$ROI";

    # Extract data from jld2 file
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

    # Postetrior
    plot!(p, x_new, y_MCMC, ribbon=(lower, upper));#, label="Estimation");

    # Add info.
    title!("ROI = $ROI, Model: $model");
    xlabel!("Direction of moving bar");
    ylabel!("z-scored dF/F");

    return p
end



function Plot_6posteriors(s_data, ang, peak, test_ang, startROI)
    plot_list = [];

    mouseid = s_data.mouseid;
    date = s_data.date;
    n_rec = s_data.n_rec;

    for ROI in startROI:startROI+5

        # Load MCMC result
        model, m, l, u = Load_MCMC(mouseid, date, n_rec, ROI);

        # Raw data
        p_new = scatter(ang, peak[:,ROI], ms=2);
        plot!(p_new, test_ang, m, ribbon=(l, u));

        # Add into list
        push!(plot_list, p_new)
    end

    #plot_list... の意味は
    plot(plot_list..., layout=(3,2), legend=false)
end


function Compare_plot(s1, s2, startROI)

    # Load data
    ang1, peak1 = Load_Data(s1.mouseid, s1.date, s1.n_rec);
    ang2, peak2 = Load_Data(s2.mouseid, s2.date, s2.n_rec);

    plot_list = [];
    test_ang = range(0, 2π, 300);

    for ROI in startROI:startROI+5

        #p_new = plot();

        p_new = Plot_Raw_MCMC(ang1, peak1, s1, ROI, test_ang, [])
        p_new = Plot_Raw_MCMC(ang2, peak2, s2, ROI, test_ang, p_new)

        push!(plot_list, p_new)
    end
    plot(plot_list..., layout=(3,2), legend=false)
end

function Plot_Raw_MCMC(ang, peak, s, ROI, test_ang, p)
    
    if isempty(p)
        # generate empty plot
        p = plot();
    end
    #Load Model
    model, m, l, u = Load_MCMC(s.mouseid, s.date, s.n_rec, ROI);
    scatter!(p, ang, peak[:,ROI], ms=2);
    plot!(p, test_ang, m, ribbon=(l, u));
    
    return p

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