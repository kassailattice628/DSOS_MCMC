#include("./code/My_MCMC_funcs.jl")
#using .My_MCMC_funcs
using CSV
using DataFrames
using Query: @filter

push!(LOAD_PATH, "./code")
using My_MCMC_funcs

f_csv = "./data/list_DSOS_20231013.csv";
@time df = DataFrame(CSV.File(f_csv)); #これが一番早かったので
#@time df = CSV.read(f_csv, DataFrame); #読み込んだテーブルの結果は同じ
#@time df3 = CSV.File(f_csv) |> DataFrame;

####### Example ##########
# mouseid = 1441;
# date = 20230731;
# n_rec = 1;
##########################
n_warm_up = 1000;
n_sampling = 5000;
n_threads = 4;
seed = 628;
s_MCMC = Setting_MCMC(n_warm_up, n_sampling, n_threads, seed);

mouse_list = unique(df.mouse);
println(mouse_list)

function ReDo_DF(df=df)
    mouse_= [1399, 1402, 1433, 1452];
    date_ = [[20230529], [20230512], [20230606], [20230622]];
    n_ = [[2,4,5], [2,3], [2], [1]];

    #df_ = filter(:mouse => n -> n == 0, df);
    df_new = empty(df); #df と同じ列名を持つからの DF を作成

    for i in eachindex(mouse_)
        for ii in eachindex(date_[i])
            for iii in eachindex(n_[i])
                x = df |> 
                    @filter(
                        _.mouse == mouse_[i] &&
                        _.date == date_[i][ii] && 
                        _.n_rec == n_[i][iii]) |>
                    DataFrame
                append!(df_new, x);
            end
        end
    end
    return df_new
end


function Run_batch(mouse_list, s_MCMC, df=df; ForceReDo=false::Bool)
########## Run ##########
for mouse_id in mouse_list
    println("Mouse: $mouse_id")

    df_mouse = filter(:mouse => n -> n == mouse_id, df);
    date_list = unique(df_mouse.date);

    for date in date_list
        df_date = filter(:date => n -> n == date, df_mouse);

        for n_rec in df_date.n_rec
        #Threads.@threads for n_rec in df_date.n_rec
            s_DATA = Setting(mouse_id, date, n_rec)
            println("Mouse:$mouse_id, Date:$date, Rec:$n_rec")
            Run_main(s_DATA, s_MCMC, "All", ForceReDo)
        end
    end
end
end


#df = ReDo_DF(df);
#mouse_list = [1399, 1402, 1433, 1452];
pop!(mouse_list)
Run_batch(mouse_list, s_MCMC, df, ForceReDo=false);