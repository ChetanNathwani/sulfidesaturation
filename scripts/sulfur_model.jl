# Generate isobar LLDs for randomly generated fO2 and pressure for the deep crust

using MAGEMin_C, DataFramesMeta, CSV, CairoMakie, JLD2, DelimitedFiles, StatsBase

function compile_arcMagma_db(n_X, n_T)

    # n_X  Number of starting compositions you want in the subsampled dataframe
    # n_T  Number of temperature points to run for each composition

    # Load compositions to data frame
    target_dir       = "./data"
    data_path        = "./data/Schmidt_Jagoutz_cleaned.csv"
    data, header     = readdlm(data_path, ',', header=true)
    df_orig          = DataFrame(data, vec(header))
    df_orig          = df_orig[df_orig[!, "SIO2(WT%)"] .< 55, :]
    df_orig[!, "CR2O3(WT%)"] .= 0.01
    df_orig[!, "O(WT%)"] .= 1.49 # Generate oxygen column to saturate O and activate the buffer

    subsampled_df = df_orig[sample(1:nrow(df_orig), n_X, replace=true), :] # Subsample with replacement

    # Generate column for random fO2 values (delta FMQ)
    min_deltaFMQ = 0.0
    max_delta_FMQ = 3.0
    subsampled_df[!, "del_FMQ"] .= min_deltaFMQ .+ (max_delta_FMQ - min_deltaFMQ) .* rand(n_X)

    # Generate column with either 0, 1 or 2 dQFM 
    # subsampled_df[!, "del_FMQ"] .= Float64.(rand(0:2, nrow(subsampled_df)))
    print(subsampled_df)

    # Generate column for pressure to run isobaric cooling (kbar)
    min_P = 4
    max_P = 9
    subsampled_df[!, "P"] .= min_P .+ (max_P - min_P) .* rand(n_X)

    # Generate column for H2O contents
    min_H2O = 0
    max_H2O = 4
    subsampled_df[!, "H2O(WT%)"] .= min_H2O .+ (max_H2O - min_H2O) .* rand(n_X)
    

    # Duplicate dataframe so each composition is duplicated n_T times for the number of temperature points to run for the minimisation
    
    subsampled_df = DataFrame(vcat([[row for _ in 1:n_T] for row in eachrow(subsampled_df)]...))

    # Generate temperatures for minimisation
    min_T = 700
    max_T = 1400
    Ts = collect(range(max_T, stop=min_T, length=n_T))
    Ts = repeat(Ts, n_X) # Repeat this for each of the compositions

    subsampled_df[!, "T"] .= Ts

    return subsampled_df
end

function runMAGEMin(df)
    # Prepare dataframe for MAGEMin 
    columns_to_select = [:"SIO2(WT%)", :"AL2O3(WT%)", :"CAO(WT%)", :"MGO(WT%)", :"FEOT(WT%)", :"K2O(WT%)", :"NA2O(WT%)", :"TIO2(WT%)", "CR2O3(WT%)", :"H2O(WT%)",:"O(WT%)"]
    X_all = df[!, columns_to_select]
    X_all = rename!(X_all, [:"SiO2", :"Al2O3", :"CaO", :"MgO", :"FeO", :"K2O", :"Na2O", :"TiO2", :"Cr2O3", :"H2O", :"O"])
    # X_all = [collect(X_all.SiO2), collect(X_all.Al2O3), collect(X_all.CaO), collect(X_all.MgO), collect(X_all.FeO), collect(X_all.K2O), collect(X_all.Na2O), collect(X_all.TiO2), collect(X_all.Cr2O3), collect(X_all.H2O), collect(X_all.O)] # Convert to nested vector for MAGEMin
    
    # Convert DataFrame rows to a nested vector (each row as a separate vector)
    X_all = [collect(row) for row in eachrow(X_all)]
    
    T_all = df[!, "T"]
    P_all = df[!, "P"]
    B_all = df[!, "del_FMQ"]

    println(length(P_all))
    
    # MAGEMin
    data           = Initialize_MAGEMin("ig", verbose = false, buffer = "qfm", solver=2)
    oxides_MAGEMin = ["SiO2", "Al2O3", "CaO", "MgO", "FeO", "K2O", "Na2O", "TiO2", "Cr2O3", "H2O", "O"]
    sys_in         = "wt"


    
    
    # GC.gc()
    # GC.enable(false)
    Out_all = multi_point_minimization(P_all, T_all, data, X = X_all, B = B_all, Xoxides = oxides_MAGEMin, sys_in = sys_in, progressbar = true)
    return Out_all
end

function postprocess_MAGEMin(Out_PT, n_X, n_T)

    n = n_X*n_T

    T_mage = Vector{Float64}(undef, n)
    P_mage = Vector{Float64}(undef, n)
    dQFM_mage = Vector{Float64}(undef, n)
    melt_fraction_mage = Vector{Float64}(undef, n)

    SiO2_liq_mage  = Vector{Float64}(undef, n)
    Al2O3_liq_mage  = Vector{Float64}(undef, n)
    CaO_liq_mage  = Vector{Float64}(undef, n)
    MgO_liq_mage  = Vector{Float64}(undef, n)
    FeO_liq_mage  = Vector{Float64}(undef, n)
    K2O_liq_mage  = Vector{Float64}(undef, n)
    Na2O_liq_mage  = Vector{Float64}(undef, n)
    TiO2_liq_mage  = Vector{Float64}(undef, n)
    H2O_liq_mage  = Vector{Float64}(undef, n)

    SiO2_sol_mage  = Vector{Float64}(undef, n)
    Al2O3_sol_mage  = Vector{Float64}(undef, n)
    CaO_sol_mage  = Vector{Float64}(undef, n)
    MgO_sol_mage  = Vector{Float64}(undef, n)
    FeO_sol_mage  = Vector{Float64}(undef, n)
    K2O_sol_mage  = Vector{Float64}(undef, n)
    Na2O_sol_mage  = Vector{Float64}(undef, n)
    TiO2_sol_mage  = Vector{Float64}(undef, n)
    H2O_sol_mage  = Vector{Float64}(undef, n)

        
    # end
    for (iOut, out_point) in enumerate(Out_PT)
        T_mage[iOut]  = Out_PT[iOut].T_C
        P_mage[iOut]  = Out_PT[iOut].P_kbar
        dQFM_mage[iOut] = Out_PT[iOut].dQFM

        # Major oxides liquid
        SiO2_liq_mage[iOut]  = Out_PT[iOut].bulk_M_wt[1]
        Al2O3_liq_mage[iOut] = Out_PT[iOut].bulk_M_wt[2]
        CaO_liq_mage[iOut]   = Out_PT[iOut].bulk_M_wt[3]
        MgO_liq_mage[iOut]   = Out_PT[iOut].bulk_M_wt[4]
        FeO_liq_mage[iOut]   = Out_PT[iOut].bulk_M_wt[5]
        K2O_liq_mage[iOut]   = Out_PT[iOut].bulk_M_wt[6]
        Na2O_liq_mage[iOut]  = Out_PT[iOut].bulk_M_wt[7]
        TiO2_liq_mage[iOut]  = Out_PT[iOut].bulk_M_wt[8]
        H2O_liq_mage[iOut]   = Out_PT[iOut].bulk_M_wt[11]
        # Major oxides solid
        SiO2_sol_mage[iOut]  = Out_PT[iOut].bulk_S_wt[1]
        Al2O3_sol_mage[iOut] = Out_PT[iOut].bulk_S_wt[2]
        CaO_sol_mage[iOut]   = Out_PT[iOut].bulk_S_wt[3]
        MgO_sol_mage[iOut]   = Out_PT[iOut].bulk_S_wt[4]
        FeO_sol_mage[iOut]   = Out_PT[iOut].bulk_S_wt[5]
        K2O_sol_mage[iOut]   = Out_PT[iOut].bulk_S_wt[6]
        Na2O_sol_mage[iOut]  = Out_PT[iOut].bulk_S_wt[7]
        TiO2_sol_mage[iOut]  = Out_PT[iOut].bulk_S_wt[8]
        H2O_sol_mage[iOut]   = Out_PT[iOut].bulk_S_wt[11]

        liq_id = findfirst(x-> x .== "liq", out_point.ph) # Find index of liquid
        if isnothing(liq_id)
            melt_fraction_mage[iOut] = 0
            continue
        end
        melt_fraction_mage[iOut] = Out_PT[iOut].ph_frac[liq_id] # Assign melt fraction wt% to df

    end

    # Store MAGEMin data
    df_out = DataFrame()
    df_out[:, "T_C"]              = T_mage
    df_out[:, "P_kbar"]           = P_mage
    df_out[:, "delta_QFM"]        = dQFM_mage
    df_out[:, "melt_frac_wt"]     = melt_fraction_mage
    df_out[:, "SiO2_liq_mage"]    = SiO2_liq_mage
    df_out[:, "Al2O3_liq_mage"]   = Al2O3_liq_mage
    df_out[:, "CaO_liq_mage"]     = CaO_liq_mage
    df_out[:, "MgO_liq_mage"]     = MgO_liq_mage
    df_out[:, "FeO_liq_mage"]     = FeO_liq_mage
    df_out[:, "K2O_liq_mage"]     = K2O_liq_mage
    df_out[:, "Na2O_liq_mage"]    = Na2O_liq_mage
    df_out[:, "TiO2_liq_mage"]    = TiO2_liq_mage
    df_out[:, "H2O_liq_mage"]     = H2O_liq_mage
    df_out[:, "SiO2_sol_mage"]    = SiO2_sol_mage
    df_out[:, "Al2O3_sol_mage"]   = Al2O3_sol_mage
    df_out[:, "CaO_sol_mage"]     = CaO_sol_mage
    df_out[:, "MgO_sol_mage"]     = MgO_sol_mage
    df_out[:, "FeO_sol_mage"]     = FeO_sol_mage
    df_out[:, "K2O_sol_mage"]     = K2O_sol_mage
    df_out[:, "Na2O_sol_mage"]    = Na2O_sol_mage
    df_out[:, "TiO2_sol_mage"]    = TiO2_sol_mage
    df_out[:, "H2O_sol_mage"]     = H2O_sol_mage

    return df_out
end

n_X = 1000 # Number of starting compositions to sample (with replacement)
n_T = 25 # Number of temperature points to run
df = compile_arcMagma_db(n_X, n_T)
Out_all = runMAGEMin(df)
df_out = postprocess_MAGEMin(Out_all, n_X, n_T)
CSV.write("./data/S_model_LLDs_delFMQ0-2discrete.csv", df_out)