# GoingNUTS

Codes for project "Going NUTS: the regional impact of extreme climate events over the medium term"

Input files (Data preperation): 
             "final_data_population_nuts3.csv" for weather 
             "Macrodata_NUTS3"  for macro data
             "GVA_NACE_sector"  Gross value added for sectoral analysis


Output files (for data analysis)
              "Final_data.dta" 



Step 1 (Optional): run do-file "prepare_data.do" 
        prepare_data.do cleans the raw data, creates the event variables and subgroups.

Step 2 (Analysis) : run "main_analysis.do" for main results and sectoral analysis
       main_analysis.do runs the main analysis (Figures 1, Figure 4 and Figure 5)
       and also conducts sectoral analysis for all three types of events. 