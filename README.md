# GoingNUTS

This repository includes replication codes for the project "Going NUTS: the regional impact of extreme climate events over the medium term"

This replication package provides:

Input files (Data preperation): 
             "final_data_population_nuts3.csv" for weather 
             "Macrodata_NUTS3"  for macro data
             "GVA_NACE_sector"  Gross value added for sectoral analysis


Output files (for data analysis)
              "Final_data.dta" 



Step 1 (Optional): Run the do-file "prepare_data.do". This file cleans the raw data and creates the event variables and subgroups for further analysis.

Step 2 (Analysis): run "main_analysis.do" runs the main analysis (Figures 1, Figure 4 and Figure 5)
       and also conducts sectoral analysis for all three types of events. 

# Set the path
cap ado uninstall setroot
net install setroot, from("https://github.com/MileslParker/GoingNUTS")

