clear all
cap log close 
set more 1 

log using LFP_correlations_from_model_output.log, replace 

cd "C:\Users\serge\Dropbox\MyResearchProjects\Demography Based Taxation\Code\UPDATED\Joint_taxation_Laffer_curve_v6_G_IRIDIS"
infile age gender ID weight maritalstatus assetholdings householdlaborincome Household_Labor_Income_Tax_Paid  Household_consumption_Tax_Paid ability hours earnings using "Results\Simulation_output_iridis.txt"

gen mod_id = ID*10 + gender

sort mod_id age
tsset mod_id age

gen particip=0
replace particip=1 if hours > 0.001
gen l_particip=l.particip
gen l2_particip=l.l_particip
gen l3_particip=l.l2_particip
gen l4_particip=l.l3_particip
gen l5_particip=l.l4_particip

sort gender maritalstatus
by gender maritalstatus: correlate particip l_particip l2_particip l3_particip l4_particip l5_particip

log close