nsims=20

source("R/fonction_simu_attrition.R")

mat.att <- matrix(nrow = 0, ncol = 0)
mat.se <- matrix(nrow = 0, ncol = 0)

for (simu_i in 1:nsims){
    data=fonction_simu_attrition(theta2_alpha_Gg=0, lambda1_alpha_St=0, sigma_alpha=2, sigma_epsilon=0.5, alpha_percentile=0.9)
    
    chained.results=chained(
                yname="Y1_chaine",
                tname="annee",
                idname="id",
                gname="annee_G",
                xformla=~X, 
                propensityformla=c("X"), 
                data=data,      
                anticipation=0,      
                weightsname=c("P_Y1_chaine"), #St   
                weight_assumption=NULL,
                bstrap=FALSE, 
                biters=1000,
                debT=3,
                finT=8,
                deb=1,
                fin=8,
                select='select',
                treated='traite_G',
                pl=FALSE,
                cores=1,
                cband=TRUE,
                clustervars=NULL)


    agg <- aggte(MP = chained.results, type = 'dynamic')
    if (nrow(mat.att) == 0) {
    mat.att <- matrix(nrow = length(agg$att.egt), ncol = 0)} #we adjust the size of the matrix for cbind.
    if (nrow(mat.se) == 0) {
    mat.se <- matrix(nrow = length(agg$se.egt), ncol = 0)} #we adjust the size of the matrix for cbind.
    
    
    mat.att <- cbind(mat.att, agg$att.egt)
    mat.se <- cbind(mat.se, unlist(agg$se.egt))
}

row_averages = rowMeans(mat.att, na.rm = TRUE)
results.att <- matrix(row_averages, ncol = 1)
row_averages = rowMeans(mat.se, na.rm = TRUE)
results.se <- matrix(row_averages, ncol = 1)
results.final <- cbind(results.att, results.se)

num_rows <- nrow(results.final)
tt <- c(-(num_rows %/% 2):((num_rows - 1) %/% 2))
results.final <- cbind(tt, results.final)
results.final









