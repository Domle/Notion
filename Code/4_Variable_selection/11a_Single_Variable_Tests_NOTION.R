### SPECIES DISTRIBUTION MODEL (GAM) - UNIVARIATE MODEL TESTS FOR VARIABLE RANKING (SUBSEQUENT SELECTION)

# $Date: 2016-01-07, updated 2020-06-10
#
# Authors: Damiano Righetti, rdamiano@ethz.ch
#
# Description: Single glm-based variable test. 
# Way one: use group-specific target groups for pseudo-absences (overlapping background)
# Way two: use random pseudo-absences (overlapping background)
#
# Variable importance is ranked per taxon and thinning strategy, T-MLD1 stratification is standard
# Questions are e.g.: does Si or Si* have a higher (single) predicive power for pres/abs of plankton?
# Input: lists with p/a data frames (taxa), output: variable ranking
# Note, taxa <15 obs have been excluded perviously in P/A generation step

### =========================================================================
### Initialize system
### =========================================================================

# Set locale
setwd("/Users/rdamiano/Desktop/Data/8_Variable_Selection")
# kryo or clusters
# if (getwd() == "/home/rdamiano/") {setwd("/home/rdamiano/Data/8_Variable_Selection")}

# Input and output
input.dir <- "~/Desktop/Data_notion/4_Generate_absences/Output_1"
output.dir <- "./Output_11"

# Libraries
require(doParallel); require(doBy); require(mgcv); require(randomForest)

# Functions
source("../COMMONS/common_functions.R", chdir = TRUE)

# Function to test goodness of fit
setwd ("..")
# gam uses summary statistics for goodness of fit => summary(gam_full)$r.sq
source(paste0("..", "/PhD/BIO232/4.R_Stats/0.Functions/adj.D2.glm.fun")) # glm, goodness of fit function adj. D-squared, from ecospat package
# rf uses oob mean out of bag error as a "goodness of fit" => mean(rf_full$er[, "OOB"])
setwd("/Users/rdamiano/Desktop/Data/8_Variable_Selection")
 
### =========================================================================
### Load lists with presence-pseudoabsence data, split by species (takes time)
### =========================================================================

# Step one: define index 
list.files(input.dir)


# Concept: 3 background selection strategies * 2 data type strategies
nms.sets <- c(
"gridded_pres_abs_hom,tot_bg", # implementation of a target group (i.e., diazotroph specific) background
"gridded_counts_abs_hom,tot_bg", # implementation of a target group (i.e., diazotroph specific) background
"gridded_pres_abs_hom,gr_bg", # implementation of a target group (i.e., diazotroph specific) background
"gridded_counts_abs_hom,gr_bg", # implementation of a target group (i.e., diazotroph specific) background
"gridded_pres_abs_hom,cr_bg",  # implementation of a cruise-specific (target-cruise) background, this refinement accounts for differences in sampling method
"gridded_counts_abs_hom,cr_bg"  # implementation of a cruise-specific (target-cruise) background, this refinement accounts for differences in sampling method
)

ind <- c(1,3) # ** CRITICAL DECISION: I use the first and third background * data scenario as the basis for the variable test

# For other papers (LETTER)
# ind <- grep("(T_MLD1)", list.files(input.dir)) # define files: Prepared species with backgrounds.
# ind <- c(4,9,11,15,16,20,21,6,7,8,10,25,42) # 7 N-Sal; 8 P-Wind; 10 sit sep T-MLD
# ind <- c(25, 30, 42, 45) # 25 - mt 300 pts overl; 27 - mt 300 sit overl; 42 - mt 600 pts overl; 45 - mt 600 sit overl
# ind <- c(1,2,3)
# ind <- c(9)

# Load data
taxa.lists <- list()
for (l in ind){taxa.lists[[length(taxa.lists)+1]]<-get(load(paste0(input.dir,"/",  list.files(input.dir)[l]  )))} 

# Define vector of variables to test (order following logic in of Sci. Adv. Letter, Righetti et al., 2019)
vc.vars<-c("x","y","T","Sal","N","P","Si",
			"MLD1", "MLD2","PAR","Chl","Wind.CCMP","pCO2",
			"MLPAR1", "MLPAR2", "Nstar", "Sistar",
			"dT_dt", "dN_dt", "dP_dt", "dSi_dt","dMLD1_dt",
			"logMLD1", "logMLD2", "logChl", "logN", "logP", "logSi")
			# I remove NPP, as well reflected in Chl

# Select one data version (not containing thinned data), to get the list of all taxa contained!
taxa.ref.list <- get(load(paste0(input.dir,"/",   list.files(input.dir)   [1]))) # gridded, tot-bg, homogenised taxa, overlapping sites, T-MLD1 => Standard approach 
 
# Define the minimum of presences (obs) required per taxon 
min.obs <- 16
# min.obs <- 50

# Get an index for the selected taxa
index.useful <- which(as.numeric(unlist(sapply(taxa.ref.list, function(x){x1<-NA; if(!is.null(x)){x1<-nrow(x[which(x$obs==1), ])};return(x1)})))>=min.obs)
length(index.useful) # 30 

### =========================================================================
###  Single variable test
### =========================================================================

# Loop across data lists (10 strategies)
# for(d in 1:length(taxa.lists)){
for(d in 1:length(ind) ){

# Get list with taxa
taxa.ls <- taxa.lists[[d]]

# Exclude elements below min.obs based on observation points in the standard taxa list
taxa.list <- taxa.ls[ index.useful ]

		# Parallel computing across taxa
		n.cores <- detectCores()-2; if (getwd()=="/UP_home/rdamiano/Data/8_Variable_Selection"){n.cores<-detectCores()-10}
		cl = makeCluster(n.cores, outfile="")	
		registerDoParallel(cl)
		list.cl<-foreach(i=c(1:length(taxa.list)), .packages = c("doBy", "mgcv", "randomForest"))%dopar%{

				# Single variable glm test: case taxon has zero obs, create NA's
				if( length(taxa.list[[i]])==0 ){
				vec.glm <- rep(NA,length(vc.vars))
				vec.gam <- rep(NA,length(vc.vars))
				vec.rf <- rep(NA,length(vc.vars))
				ranks.glm <- rep(NA,length(vc.vars))
				ranks.gam <- rep(NA,length(vc.vars))
				ranks.rf <- rep(NA,length(vc.vars))
				obs <- c(NA, NA)
				}else{
				
				# Filter taxa with obs < minimum obs
				if(  length(taxa.list[[i]])>0 & nrow(taxa.list[[i]][which(taxa.list[[i]]$obs==1),]) < min.obs ) {
				vec.glm <- rep(NA,length(vc.vars))
				vec.gam <- rep(NA,length(vc.vars))
				vec.rf <- rep(NA,length(vc.vars))
				ranks.glm <- rep(NA,length(vc.vars))
				ranks.gam <- rep(NA,length(vc.vars))
				ranks.rf <- rep(NA,length(vc.vars))
				obs <- c(NA, NA)
				} else {
				
				# Useful taxa: vector to insert value for each variable == k
				
				#
				#
				#
		
				# Loop across variables
				vec.glm <- rep(NA,length(vc.vars))
				vec.gam <- rep(NA,length(vc.vars))
				vec.rf <- rep(NA,length(vc.vars))
				#
				for(k in 1:length(vc.vars)){
				# get data: removing NA's with regard to the variable tested
				dat <- taxa.list[[i]][ which(is.na(taxa.list[[i]][, vc.vars[k]])==F), ]
				if(    nrow(dat[which(dat$obs==1),]) < 14 ){ # we want at least 14 points to fit a model
				vec.gam[k] <- NA
				vec.glm[k] <- NA
				vec.rf[k] <- NA
				next
				}
				# Re-adjust the weighting of the data: not entirely clear to me if weighting is really needed. As adjusted D-squared and adjusted R-squared take into account the sample size?
				dat$weights_dec <- 1
				weight_abs <-length(which(dat$obs == 1))/length(which(dat$obs == 0)) # New ratio of presence count divided by absences counts, in case some presences or absences were lost due to variable constraints
				dat$weights_dec[dat$obs == 0] <- weight_abs 
				# fit models: weighted (note: the weighted models yield more comparable and similar results between gam and glm, tested at the example of species i = 355)
				if(i == 43 & k == 17){vec.gam[k] <- NA; vec.rf[k] <- NA; next}
				mod.gam <- gam(formula(paste0("obs~s(",vc.vars[k] ,",k=5)")), data = dat, na.action=na.omit, family=binomial, weights=weights_dec)
				mod.glm <- glm( dat$obs ~ dat[  ,  vc.vars[k] ] + I(  (dat[,  vc.vars[k] ])^2), data=dat, na.action=na.omit, family=binomial, weights=weights_dec)
				# fit models: unweighted (previously)
				# mod.gam <- gam(formula(paste0("obs~s(",vc.vars[k] ,",k=5)")), data = dat, na.action=na.omit, family=binomial)
				# mod.glm <- glm( dat$obs ~ dat[  ,  vc.vars[k] ] + I(  (dat[,  vc.vars[k] ])^2), data=dat, na.action=na.omit, family=binomial, weights=weights_dec)
				# fit random forest
				if(weight_abs > 1){mod.rf <- randomForest( eval(parse(text=paste0("as.factor(obs)~", vc.vars[k]))), data= dat, importance=T, nodesize=1, ntree=300, sampsize=c( nrow(dat[which(dat $obs!=1),]) , nrow(dat[which(dat $obs==1),])))} else {
				mod.rf <- randomForest( eval(parse(text=paste0("as.factor(obs)~", vc.vars[k]))), data= dat, importance=T, nodesize=1, ntree=300, sampsize=c( nrow(dat[which(dat $obs==1),]) , nrow(dat[which(dat $obs==1),])))}
				# evaluate models (fit)
				vec.gam[k] <- round(summary(mod.gam)$r.sq , 4)
				vec.glm[k] <- round(adj.D2.glm(mod.glm), 4)
				vec.rf[k] <- round(1-mean(mod.rf$er[, "OOB"]), 4)
				}

				#
				#
				#
			
				#
				ranks.glm <- c(NA,NA, rank(-vec.gam[ 3:length(vec.gam)], ties.method=c("max")))
				ranks.gam <- c(NA,NA, rank(-vec.glm[ 3:length(vec.glm)], ties.method=c("max")))
				ranks.rf <- c(NA,NA, rank(-vec.rf[ 3:length(vec.rf)], ties.method=c("max")))
				obs <- rep(      nrow(          dat[which(dat$obs>0),]               ),2)
				}
				}
						
				# Merge test metrics (adjusted R-squared of model), ranking of test metric, obs, stratification-way, taxon name, group name and NA's row for better structure
				df.gam<-data.frame(
				rbind(vec.gam, ranks.gam, c(rep(NA, length(vec.gam)))),
				"obs"= c(obs, NA), 
				"pa_strat"= c(rep(list.files(input.dir)[ind[d]],2),"NA"),
				"value" = c("adj.Rsq", "rank", "NA"),
				"taxon"= c(rep(names(taxa.list)[i],2),"NA"),
				"group"= c(rep(as.character(taxa.list[[i]]$group[1]),2),"NA")
				)
				colnames(df.gam)[1:length(vc.vars)] <- vc.vars; print(paste(i, "gam"))
			
				df.glm<-data.frame(
				rbind(vec.glm, ranks.glm, c(rep(NA, length(vec.glm)))),
				"obs"= c(obs, NA), 
				"pa_strat"= c(rep(list.files(input.dir)[ind[d]],2),"NA"),
				"value" = c("adj.Dsq", "rank", "NA"),
				"taxon"= c(rep(names(taxa.list)[i],2),"NA"),
				"group"= c(rep(as.character(taxa.list[[i]]$group[1]),2),"NA")
				)
				colnames(df.glm)[1:length(vc.vars)] <- vc.vars; print(paste(i, "glm"))
			
				df.rf<-data.frame(
				rbind(vec.rf, ranks.rf, c(rep(NA, length(vec.rf)))),
				"obs"= c(obs, NA), 
				"pa_strat"= c(rep(list.files(input.dir)[ind[d]],2),"NA"),
				"value" = c("1-OOB.error", "rank", "NA"),
				"taxon"= c(rep(names(taxa.list)[i],2),"NA"),
				"group"= c(rep(as.character(taxa.list[[i]]$group[1]),2),"NA")
				)
				colnames(df.rf)[1:length(vc.vars)] <- vc.vars; print(paste(i, "rf"))

				# Define output
				list.mod <- list()
				list.mod[[1]] <- df.gam
				list.mod[[2]] <- df.glm
				list.mod[[3]] <- df.rf
				list.mod
		} #Close parallel computing across taxa
	
		# Stop cluster
		stopCluster(cl)
	
	#-----------------------------------------------------------------------------------------------------------------------------
	# Save species results: gam
	df.new <- as.data.frame(do.call("rbind",   lapply(list.cl, "[[", 1)   ))
	df.new<- cbind(df.new[, 1:2], round(df.new[ ,3:(ncol(df.new)-4)], 3), df.new[,(ncol(df.new)-3):ncol(df.new)]) # round values (optional)
	df.new[is.na(df.new)] <- "________________"
	df.new <- df.new[  , c("x", "y", "T", "pCO2", "P", "logP", "Nstar", "Wind.CCMP", "logChl", "Si", "N", "Sal", "Sistar", "dMLD1_dt", "Chl", "logN", "logSi", "PAR", "MLPAR1", "MLPAR2", "logMLD1", "MLD1", 
	"dN_dt", "logMLD2", "MLD2", "dT_dt", "dP_dt", "dSi_dt", 
	 "obs",  "taxon", "group", "pa_strat", "value")]
	fln <- paste0(output.dir,"/Taxa_Gam_",gsub(".RData","",list.files(input.dir)[ind[d]]),"_minobs_",min.obs,".csv")
	write.csv(df.new,file=fln,row.names=F)
	#
	# Save species results: glm
	df.new <- as.data.frame(do.call("rbind",   lapply(list.cl, "[[", 2)   ))
	df.new<- cbind(df.new[, 1:2], round(df.new[ ,3:(ncol(df.new)-4)], 3), df.new[,(ncol(df.new)-3):ncol(df.new)]) # round values (optional)
	df.new[is.na(df.new)] <- "________________"
	df.new <- df.new[  , c("x", "y", "T", "pCO2", "P", "logP", "Nstar", "Wind.CCMP", "logChl", "Si", "N", "Sal", "Sistar", "dMLD1_dt", "Chl", "logN", "logSi", "PAR", "MLPAR1", "MLPAR2", "logMLD1", "MLD1", 
	"dN_dt", "logMLD2", "MLD2", "dT_dt", "dP_dt", "dSi_dt", 
	 "obs",  "taxon", "group", "pa_strat", "value")]
	fln <- paste0(output.dir,"/Taxa_Glm_",gsub(".RData","",list.files(input.dir)[ind[d]]),"_minobs_",min.obs,".csv")
	write.csv(df.new,file=fln,row.names=F)
	#
	# Save species results: rf
	df.new <- as.data.frame(do.call("rbind",   lapply(list.cl, "[[", 3)   ))
	df.new<- cbind(df.new[, 1:2], round(df.new[ ,3:(ncol(df.new)-4)], 3), df.new[,(ncol(df.new)-3):ncol(df.new)]) # round values (optional)
	df.new[is.na(df.new)] <- "________________"
	df.new <- df.new[  , c("x", "y", "T", "pCO2", "P", "logP", "Nstar", "Wind.CCMP", "logChl", "Si", "N", "Sal", "Sistar", "dMLD1_dt", "Chl", "logN", "logSi", "PAR", "MLPAR1", "MLPAR2", "logMLD1", "MLD1", 
	"dN_dt", "logMLD2", "MLD2", "dT_dt", "dP_dt", "dSi_dt", 
	 "obs",  "taxon", "group", "pa_strat", "value")]
	fln <- paste0(output.dir,"/Taxa_Rf_",gsub(".RData","",list.files(input.dir)[ind[d]]),"_minobs_",min.obs,".csv")
	write.csv(df.new,file=fln,row.names=F)
	#-----------------------------------------------------------------------------------------------------------------------------

	
	#-----------------------------------------------------------------------------------------------------------------------------
	# Save group results (mean and SD of the species per group): gam
	df.all <- as.data.frame(do.call("rbind",   lapply(list.cl, "[[", 1)   ))
	lisi.rows <- list()
	no.taxa<-list()
	tt <-df.all[which(df.all$value == "adj.Rsq"),]
	ss <-df.all[which(df.all$value == "rank"),]
	no.taxa[[length(no.taxa)+1]] <- c(rep(nrow(tt[which(tt$obs>0), ]),4), NA)
	vec <- c(apply( tt[ , 1:(length(tt)-4)], 2,function(x){mean(as.numeric(x[which(x>0)]), na.rm = T)})) # mean across species' adj R2
	vec2 <- c(apply( ss[ , 1:(length(ss)-4)], 2,function(x){mean(as.numeric(x[which(x>0)]))})) # mean across species' rank
	lisi.rows[[length(lisi.rows)+1]] <- vec
	lisi.rows[[length(lisi.rows)+1]] <- c(apply( tt[ , 1:(length(tt)-4)], 2,function(x){sd(as.numeric(x[which(x>0)]))/mean(as.numeric(x[which(x>0)]))})) # relative standard deviation (sd/mean)
	lisi.rows[[length(lisi.rows)+1]] <- if(nrow(tt)>1) {c(rep(NA,2), rank(-vec[3:(length(vec)-1)], ties.method=c("max")),NA)} else {rep(NA, length(vec))} # ranking of mean adj R2 (we skip x, y, obs for ranking, replacing by NA)
	lisi.rows[[length(lisi.rows)+1]] <- if(nrow(tt)>1) {c(rep(NA,2), rank(vec2[3:(length(vec2)-1)], ties.method=c("max")),NA)} else {rep(NA, length(vec))} # ranking of mean rank (we skip x, y, obs for ranking, replacing by NA
	lisi.rows[[length(lisi.rows)+1]] <- rep(NA,length(vec))# NA filler
	df.total <- data.frame(
	round(do.call(rbind, lisi.rows), 3),
	"no_taxa" = do.call(c, no.taxa),
	"group"= c(
	rep("total", 4), "NA"),
	"pa_strat"=  rep("(T_MLD1)",5),
	"value"= (rep(c("mean.adjRsq", "sd/mean","rank_of_mean.adjR2","rank_of_mean.rank","NA"), 1)))
	df.total[is.na(df.total)] <- "________________"
	df.total <- df.total[  , c("x", "y", "T", "pCO2", "P", "logP", "Nstar", "Wind.CCMP", "logChl", "Si", "N", "Sal", "Sistar", "dMLD1_dt", "Chl", "logN", "logSi", "PAR", "MLPAR1", "MLPAR2", "logMLD1", "MLD1", 
	"dN_dt", "logMLD2", "MLD2", "dT_dt", "dP_dt", "dSi_dt", 
	 "obs",  "no_taxa",  "group", "pa_strat", "value")] 
	fln <- paste0(output.dir,"/Total_Gam_",gsub(".RData","",list.files(input.dir)[ind[d]]),"_minobs_",min.obs,".csv")
	fln <- write.csv(df.groups,file=fln,row.names=F)
	#
	# ---------------------------------
	#
	# Save group results (mean and SD of the species per group): glm
	df.all <- as.data.frame(do.call("rbind",   lapply(list.cl, "[[", 2)   ))
	lisi.rows <- list()
	no.taxa<-list()
	tt <-df.all[which(df.all$value == "adj.Dsq"),]
	ss <-df.all[which(df.all$value == "rank"),]	
	no.taxa[[length(no.taxa)+1]] <- c(rep(nrow(tt[which(tt$obs>0), ]),4), NA)
	vec <- c(apply( tt[ , 1:(length(tt)-4)], 2,function(x){mean(as.numeric(x[which(x>0)]), na.rm = T)})) # mean across species' adj R2
	vec2 <- c(apply( ss[ , 1:(length(ss)-4)], 2,function(x){mean(as.numeric(x[which(x>0)]))})) # mean across species' rank
	lisi.rows[[length(lisi.rows)+1]] <- vec
	lisi.rows[[length(lisi.rows)+1]] <- c(apply( tt[ , 1:(length(tt)-4)], 2,function(x){sd(as.numeric(x[which(x>0)]))/mean(as.numeric(x[which(x>0)]))})) # relative standard deviation (sd/mean)
	lisi.rows[[length(lisi.rows)+1]] <- if(nrow(tt)>1) {c(rep(NA,2), rank(-vec[3:(length(vec)-1)], ties.method=c("max")),NA)} else {rep(NA, length(vec))} # ranking of mean adj R2 (we skip x, y, obs for ranking, replacing by NA)
	lisi.rows[[length(lisi.rows)+1]] <- if(nrow(tt)>1) {c(rep(NA,2), rank(vec2[3:(length(vec2)-1)], ties.method=c("max")),NA)} else {rep(NA, length(vec))} # ranking of mean rank (we skip x, y, obs for ranking, replacing by NA
	lisi.rows[[length(lisi.rows)+1]] <- rep(NA,length(vec))# NA filler
	df.total <- data.frame(
	round(do.call(rbind, lisi.rows), 3),
	"no_taxa" = do.call(c, no.taxa),
	"group"= c(
	rep("total", 4), "NA"),
	"pa_strat"=  rep("(T_MLD1)",5),
	"value"= (rep(c("mean.adjDsq", "sd/mean","rank_of_mean.adjD2","rank_of_mean.rank","NA"), 1)))
	df.total[is.na(df.total)] <- "________________"
	df.total <- df.total[  , c("x", "y", "T", "pCO2", "P", "logP", "Nstar", "Wind.CCMP", "logChl", "Si", "N", "Sal", "Sistar", "dMLD1_dt", "Chl", "logN", "logSi", "PAR", "MLPAR1", "MLPAR2", "logMLD1", "MLD1", 
	"dN_dt", "logMLD2", "MLD2", "dT_dt", "dP_dt", "dSi_dt", 
	 "obs",  "no_taxa",  "group", "pa_strat", "value")] 
	fln <- paste0(output.dir,"/Total_Glm_",gsub(".RData","",list.files(input.dir)[ind[d]]),"_minobs_",min.obs,".csv")
	fln <- write.csv(df.groups,file=fln,row.names=F)
	#
	# ---------------------------------
	#
	# Save group results (mean and SD of the species per group): rf
	df.all <- as.data.frame(do.call("rbind",   lapply(list.cl, "[[", 3)   ))
	lisi.rows <- list()
	no.taxa<-list()
	tt <-df.all[which(df.all$value == "1-OOB.error"),]
	ss <-df.all[which(df.all$value == "rank"),]	
	no.taxa[[length(no.taxa)+1]] <- c(rep(nrow(tt[which(tt$obs>0), ]),4), NA)
	vec <- c(apply( tt[ , 1:(length(tt)-4)], 2,function(x){mean(as.numeric(x[which(x>0)]), na.rm = T)})) # mean across species' adj R2
	vec2 <- c(apply( ss[ , 1:(length(ss)-4)], 2,function(x){mean(as.numeric(x[which(x>0)]))})) # mean across species' rank
	lisi.rows[[length(lisi.rows)+1]] <- vec
	lisi.rows[[length(lisi.rows)+1]] <- c(apply( tt[ , 1:(length(tt)-4)], 2,function(x){sd(as.numeric(x[which(x>0)]))/mean(as.numeric(x[which(x>0)]))})) # relative standard deviation (sd/mean)
	lisi.rows[[length(lisi.rows)+1]] <- if(nrow(tt)>1) {c(rep(NA,2), rank(-vec[3:(length(vec)-1)], ties.method=c("max")),NA)} else {rep(NA, length(vec))} # ranking of mean adj R2 (we skip x, y, obs for ranking, replacing by NA)
	lisi.rows[[length(lisi.rows)+1]] <- if(nrow(tt)>1) {c(rep(NA,2), rank(vec2[3:(length(vec2)-1)], ties.method=c("max")),NA)} else {rep(NA, length(vec))} # ranking of mean rank (we skip x, y, obs for ranking, replacing by NA
	lisi.rows[[length(lisi.rows)+1]] <- rep(NA,length(vec))# NA filler
	df.total <- data.frame(
	round(do.call(rbind, lisi.rows), 3),
	"no_taxa" = do.call(c, no.taxa),
	"group"= c(
	rep("total", 4), "NA"),
	"pa_strat"=  rep("(T_MLD1)",5),
	"value"= (rep(c("mean.1-OOB.error", "sd/mean","rank_of_mean.1-OOB.error","rank_of_mean.rank","NA"), 1)))
	df.total[is.na(df.total)] <- "________________"
	df.total <- df.total[  , c("x", "y", "T", "pCO2", "P", "logP", "Nstar", "Wind.CCMP", "logChl", "Si", "N", "Sal", "Sistar", "dMLD1_dt", "Chl", "logN", "logSi", "PAR", "MLPAR1", "MLPAR2", "logMLD1", "MLD1", 
	"dN_dt", "logMLD2", "MLD2", "dT_dt", "dP_dt", "dSi_dt", 
	 "obs",  "no_taxa",  "group", "pa_strat", "value")] 
	fln <- paste0(output.dir,"/Total_Rf_",gsub(".RData","",list.files(input.dir)[ind[d]]),"_minobs_",min.obs,".csv")
	fln <- write.csv(df.groups,file=fln,row.names=F)
	#-----------------------------------------------------------------------------------------------------------------------------

	
# Close loop over data sets
}






##################### PS: *** NOT IMPLEMENTED SO FAR *** Rearrange output for: comparison of variable ranking across thinning/bg strategies : Value == RANK; applied only to absence stratification by T and MLD1 #####################



#

#

## Loop across taxa
# Preparatory
input.dir <- "../7_Generate_Absences/Output"
output.dir <- "./Output_vartest"
df <- read.csv(paste0(output.dir,"/Taxa_(gridded,tot_bg).csv"), header = T, sep = ",")
vec.taxa <- as.character(df[which(!is.na(df$taxon) & df$value == "adj.Rsq" ), ]$taxon)
nms.sets <- c(
"gridded,gr_bg",
"gridded,tot_bg",
"mt_300,gr_bg", 
"mt_300,tot_bg", 
"mt_600,gr_bg", 
"mt_600,gr_bg")
file.list <- list(); for (d in 1:length(nms.sets)){file.list[[length(file.list)+1]]<- read.csv(paste0(output.dir,"/Taxa_(",nms.sets[d],").csv"),header=T, sep=",")}
# 
lisi.rows <- list ()
for (i in vec.taxa){
	# Loop across files: thinning strategies
	for (j in 1:length(file.list)){
	# lisi.rows[[length(lisi.rows)+1]] <- file.list[[j]][which( file.list[[j]]$taxon==i & file.list[[j]]$pa_strat=="(T_MLD1)"),][1,] # adj. Rsq
	lisi.rows[[length(lisi.rows)+1]] <- file.list[[j]][which( file.list[[j]]$taxon==i & file.list[[j]]$pa_strat=="(T_MLD1)"),][2,] # rank adj. Rsq
	# Close loop across files
	} 
	lisi.rows[[length(lisi.rows)+1]] <-rep(NA, 32) # Add filler with NAs
# Close loop across taxa
}
df.tax <- do.call(rbind, lisi.rows)
df.tax$sets <- rep(c(nms.sets, "NA"), length(vec.taxa))
df.tax[is.na(df.tax)] <- "________________"

# Save file
fln <- paste0(output.dir,"/Taxa_(T_MLD1)_vs_thinning.csv")
fln <- write.csv(df.tax,file=fln,row.names=F)

#

#

## Loop across groups
# Preparatory
input.dir <- "../7_Generate_Absences/Output"
output.dir <- "./Output_vartest"
vec.taxa <- c("bacillariophyceae", "chlorophyta", "chrysophyceae", "cryptophyta", "cyanobacteria", "dinoflagellata", "euglenoidea","haptophyta", "raphidophyceae","total")
nms.sets <- c(
"gridded,gr_bg",
"gridded,tot_bg",
"mt_300,gr_bg", 
"mt_300,tot_bg", 
"mt_600,gr_bg", 
"mt_600,gr_bg")
file.list <- list(); for (d in 1:length(nms.sets)){file.list[[length(file.list)+1]]<- read.csv(paste0(output.dir,"/Groups_(",nms.sets[d],").csv"),header=T, sep=",")}
#
lisi.rows <- list ()
for (i in vec.taxa){
	# Loop across files: thinning strategies
	for (j in 1:length(file.list)){
	# lisi.rows[[length(lisi.rows)+1]] <- file.list[[j]][which(file.list[[j]]$group==i & file.list[[j]]$pa_strat=="(T_MLD1)"),][1,] # mean adj. Rsq
	lisi.rows[[length(lisi.rows)+1]] <- file.list[[j]][which(file.list[[j]]$group==i & file.list[[j]]$pa_strat=="(T_MLD1)"),][4,] # rank of mean of ranks of adj. Rsq
	# Close loop across files
	} 
	lisi.rows[[length(lisi.rows)+1]] <-rep(NA, 32) # Add filler with NAs
# Close loop across taxa
}
df.tax <- do.call(rbind, lisi.rows)
df.tax$sets <- rep(c(nms.sets, "NA"), length(vec.taxa))
df.tax[is.na(df.tax)] <- "________________"

# Save file
fln <- paste0(output.dir,"/Groups_(T_MLD1)_vs_thinning_(adjRsq).csv")
fln <- write.csv(df.tax,file=fln,row.names=F)


## Loop across groups Value == relative sd (of mean ad. Rsq)

# Preparatory
input.dir <- "../7_Generate_Absences/Output"
output.dir <- "./Output_vartest"
vec.taxa <- c("bacillariophyceae", "chlorophyta", "chrysophyceae", "cryptophyta", "cyanobacteria", "dinoflagellata", "euglenoidea","haptophyta", "raphidophyceae","total")
nms.sets <- c(
"gridded,gr_bg",
"gridded,tot_bg",
"mt_300,gr_bg", 
"mt_300,tot_bg", 
"mt_600,gr_bg", 
"mt_600,gr_bg")
file.list <- list(); for (d in 1:length(nms.sets)){file.list[[length(file.list)+1]]<- read.csv(paste0(output.dir,"/Groups_(",nms.sets[d],").csv"),header=T, sep=",")}
#
lisi.rows <- list ()
for (i in vec.taxa){
	# Loop across files: thinning strategies
	for (j in 1:length(file.list)){
	# lisi.rows[[length(lisi.rows)+1]] <- file.list[[j]][which(file.list[[j]]$group==i & file.list[[j]]$pa_strat=="(T_MLD1)"),][1,] # mean adj. Rsq
	lisi.rows[[length(lisi.rows)+1]] <- file.list[[j]][which(file.list[[j]]$group==i & file.list[[j]]$pa_strat=="(T_MLD1)"),][2,] # rel sd of mean adj. Rsq
	lisi.rows[[length(lisi.rows)+1]] <- file.list[[j]][which(file.list[[j]]$group==i & file.list[[j]]$pa_strat=="(T_MLD1)"),][4,] # rank of mean of ranks of adj. Rsq
	# Close loop across files
	} 
	lisi.rows[[length(lisi.rows)+1]] <-rep(NA, 32) # Add filler with NAs
# Close loop across taxa
}
df.tax <- do.call(rbind, lisi.rows)
df.tax$sets <- rep(c(rep(nms.sets, each = 2), "NA"), length(vec.taxa))
df.tax[is.na(df.tax)] <- "________________"

# Save file
fln <- paste0(output.dir,"/Groups_(T_MLD1)_vs_thinning_(SD).csv")
fln <- write.csv(df.tax,file=fln,row.names=F)






### PS ###########################################################################################################################



# Variables to test (desired order is according to ranking and correlations among variables from previous tests across taxa)
vc.vars<-c("x","y","T","N","P","logN","logP","Si","logSi","Sal","pCO2","Nstar","Pstar","Sistar","PAR","Wind.CCMP","monsd.Wind.CCMP", 
"MLPAR1","MLPAR2","MLD1","MLD2","logMLD1","logMLD2","Chl", "dMLD1_dt")


# Define functions:
# these function work with # in_data=data-frame, spp_name="obs", variables=vector with variables
# use for reproducible example: 
in_data <- taxa.ref.list[[4]]
spp_name <- "obs"
variables <- vc.vars

### 1. GLM
vartest0.glm_ap_dr <- function(in_data,spp_name,variables){ # in_data=data-frame, spp_name="obs", variables=vector with variables
ql.matrix <- data.frame(matrix(ncol=2,nrow=length(variables)))
names(ql.matrix) <- c("variable", "adjR2")
  for (i in 1:length(variables)){
  	var <- variables[i]
    tmp <- glm(in_data[,spp_name]>0 ~ in_data[,var] + I((in_data[,var])^2), na.action=na.omit, family=binomial)
    ql.matrix[(i),2] <- (1-(tmp$deviance/tmp$null.deviance)) # the D-squared  ?! I guess. Not adjusted D-squared!!!
    ql.matrix[(i),1] <- names(in_data)[var] 
    if (class(var)=="character"){ql.matrix[(i),1] <- var} # added in comparison to original formula (Achilleas)
    }
return(ql.matrix)}

# test
vartest0.glm_ap_dr(in_data, spp_name, variables)




