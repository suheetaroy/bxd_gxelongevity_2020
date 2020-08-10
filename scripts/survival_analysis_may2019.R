### Final figures
suppressMessages(library(survival))
suppressMessages(library(survminer))
suppressMessages(library(coxme))
suppressMessages(library(ggplot2))
suppressMessages(library(cowplot))
suppressMessages(library(ggrepel))
suppressMessages(library(plotly))
suppressMessages(library(lubridate))
suppressMessages(library(DT))
suppressMessages(library(openxlsx))
theme_set(theme_cowplot())
library(lme4)
library(lmerTest)
library(zoo)


strainReformat <- function(x){
  before <- as.character(x)
  wbxd <- which(grepl("BXD", before))
  after <- before
  ## add buffering zeros
  after[wbxd] <- gsub("BXD", "", after[wbxd])
  after[wbxd] <- sapply(after[wbxd], function(y){
    paste0("BXD", paste(rep("0", 3-nchar(gsub("[a-z]", "", y)) ), collapse = ""), y)
  })
  bxdlevs <- sort(unique(after[wbxd]))
  nonbxdlevs <- c("C57BL/6J", "DBA/2J", "B6D2F1", "D2B6F1")
  after <- factor(after, levels = c(nonbxdlevs, bxdlevs), ordered = T)
  after
}

### read and format the input data
dt <- read.table("./Data/AgingBXD_Longevitydata_Nov2018.natural.deaths.txt", sep = "\t", comment.char = "", header = T, stringsAsFactors = F)
if(any(!dt$Sex)){dt$Sex <- substr(as.character(dt$Sex),1,1)} # R is reading sex as a boolean
dt$DateOfDeath <- as.character(dt$DateOfDeath)
dtalive <- read.table("./Data/AgingBXD_Longevitydata_Nov2018.alive.txt", sep = "\t", comment.char = "", header = T)
dtalive <- dtalive[,colnames(dtalive) %in% colnames(dt)]
dtalive$Sex <- as.character(dtalive$Sex)
lastFollowupChar <- "11.15.2018"
lastFollowup <- mdy(lastFollowupChar)
# add missing columns to dtalive
miscol <- setdiff(colnames(dt), colnames(dtalive))
toadd <- as.data.frame(matrix(ncol = length(miscol), nrow = nrow(dtalive)))
colnames(toadd) <- miscol
dtalive <- cbind(dtalive, as.data.frame(toadd))
dtalive <- dtalive[,colnames(dt)]
dtalive$DateOfDeath <- lastFollowupChar
dtalive$CauseOfDeath <- "Alive"
dt <- rbind(dt, dtalive)
dt$DateOfBirth <- mdy(dt$DateOfBirth)
dt$DateOfDeath <- mdy(dt$DateOfDeath)
dt <- dt[dt$Sex == "F",] # remove all males
# dt$DateOfDeath[dt$CauseOfDeath == "Alive"] <- lastFollowup #### for the censored cases, this is not the real date of death
dt$DateDietStart <- mdy(dt$DateDietStart)
dt$AgeDietStart <- dt$DateDietStart - dt$DateOfBirth 
dt$AgeAtDeath <- dt$DateOfDeath - dt$DateOfBirth
dt$DietCode[dt$DietCode == "HF"] <- "HFD"
dt$Diet <- dt$DietCode
dt$StrainNameCurrent <- strainReformat(dt$StrainNameCurrent)
dt$Strain <- dt$StrainNameCurrent
dt$StrainDiet <- paste0(dt$StrainNameCurrent, "-", dt$Diet)
dt$Surv <- Surv(dt$AgeAtDeath, event = as.numeric(dt$CauseOfDeath != "Alive") )
### some exclusion criteria as per Pooja's remarks
# 1.	All males were excluded (as per discussion with Suheeta).
# 2.	DBA/2J-Gpmb were excluded.
dt <- dt[dt$StrainNameCurrent != "DBA/2J-Gpmb",]
# 3.	All mice that entered the colony in 2017 in the worksheet “Alive” were excluded”. Reason for exclusion: (i) There longevity so far was less than their corresponding dead mice (ii) New strains so we cannot know what their actual longevity would be.
toremove <- (dt$CauseOfDeath == "Alive") & (year(dt$DateOfBirth) == 2017)
dt <- dt[-which(toremove),]
# 4.	BXD 168 and 157 in the worksheet “Alive” were also excluded because they entered the colony between Oct-Dec 2016 and their lifespan so far is less than their corresponding mice that have already died.
toremove2 <- (dt$CauseOfDeath == "Alive") & (dt$StrainNameCurrent == "BXD168") & (dt$DateOfBirth > "2016-10-1") & (dt$DateOfBirth < "2016-12-31")
dt <- dt[-which(toremove2),]
toremove3 <- (dt$CauseOfDeath == "Alive") & (dt$StrainNameCurrent == "BXD157") & (dt$DateOfBirth > "2016-10-1") & (dt$DateOfBirth < "2016-12-31")
dt <- dt[-which(toremove3),]
# Also remove mice where the AgeAtDietStart > 365
dt <- dt[dt$AgeDietStart <= 365,]

# remove mice that lived less than 100 days
dt <- dt[dt$AgeAtDeath >= 100,]


# Color declaration and figure size variables
dietColors <- c("CD" = "#E30613", "HFD" = "#312783")
dietColorsStrata <- setNames(dietColors, paste0("Diet=", names(dietColors)))
figwidth <- figheight <- 4


# add weight data
bw <- read.xlsx("./Data/AgingBXD_BodyWts_OCT2018_forMBSfig.xlsx", startRow = 9)
head(bw)
bwBaseline <- read.xlsx("./Data/AgingBXD_baselineweights.xlsx", sheet = "Naturaldeath")
bw$Weight00 <- bwBaseline$Weight00Baseline[match(bw$EarTagNumberCurrent, bwBaseline$EarTagNumberCurrent)]

bwcolskeep <- grepl("Weight", colnames(bw))
bw$StrainNameCurrent <- strainReformat(bw$StrainNameCurrent)
dt <- cbind(dt, bw[match(dt$IDFemaleAgingColony, bw$IDFemaleAgingColony),bwcolskeep])
dt$AgeAtWeight00 <- dt$AgeDietStart
dt$DaysOnDietAtWeight00 <- 0

wv <- c("00", paste0("0", 1:9), 10:17)
wmv <- paste0("Weight", wv)
awv <- paste0("AgeAtWeight", wv)
adv <- paste0("DaysOnDietAtWeight", wv)
weightInt <- seq(from = 50, to = 1200, by = 25)

dt <- dt[-which(is.na(dt$IDFemaleAgingColony)),]

## get coef for slope calculation (body weight change - note that this is calculated on a monthly basis)
Coef <- function(Z) coef(lm(as.numeric(weight) ~ as.numeric(ageAtWeightMonths), as.data.frame(Z, stringsAsFactors = F)))    

trajectories <- lapply(dt$IDFemaleAgingColony, function(x){
  cur <- dt[dt$IDFemaleAgingColony == x,]
  
  out <- data.frame(weightID = wv, ageAtWeight = as.numeric(cur[,awv]),
                    daysOnDiet = as.numeric(cur[,adv]),
                    weight = as.numeric(cur[,wmv]))
  if(sum(!is.na(out$weight)) < 3){return(NULL)}
  out$IDFemaleAgingColony <- x
  out$StrainNameCurrent <- cur$StrainNameCurrent
  out$DietCode <- cur$DietCode
  out$Sex <- cur$Sex
  out$AgeAtDeath <- cur$AgeAtDeath
  out$CauseOfDeath <- cur$CauseOfDeath
  out$AgeAtDietStart <- cur$AgeDietStart
  out <- out[out$ageAtWeight < out$AgeAtDeath,]
  out <- out[!is.na(out$weight),]
  
  # add weight at age 0, and set to 0 grams
  zerogram <- out[1,]
  zerogram$weightID <- "-1"
  zerogram$ageAtWeight <- 0
  zerogram$daysOnDiet <- 0
  zerogram$weight <- 0
  out <- rbind(zerogram, out)
  
  # estimated weights at regular intervals
  weightIntCur <- weightInt[(weightInt < max(out$ageAtWeight) ) & (weightInt > min(out$ageAtWeight)) & (!weightInt %in% out$ageAtWeight)]
  real <- zoo(x = out$weight, order.by = out$ageAtWeight)
  cmbn <- c(real, zoo(rep(NA, length( weightIntCur )), order.by = weightIntCur))
  interpolated <- na.approx(cmbn)
  toadd <- interpolated[as.character(weightIntCur)]
  dftoadd <- out[rep(1, length(toadd)),]
  dftoadd$weightID <- paste0("Interpolated.", weightIntCur, ".day")
  dftoadd$ageAtWeight <- weightIntCur
  dftoadd$weight <- as.numeric(toadd)
  dftoadd$daysOnDiet <- dftoadd$ageAtWeight - dftoadd$AgeAtDietStart
  dftoadd$daysOnDiet[dftoadd$daysOnDiet<0] <- 0
  out <- rbind(out, dftoadd)
  rownames(out) <- NULL
  out <- out[order(out$ageAtWeight, decreasing = F),]
  out$ageAtWeightMonths <- out$ageAtWeight/30
  
  slopeDays <- 90
  # for each age at weight, calculate the slope for the next "slopeDays")
  rollingIntervals <- lapply(out$ageAtWeight, function(aaw){
    c(min(which(out$ageAtWeight >= (aaw - (slopeDays/2) ))), max(which(out$ageAtWeight <= (aaw + (slopeDays/2)))))
  })
  
  out$rollingIntervals <- rollingIntervals
  out$rollingAges <- lapply(out$rollingIntervals, function(ri){
    out[c(ri[1],ri[2]),"ageAtWeight"]
  })
  out$rollingSlope <- sapply(rollingIntervals, function(ri){
    Coef(out[ri[1]:ri[2],])[2]
  })
  out
})


# ggplot(trajectories[[167]], aes(x = ageAtWeight, y = weight)) + 
#   geom_point() +
#   geom_line(aes(y = rollingSlope*10)) +
#   geom_hline(yintercept = 0, color = "red")
# 

# paramByMouse <- lapply(trajectories, function(x){
#   maxWeight <- max(x$weight)
#   maxWeightAge <- x$ageAtWeight[which.max(x$weight)]
#   maxWeightDaysOnDiet <- x$daysOnDiet[which.max(x$weight)]
#   ### to complete
#   
# })

trajectories <- do.call(rbind, trajectories)
trajectories$Interpolated <- grepl("Interpolated", trajectories$weightID)
ma <- mean(dt$AgeAtDeath)
trajectories <- trajectories[order(trajectories$AgeAtDeath, decreasing = F),]








### plot the weight evolution per strain
trajectories.strain <- split(trajectories, trajectories$StrainNameCurrent)

trajetories.strain.plots <- lapply(trajectories.strain, function(x){
  if(nrow(x) == 0){return(NULL)}
  ds <- x[match(unique(x$IDFemaleAgingColony), x$IDFemaleAgingColony),]
  de <- x[order(x$ageAtWeight, decreasing = T),]
  de <- de[match(unique(de$IDFemaleAgingColony), de$IDFemaleAgingColony),]
  ggplot(x, aes(x=ageAtWeight, y = weight, group = IDFemaleAgingColony, color = AgeAtDeath)) + 
    geom_line() + facet_wrap(~DietCode) +
    scale_color_gradient2(midpoint = as.numeric(ma), low = "red", mid = "grey", high = "green" ) +
    geom_point(aes(shape = Interpolated)) +
    geom_segment(data = ds,aes(x=AgeAtDietStart, xend = ageAtWeight, y = 0, yend = weight), color = "grey", lty = 2 , size = 0.3) +
    geom_segment(data = de,aes(x=AgeAtDeath, xend = ageAtWeight, y = 0, yend = weight), color = "grey", lty = 2 , size = 0.3) +
    xlim(0, 1100) +
    ylim(0, max(trajectories$weight)) +
    ggtitle(paste0("Strain = ", x$StrainNameCurrent[1]))
})
# trajetories.strain.plots$`DBA/2J`
# trajetories.strain.plots$BXD79
# 
# trajetories.strain.plots$BXD65b
pdf("./Plots/trajectories.per.strain.pdf", width = 8, height = 6)

lapply(trajetories.strain.plots, function(x){
  print(x)
  
})

dev.off()

ds <- trajectories[match(unique(trajectories$IDFemaleAgingColony), trajectories$IDFemaleAgingColony),]
de <- trajectories[order(trajectories$ageAtWeight, decreasing = T),]
de <- de[match(unique(de$IDFemaleAgingColony), de$IDFemaleAgingColony),]

facetedByStrain <- ggplot(trajectories, aes(x=ageAtWeight, y = weight, group = IDFemaleAgingColony, color = AgeAtDeath, lty = DietCode)) + 
  geom_line(size = 0.5) + facet_wrap(~StrainNameCurrent) +
  scale_color_gradient2(midpoint = as.numeric(ma), low = "red", mid = "grey", high = "green" ) +
  # geom_segment(data = ds,aes(x=AgeAtDietStart, xend = ageAtWeight, y = 0, yend = weight), color = "grey", lty = 2 , size = 0.3) +
  # geom_segment(data = de,aes(x=AgeAtDeath, xend = ageAtWeight, y = 0, yend = weight), color = "grey", lty = 2 , size = 0.3) +
  xlim(0, 1100) +
  ylim(0, max(trajectories$weight))

facetedByStrain + ggsave("./Plots/trajectories.per.strain.faceted.pdf", width = 15, height = 15)


#First, a look at the overall diet effect.
fit <- survfit(Surv ~ Diet, data= dt)
sfit <- summary(fit)
# coxnormal <- coxph(Surv ~ Diet, dt)
cme <- coxme(Surv ~ Diet + (1|Strain), dt)
scme <- summary(cme)
beta <- cme$coefficients
nvar <- length(beta)
nfrail<- nrow(cme$var) - nvar
se <- sqrt(diag(cme$var)[nfrail+1:nvar])
pval <- 1 - pchisq((beta / se)^2, 1)
md <- summary(fit)$table[,"median"]
names(md) <- gsub("Diet=", "", names(md))
ptext <- paste0("p = ", signif(pval,2) ,"\nMedian: ", md["CD"], " (CD), ", md["HFD"], " (HFD)")

ggsurvplot(fit, 
           risk.table = F, 
           surv.median.line = "h",
           pval = ptext,
           palette = dietColorsStrata,
           pval.size = 4) +
  ggsave("./Plots/overall.diet.effect.surv.filtered.pdf", width = figwidth, height = figheight)

sink("./Plots/overall.diet.effect.surv.filtered.txt")
print(cme)
sink()

## test proportional hazard assumption, cannot find a way for coxme, so used the frailty method within coxph
# coxph <- coxph(Surv ~ Diet + frailty(Strain), dt)
coxph <- coxph(Surv ~ Diet, dt) # not taking strain into account....
zp <- cox.zph(coxph)

ptextcumhaz <- paste0("cox.zph p = ", signif(zp$table[,"p"],2))
ggsurvplot(fit,
           fun = "cumhaz",
           risk.table = F,
           pval = ptextcumhaz,
           pval.size = 4,
           palette = dietColorsStrata,
           pval.coord = c(0,0.5)) +
  ggsave("./Plots/overall.diet.effect.hazard.filtered.pdf", width = figwidth, height = figheight)

## BL6 and dba survival
b6dba <- dt[dt$Strain %in% c("DBA/2J", "C57BL/6J"),]
tb <- table(b6dba$StrainDiet)
fit <- survfit(Surv ~ Strain + Diet, data= b6dba)
md <- summary(fit)$table[,"median"]
caption <- paste(paste0(names(md), md), collapse = "\n")
caption <- gsub("Strain=|Diet=|,", "", caption)
ggsurvplot(fit, linetype = "Strain", color = "Diet", palette = dietColors, pval = caption) +
  ggsave("./Plots/b6dba.surv.filtered.pdf", width = figwidth, height = figheight)

# Plot differences in longevity between two diets at single strain level
fit <- survfit(Surv ~ Strain + Diet, data = dt)
sfit <- summary(fit)
median.table <- as.data.frame(sfit$table)
median.table <- median.table[!is.na(median.table$median),]
median.table <- median.table[median.table$records > 4,]
median.table <- median.table[median.table$events > 4,]
median.table$ID <- rownames(median.table)
spl <- strsplit(median.table$ID, ", ", fixed = T)
spl <- lapply(spl, function(x){
  x <- gsub("Strain=", "", x)
  x <- gsub("Diet=", "", x)
  x <- gsub(" ", "", x)
  x
})
spl <- do.call(rbind, spl)
median.table$Strain <- spl[,1]
median.table$Diet <- spl[,2]

#remove strains that appear once
tb <- table(median.table$Strain)
# rm <- names(tb)[tb < 2]
# median.table <- median.table[!median.table$Strain %in% rm,]
diff <- unlist(lapply(split(median.table, median.table$Strain),function(x){
  x$median[x$Diet == "CD"] - x$median[x$Diet == "HFD"]
}))

dt$Strain <- factor(dt$StrainNameCurrent, levels = names(sort(diff)), ordered = T)
dt$diff <- diff[as.character(dt$StrainNameCurrent)]
dtplot <- droplevels(dt[!is.na(dt$Strain) ,])

ggplot(dtplot, aes(x = Strain, y = AgeAtDeath, color = Diet)) +
  # geom_dotplot(aes(group = interaction(Strain, Diet)),binaxis = "y", method = "histodot", binwidth = 0.3,  stackdir="center")
  # geom_boxplot(outlier.size = -1) +
  geom_point(position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.1), size = 0.5) +
  geom_crossbar(data=droplevels(median.table[median.table$Strain %in% dtplot$Strain,]),aes(ymin=median, ymax=median,y=median), width = 0.5) +
  # stat_summary(data= droplevels(median.table[median.table$Strain %in% dtplot$Strain,]), aes(y = median, yend = median, xend = Strain, group = 1), geom = "segment") +
  theme_cowplot() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1), legend.position = "none") +
  xlab("") + ylab("Age at death (days)") + 
  scale_y_continuous(breaks = seq(signif(min(diff), -2) - 100, max(as.numeric(dt$AgeAtDeath)), 100)) +
  scale_color_manual(values = c("#E30613", "#312783")) +
  geom_col(data = dtplot[match(names(diff), dtplot$Strain),], aes(y= diff), color = "#FFFFFF00", fill = "grey", position= "identity", alpha = 0.5) +
  ggsave("./Plots/death.boxplot.sorted.by.hfdresponse.pdf", width = figwidth*3, height = figheight, useDingbats = F)

write.table(median.table, file = './Plots/death.boxplot.sorted.by.hfdresponse.medians.txt')
write.table(diff, file = './Plots/death.boxplot.sorted.by.hfdresponse.median.diff.txt')
write.table(dtplot, file = './Plots/death.boxplot.sorted.by.hfdresponse.dtplot.txt', sep = "\t")
write.table(do.call(rbind, trajectories), file = './Data/trajectories.txt', sep = "\t")



# stats for weight and longevity along time (weightInt)
weightsForModel <- weightInt[weightInt<950]
lmerWeight <- lapply(weightsForModel, function(wi){
  sub <- trajectories[trajectories$ageAtWeight == wi,]
  model <- lmer(weight ~ DietCode + daysOnDiet  + (1|StrainNameCurrent), data = sub)
  print(wi)
  
  modelCD <- lmer(weight ~ daysOnDiet + (1|StrainNameCurrent), data = sub[sub$DietCode == "CD",])

  modelHFD <- lmer(weight ~ daysOnDiet + (1|StrainNameCurrent), data = sub[sub$DietCode == "HFD",])

  sub$Surv <- Surv(sub$AgeAtDeath, rep(1, nrow(sub)))
  cmodel <- coxph(Surv ~  DietCode + daysOnDiet + weight + frailty(StrainNameCurrent), data = sub)
  cmodel$weightMin <- min(sub$weight, na.rm = T)
  cmodel$weightMax <- max(sub$weight, na.rm = T)
  cmodel$rollingSlopeMin <- NA
  cmodel$rollingSlopeMax <- NA
  cmodel$data <- sub
  
  cmodelCD <- coxph(Surv ~  daysOnDiet + weight + frailty(StrainNameCurrent), data = sub[sub$DietCode == "CD",])
  cmodelCD$weightMin <- min( sub[sub$DietCode == "CD",]$weight, na.rm = T)
  cmodelCD$weightMax <- max( sub[sub$DietCode == "CD",]$weight, na.rm = T)
  cmodelCD$rollingSlopeMin <- NA
  cmodelCD$rollingSlopeMax <- NA
  cmodelCD$data <- sub[sub$DietCode == "CD",]
  
  
  cmodelHFD <- coxph(Surv ~  daysOnDiet + weight+ frailty(StrainNameCurrent), data = sub[sub$DietCode == "HFD",])
  cmodelHFD$weightMin <- min( sub[sub$DietCode == "HFD",]$weight, na.rm = T)
  cmodelHFD$weightMax <- max( sub[sub$DietCode == "HFD",]$weight, na.rm = T)
  cmodelHFD$rollingSlopeMin <- NA
  cmodelHFD$rollingSlopeMax <- NA
  cmodelHFD$data <- sub[sub$DietCode == "HFD",]
  
  cmodelslope <- coxph(Surv ~  DietCode + weight + daysOnDiet + rollingSlope + frailty(StrainNameCurrent), data = sub)
  cmodelslope$weightMin <- min( sub$weight, na.rm = T)
  cmodelslope$weightMax <- max( sub$weight, na.rm = T)
  cmodelslope$rollingSlopeMin <- min( sub$rollingSlope, na.rm = T)
  cmodelslope$rollingSlopeMax <- max( sub$rollingSlope, na.rm = T)
  cmodelslope$data <- sub
  
  cmodelCDslope <- coxph(Surv ~  daysOnDiet + weight + rollingSlope + frailty(StrainNameCurrent), data = sub[sub$DietCode == "CD",])
  cmodelCDslope$weightMin <- min( sub[sub$DietCode == "CD",]$weight, na.rm = T)
  cmodelCDslope$weightMax <- max( sub[sub$DietCode == "CD",]$weight, na.rm = T)
  cmodelCDslope$rollingSlopeMin <- min( sub[sub$DietCode == "CD",]$rollingSlope, na.rm = T)
  cmodelCDslope$rollingSlopeMax <- max( sub[sub$DietCode == "CD",]$rollingSlope, na.rm = T)
  cmodelCDslope$data <- sub[sub$DietCode == "CD",]
  
  cmodelHFDslope <- coxph(Surv ~  daysOnDiet + weight + rollingSlope + frailty(StrainNameCurrent), data = sub[sub$DietCode == "HFD",])
  cmodelHFDslope$weightMin <- min( sub[sub$DietCode == "HFD",]$weight, na.rm = T)
  cmodelHFDslope$weightMax <- max( sub[sub$DietCode == "HFD",]$weight, na.rm = T)
  cmodelHFDslope$rollingSlopeMin <- min( sub[sub$DietCode == "HFD",]$rollingSlope, na.rm = T)
  cmodelHFDslope$rollingSlopeMax <- max( sub[sub$DietCode == "HFD",]$rollingSlope, na.rm = T)
  cmodelHFDslope$data <- sub[sub$DietCode == "HFD",]
  
  cmodel$data <- sub # stored the data in cmodel
  return(list(model=model, 
              modelCD = modelCD, 
              modelHFD = modelHFD, 
              cmodel = cmodel, 
              cmodelCD = cmodelCD, 
              cmodelHFD = cmodelHFD, 
              cmodelslope = cmodelslope, 
              cmodelCDslope = cmodelCDslope, 
              cmodelHFDslope = cmodelHFDslope ))
})

names(lmerWeight)<-weightsForModel
summary(lmerWeight$`50`$cmodel)$coef

lmerWeightReshaped <- do.call(rbind, mapply(x=lmerWeight,y=names(lmerWeight), function(x,y){
  do.call(rbind, mapply(z=x, w = names(x), function(z,w){
    
    out <- reshape::melt(summary(z)$coef)
    
    out$model <- w
    out$ageAtWeight <- as.numeric(y)
    
    if(grepl("cmodel", w)){
      
      out$weightMin = z$weightMin
      out$weightMax = z$weightMax
      out$rollingSlopeMin = z$rollingSlopeMin
      out$rollingSlopeMax = z$rollingSlopeMax
    }else{
      
      out$weightMin = NA
      out$weightMax = NA
      out$rollingSlopeMin = NA
      out$rollingSlopeMax = NA
      
    }
    
    out
  }, SIMPLIFY = FALSE))
}, SIMPLIFY = FALSE))

lmerWeightReshaped$modelType <- ifelse(grepl("cmodel", lmerWeightReshaped$model), "coxph", "lmer")
lmerWeightReshaped$modelWithSlope <- ifelse(grepl("slope", lmerWeightReshaped$model), "with slope", "no slope")

lmerWeightReshaped$modelDiet <- gsub("slope", "",gsub("(cmodel)|(model)", "",lmerWeightReshaped$model))
lmerWeightReshaped$modelDiet[lmerWeightReshaped$modelDiet == ""] <- "Both"

lmerWeightReshaped$measureUnit <- ifelse(lmerWeightReshaped$X1 %in% "DietCodeHFD", "categorical", "continuous")


p <- ggplot(lmerWeightReshaped[(lmerWeightReshaped$X2 == "p")&(!grepl("frail", lmerWeightReshaped$X1))&(grepl("cmodel", lmerWeightReshaped$model)),], 
       aes(x = ageAtWeight, y = -log10(value), group = X1, color = X1)) +
  geom_line() +
  geom_hline(yintercept = -log10(0.05), color = 'grey') +
  facet_grid(modelWithSlope~modelDiet)

p + ggsave(filename = "./Plots/modelPvalues.with.withoutslope.pdf", width = 10, height = 6)

ph <- ggplot(lmerWeightReshaped[(lmerWeightReshaped$X2 == "coef")&(!grepl("frail", lmerWeightReshaped$X1))&(grepl("cmodel", lmerWeightReshaped$model)),], 
            aes(x = ageAtWeight, y = exp(value), group = X1, color = X1)) +
  geom_line() +
  geom_hline(yintercept = 1, color = 'grey') +
  facet_grid(modelWithSlope+measureUnit~modelDiet, scales = "free_y") +
  ylab("Hazard ratio")
ph + ggsave(filename = "./Plots/modelHazard.with.withoutslope.pdf", width = 10, height = 6)


lmerWeightReshaped.weight <- lmerWeightReshaped[lmerWeightReshaped$modelType == "lmer",]
pw <- ggplot(lmerWeightReshaped.weight[(lmerWeightReshaped.weight$X2 == "Pr(>|t|)")&(!grepl("Interc", lmerWeightReshaped.weight$X1))&(grepl("model", lmerWeightReshaped.weight$model)),], 
            aes(x = ageAtWeight, y = -log10(value), group = X1, color = X1)) +
  geom_line() +
  geom_hline(yintercept = -log10(0.05), color = 'grey') +
  facet_grid(~modelDiet) +
  ylab("-log10(p-value)")

pwe <- ggplot(lmerWeightReshaped.weight[(lmerWeightReshaped.weight$X2 == "Estimate")&(!grepl("Interc", lmerWeightReshaped.weight$X1))&(grepl("model", lmerWeightReshaped.weight$model)),], 
             aes(x = ageAtWeight, y = value, group = X1, color = X1)) +
  geom_line() +
  geom_hline(yintercept = 0, color = 'grey') +
  facet_grid(~modelDiet) +
  ylab("Estimate (coefficient)")

plot_grid(pw, pwe, nrow = 2) + ggsave(filename = "./Plots/weightModelPvaluesAndCoef.pdf", width = 10, height = 6)

# ggplot(lmerWeightReshaped[(lmerWeightReshaped$X2 == "coef")&(!grepl("frail", lmerWeightReshaped$X1))&(grepl("cmodel", lmerWeightReshaped$model)),], 
#        aes(x = ageAtWeight, y = exp(-value), group = X1, color = X1)) +
#   geom_line() +
#   facet_wrap(~model)


strainMean <- plyr::ddply(dt, c("DietCode", "StrainNameCurrent"), summarise, meanAgeAtDeath = mean(AgeAtDeath))
strainMean$StrainDiet <- paste(strainMean$StrainNameCurrent, strainMean$DietCode, sep = "-")
trajectories$StrainDiet <- paste(trajectories$StrainNameCurrent, trajectories$DietCode, sep = "-")
trajectories$meanAgeAtDeath <- as.numeric(strainMean$meanAgeAtDeath[match(trajectories$StrainDiet, strainMean$StrainDiet)])


pall <- ggplot(trajectories[grepl("Interpo", trajectories$weightID),], aes(x=ageAtWeight, y = weight, group = IDFemaleAgingColony, color = AgeAtDeath)) + 
  geom_line() + facet_wrap(~DietCode) +
  scale_color_gradient2(midpoint = as.numeric(ma), low = "red", mid = "grey", high = "green" )
pall


weightsToPlot <- c(100,200, 300,400,500, 600, 700)

a <- ggplot(trajectories[trajectories$ageAtWeight %in% weightsToPlot,], aes(x = weight, y = AgeAtDeath, size = weight, fill=weight)) +
  geom_point(pch = 21, alpha = 0.5) +
  geom_smooth(method = "lm", color = "red", size = 1) +
  geom_hline(aes(yintercept = ageAtWeight), color = "grey") +
  facet_grid(DietCode~ageAtWeight)+
  scale_fill_gradient2(low = "green", mid = "grey" , high = "red", midpoint = 50) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  xlab("Body weight") +
  ylab("Age at death (days)")

b <- ggplot(trajectories[trajectories$ageAtWeight %in% weightsToPlot,], aes(x = rollingSlope, y = AgeAtDeath, size = weight, fill=weight)) +
  geom_point(pch = 21, alpha = 0.5) +
  geom_vline(xintercept = 0, color = "red") +
  geom_hline(aes(yintercept = ageAtWeight), color = "darkgrey") +
  geom_smooth(method = "lm", color = "red", size = 1) +
  facet_grid(DietCode~ageAtWeight) +
  scale_size_continuous(range = c(0.1,6)) +
  scale_fill_gradient2(low = "green", mid = "grey" , high = "red", midpoint = 50) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  xlab("Body weight slope (beta)")+
  ylab("Age at death (days)")

coxsubset <- lmerWeightReshaped[(lmerWeightReshaped$X2 == "p")&
                                  (!grepl("frail", lmerWeightReshaped$X1))&
                                  (grepl("cmodel", lmerWeightReshaped$model))&
                                  (lmerWeightReshaped$modelDiet != "Both")&
                                  (lmerWeightReshaped$modelWithSlope == "with slope"),]
coxsubsetbeta <- lmerWeightReshaped[(lmerWeightReshaped$X2 == "coef")&
                                  (!grepl("frail", lmerWeightReshaped$X1))&
                                  (grepl("cmodel", lmerWeightReshaped$model))&
                                  (lmerWeightReshaped$modelDiet != "Both")&
                                  (lmerWeightReshaped$modelWithSlope == "with slope"),]
coxsubsetp <- lmerWeightReshaped[(lmerWeightReshaped$X2 == "p")&
                                      (!grepl("frail", lmerWeightReshaped$X1))&
                                      (grepl("cmodel", lmerWeightReshaped$model))&
                                      (lmerWeightReshaped$modelDiet != "Both")&
                                      (lmerWeightReshaped$modelWithSlope == "with slope"),]
coxsubsetbeta$p <- coxsubsetp$value

c <- ggplot(coxsubset[coxsubset$ageAtWeight %in% weightsToPlot,], 
            aes(y = -log10(value), fill = X1, x = X1)) +
  geom_col() +
  geom_hline(aes(yintercept = -log10(0.05)), color = 'grey') +
  facet_grid(modelDiet ~ageAtWeight) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  xlab("Cox model coefficient")+
  ylab("-log10(p-value)")
c

# coxsubsetbeta <- coxsubsetbeta[coxsubsetbeta$X1 != "rollingSlope",]

d <- ggplot(coxsubsetbeta[coxsubsetbeta$ageAtWeight %in% weightsToPlot,], 
           aes(y = exp(value), fill = X1, x = X1)) +
  geom_col() +
  geom_hline(aes(yintercept = exp(0) ), color = 'grey') +
  facet_grid(modelDiet ~ageAtWeight) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  xlab("Cox model coefficient")+
  ylab("coef")
d
plot_grid(c,d, nrow = 2)

plot_grid(a,b,c, nrow = 3, align = "hv") + ggsave(filename = "./Plots/weight.longevity.coxph.slopes.pdf", width = 10, height = 15)

pall <- ggplot(trajectories[,], aes(x=ageAtWeight, y = weight, group = IDFemaleAgingColony, color = AgeAtDeath)) + 
  geom_line() + facet_wrap(~DietCode) +
  scale_color_gradient2(midpoint = as.numeric(ma), low = "red", mid = "grey", high = "green", name = "Age at death") +
  xlim(0, 1100)+
  xlab("Age (days)") +
  ylab("Weight (g)")
pall


pallslope <- ggplot(trajectories[,], aes(x=ageAtWeight, y = rollingSlope, group = IDFemaleAgingColony, color = AgeAtDeath)) + 
  geom_line() + facet_wrap(~DietCode) +
  scale_color_gradient2(midpoint = as.numeric(ma), low = "red", mid = "grey", high = "green", name = "Age at death") +
  geom_hline(yintercept = 0, color = "black") +
  # geom_smooth(group = 1) +
  xlim(0, 1100) +
  xlab("Age (days)") +
  ylab("Weight slope\n(g/month)")
pallslope

call <- ggplot(coxsubset[(coxsubset$ageAtWeight > 50) & (coxsubset$X1 != "daysOnDiet"), ], aes(x = ageAtWeight, y = -log10(value), group = X1, color = X1)) +
  geom_line() +
  geom_point() +
  geom_hline(yintercept = -log10(0.05), color = "grey") +
  geom_hline(yintercept = -log10(0.05 / length(unique(coxsubset$ageAtWeight[(coxsubset$ageAtWeight > 50)])))) +
  scale_color_discrete(name = "Model parameter") +
  facet_wrap(~modelDiet)  +
  xlim(0, 1100)+
  xlab("Age (days)") +
  ylab(expression(-log[10](p-value)))

callbeta <- ggplot(coxsubsetbeta[(coxsubsetbeta$ageAtWeight > 50) & (coxsubsetbeta$X1 != "daysOnDiet"),], aes(x = ageAtWeight, y = exp(value), group = X1, color = X1)) +
  geom_line() +
  geom_point() +
  geom_hline(yintercept = exp(0.0), color = "grey") +
  scale_color_discrete(name = "Model parameter") +
  facet_wrap(~modelDiet)  +
  xlim(0, 1100) +
  xlab("Age (days)") +
  ylab("Hazard Ratio")


plot_grid(pall, pallslope, call, callbeta, nrow = 4, align = "hv") + 
  ggsave(filename = "./Plots/weight.longevity.coxph.model.parameters.pdf", width = 10, height = 10)


## now some effect size plots as a function of different variables
# for the weight
weightRange <- range(coxsubsetbeta$weightMax -  coxsubsetbeta$weightMin)
weightHalfRange <- round((weightRange[2]-weightRange[1])/2)
weightsEffectVector <- seq(-weightHalfRange,weightHalfRange,0.1)

coxsubsetbetaWeight <- lapply(weightsEffectVector, function(x){
  out <- coxsubsetbeta[coxsubsetbeta$X1 == "weight",]
  out$weightChange <- x
  out$hazardRatio <- exp(out$value)^x
  out
})
coxsubsetbetaWeight <- do.call(rbind, coxsubsetbetaWeight)
coxsubsetbetaWeight$weightDiff <- coxsubsetbetaWeight$weightMax - coxsubsetbetaWeight$weightMin
coxsubsetbetaWeight$weightMargin <- coxsubsetbetaWeight$weightDiff / 2

coxsubsetbetaWeight.plot <- coxsubsetbetaWeight[(coxsubsetbetaWeight$ageAtWeight < 1000) & (coxsubsetbetaWeight$ageAtWeight >49),]
coxsubsetbetaWeight.plot$withinRange <- abs(coxsubsetbetaWeight.plot$weightChange) < coxsubsetbetaWeight.plot$weightMargin

# add text annotations
coxsubsetbetaWeight.plot <- do.call(rbind, lapply(split(coxsubsetbetaWeight.plot, paste0(coxsubsetbetaWeight.plot$modelDiet, coxsubsetbetaWeight.plot$ageAtWeight)), function(x){
  first <- min(which(x$withinRange))
  last <- max(which(x$withinRange))
  x$label <- NA
  x$label[c(first, last)] <- round(x$hazardRatio[c(first, last)], 2)
  x
}))
coxsubsetbetaWeight.plot$label50 <- coxsubsetbetaWeight.plot$label
coxsubsetbetaWeight.plot$label50[(coxsubsetbetaWeight.plot$ageAtWeight %% 50) != 0] <- NA



rollingSlopeRange <- range(coxsubsetbeta$rollingSlopeMax -  coxsubsetbeta$rollingSlopeMin)
rollingSlopeHalfRange <- round((rollingSlopeRange[2]-rollingSlopeRange[1])/2)

# weightsSlopeEffectVector <- seq(-5,5,0.1)
weightsSlopeEffectVector <- seq(-rollingSlopeHalfRange,rollingSlopeHalfRange, 0.1)


  coxsubsetbetaSlope <- lapply(weightsSlopeEffectVector, function(x){
    out <- coxsubsetbeta[coxsubsetbeta$X1 == "rollingSlope",]
    out$slopeChange <- x
    out$hazardRatio <- exp(out$value)^x
    out
  })
  
  coxsubsetbetaSlope <- do.call(rbind, coxsubsetbetaSlope)
  coxsubsetbetaSlope$slopeDiff <- coxsubsetbetaSlope$rollingSlopeMax - coxsubsetbetaSlope$rollingSlopeMin
  coxsubsetbetaSlope$slopeMargin <- coxsubsetbetaSlope$slopeDiff / 2
  
  weightSlopeEffect.plot <- coxsubsetbetaSlope[(coxsubsetbetaSlope$ageAtWeight < 1000)&(coxsubsetbetaSlope$ageAtWeight >50),]
  weightSlopeEffect.plot$withinRange <- abs(weightSlopeEffect.plot$slopeChange) < weightSlopeEffect.plot$slopeMargin
  # add text annotations
  weightSlopeEffect.plot <- do.call(rbind, lapply(split(weightSlopeEffect.plot, paste0(weightSlopeEffect.plot$modelDiet, weightSlopeEffect.plot$ageAtWeight)), function(x){
    first <- min(which(x$withinRange), na.rm= T)
    last <- max(which(x$withinRange), na.rm= T)
    x$label <- NA
    x$label[c(first, last)] <- round(x$hazardRatio[c(first, last)], 2)
    x
  }))
  weightSlopeEffect.plot$label50 <- weightSlopeEffect.plot$label
  weightSlopeEffect.plot$label50[(weightSlopeEffect.plot$ageAtWeight %% 50) != 0] <- NA

  # get the ranges of the two plots to harmonize the colors
  
  prange <- sort(-log10(range(c(coxsubsetbetaWeight.plot$p, weightSlopeEffect.plot$p))), decreasing = F)
  prange[1] <- floor(prange[1])
  prange[2] <- ceiling(prange[2])
  hrange <- sort(log2(range(c(coxsubsetbetaWeight.plot$hazardRatio, weightSlopeEffect.plot$hazardRatio))), decreasing = F)
  

  weightEffect <- ggplot(coxsubsetbetaWeight.plot, 
                         aes(x = ageAtWeight, fill = log2(hazardRatio), y = weightChange, label = label50, color = -log10(p))) +
    geom_tile(color = NA) +
    facet_wrap(~modelDiet)  +
    xlim(0, 1100) +
    scale_fill_gradient2(midpoint = 0, low = "green", high = "red", mid = "grey", name =  expression(log[2](HR)), limits = hrange) +
    scale_color_gradient(low = "lightgrey", high = "black", limits = prange, name = expression(-log[10](p-value)) ) +
    geom_rug(data = coxsubsetbetaWeight.plot[coxsubsetbetaWeight.plot$weightChange == 0,] , 
             size = 2, sides ="tb",
             aes(), y = 1) +
    geom_point(data = coxsubsetbetaWeight.plot[!is.na(coxsubsetbetaWeight.plot$label50),] , shape = 21, color = "black") +
    geom_label_repel(color = "black") +
    geom_hline(yintercept = 0) +
    xlab("Age (days)")+
    ylab("Weight difference\n(grams)")
  
  weightEffect
  
  
   
weightSlopeEffect <- ggplot(weightSlopeEffect.plot, 
                           aes(x = ageAtWeight, fill = log2(hazardRatio), y = slopeChange, label = label50, color = -log10(p))) +
  geom_tile(color = NA) +
  # geom_point() +
  # geom_hline(yintercept = exp(0.0), color = "grey") +
  facet_wrap(~modelDiet)  +
  xlim(0, 1100) +
  scale_fill_gradient2(midpoint = 0, low = "green", high = "red", mid = "grey", name =  expression(log[2](HR)), limits = hrange) +
  # scale_fill_gradientn(colours = c("grey", "blue", "yellow", "green", "red"), name = "Hazard Ratio") +
  # scale_fill_viridis_c(name = "Hazard Ratio", option = 'viridis') +
  scale_color_gradient(low = "lightgrey", high = "black", limits = prange, name = expression(-log[10](p-value)) ) +
  geom_rug(data = weightSlopeEffect.plot[weightSlopeEffect.plot$slopeChange == 0,] , 
           size = 2, sides ="tb",
           aes(), y = 1) +
  geom_point(data = weightSlopeEffect.plot[!is.na(weightSlopeEffect.plot$label50),] , shape = 21, color = "black") +
  geom_label_repel(color = "black") +
  geom_hline(yintercept = 0) +
  xlab("Age (days)")+
  ylab("Weight slope change\n(grams/month)") 
  
weightSlopeEffect


weightRangePlot <- ggplot(coxsubsetbetaWeight.plot[coxsubsetbetaWeight.plot$weightChange ==0,], aes(x= ageAtWeight)) +
  geom_linerange(aes(ymin = weightMin, ymax = weightMax)) +
  geom_hline(yintercept = 0) +
  xlim(0, 1100) +
  facet_wrap(~modelDiet)+
  xlab("Age (days)")+
  ylab("Range of weights\n(grams)") 
weightRangePlot

slopeRangePlot <- ggplot(weightSlopeEffect.plot[weightSlopeEffect.plot$slopeChange ==0,], aes(x= ageAtWeight)) +
  geom_linerange(aes(ymin = rollingSlopeMin, ymax = rollingSlopeMax)) +
  geom_hline(yintercept = 0) +
  xlim(0, 1100) +
  facet_wrap(~modelDiet) +
  xlab("Age (days)")+
  ylab("Range of weight slopes\n(grams/month)") 
slopeRangePlot



plot_grid(weightEffect, weightSlopeEffect, align = "hv", nrow = 2) +
  ggsave("./Plots/effect.hazard.trajectory.pdf", width = 10, height = 7)

plot_grid(pall, pallslope, call, callbeta,weightEffect + guides(fill = "none"),weightSlopeEffect + guides(color = "none"), nrow = 6, align = "hv", labels = LETTERS[1:6],
          rel_heights = c(1,1,0.75,0.75,1.5,1.5))+
  ggsave("./Plots/weight.longevity.coxph.model.parameters.effect.hazard.trajectory.pdf", width = 10, height = 16)


# plot the cox model data.
weightsToPlot_scatter <- c(100,200,300,400,500,600,700,800,900,1000)
weightVsAge <- ggplot(trajectories[trajectories$ageAtWeight %in% weightsToPlot_scatter,], aes(x = weight, y = AgeAtDeath, color = DietCode)) +
  geom_point(alpha = 0.5) + facet_wrap(~ageAtWeight, ncol = 5) +
  scale_color_manual(values = dietColors) +
  geom_smooth(method = "lm", fill = 'white')

rollingSlopeVsAge <- ggplot(trajectories[trajectories$ageAtWeight %in% weightsToPlot_scatter,], aes(x = rollingSlope, y = AgeAtDeath, color = DietCode)) +
  geom_point(alpha = 0.5) + facet_wrap(~ageAtWeight, ncol = 5) +
  scale_color_manual(values = dietColors) +
  geom_smooth(method = "lm", fill = 'white')

weightVsRollingSlope <- ggplot(trajectories[trajectories$ageAtWeight %in% weightsToPlot_scatter,], aes(x = weight, y = rollingSlope, color = DietCode)) +
  geom_point(alpha = 0.5) + facet_wrap(~ageAtWeight, ncol = 5) +
  scale_color_manual(values = dietColors) +
  geom_smooth(method = "lm", fill = 'white')
plot_grid(weightVsAge,rollingSlopeVsAge,weightVsRollingSlope, nrow = 3) +
  ggsave("./Plots/scatterPlotsbyAge.pdf", width = 16, height = 16)

############ do something similar to both diets ####################### 

coxsubsetbetaBoth <- lmerWeightReshaped[(lmerWeightReshaped$X2 == "coef")&
                                      (!grepl("frail", lmerWeightReshaped$X1))&
                                      (grepl("cmodel", lmerWeightReshaped$model))&
                                      (lmerWeightReshaped$modelDiet == "Both")&
                                      (lmerWeightReshaped$modelWithSlope == "with slope"),]
coxsubsetpBoth <- lmerWeightReshaped[(lmerWeightReshaped$X2 == "p")&
                                   (!grepl("frail", lmerWeightReshaped$X1))&
                                   (grepl("cmodel", lmerWeightReshaped$model))&
                                   (lmerWeightReshaped$modelDiet == "Both")&
                                   (lmerWeightReshaped$modelWithSlope == "with slope"),]
coxsubsetbetaBoth$p <- coxsubsetpBoth$value

weightRange <- range(coxsubsetbetaBoth$weightMax -  coxsubsetbetaBoth$weightMin)
weightHalfRange <- round((weightRange[2]-weightRange[1])/2)
weightsEffectVector.bothDiets <- seq(-weightHalfRange,weightHalfRange,0.1)


coxsubsetbetaWeightBoth <- lapply(weightsEffectVector.bothDiets, function(x){
  out <- coxsubsetbetaBoth[coxsubsetbetaBoth$X1 == "weight",]
  out$weightChange <- x
  out$hazardRatio <- exp(out$value)^x
  out
})
coxsubsetbetaWeightBoth <- do.call(rbind, coxsubsetbetaWeightBoth)

coxsubsetbetaWeightBoth$weightDiff <- coxsubsetbetaWeightBoth$weightMax - coxsubsetbetaWeightBoth$weightMin
coxsubsetbetaWeightBoth$weightMargin <- coxsubsetbetaWeightBoth$weightDiff / 2
coxsubsetbetaWeightBoth.plot <- coxsubsetbetaWeightBoth[(coxsubsetbetaWeightBoth$ageAtWeight < 1000)&(coxsubsetbetaWeightBoth$ageAtWeight >50),]
coxsubsetbetaWeightBoth.plot$withinRange <- abs(coxsubsetbetaWeightBoth.plot$weightChange) < coxsubsetbetaWeightBoth.plot$weightMargin
# add text annotations
coxsubsetbetaWeightBoth.plot <- do.call(rbind, lapply(split(coxsubsetbetaWeightBoth.plot, paste0(coxsubsetbetaWeightBoth.plot$modelDiet, coxsubsetbetaWeightBoth.plot$ageAtWeight)), function(x){
  first <- min(which(x$withinRange), na.rm= T)
  last <- max(which(x$withinRange), na.rm= T)
  x$label <- NA
  x$label[c(first, last)] <- round(x$hazardRatio[c(first, last)], 2)
  x
}))
coxsubsetbetaWeightBoth.plot$label50 <- coxsubsetbetaWeightBoth.plot$label
coxsubsetbetaWeightBoth.plot$label50[(coxsubsetbetaWeightBoth.plot$ageAtWeight %% 50) != 0] <- NA



rollingSlopeRange <- range(coxsubsetbetaBoth$rollingSlopeMax -  coxsubsetbetaBoth$rollingSlopeMin)
rollingSlopeHalfRange <- round((rollingSlopeRange[2]-rollingSlopeRange[1])/2)

slopeEffectVector.bothDiets <- seq(-rollingSlopeHalfRange,rollingSlopeHalfRange,0.1)

coxsubsetbetaSlopeBoth <- lapply(slopeEffectVector.bothDiets, function(x){
  out <- coxsubsetbetaBoth[coxsubsetbetaBoth$X1 == "rollingSlope",]
  out$slopeChange <- x
  out$hazardRatio <- exp(out$value)^x
  out
})
coxsubsetbetaSlopeBoth <- do.call(rbind, coxsubsetbetaSlopeBoth)
coxsubsetbetaSlopeBoth$slopeDiff <- coxsubsetbetaSlopeBoth$weightMax - coxsubsetbetaSlopeBoth$weightMin
coxsubsetbetaSlopeBoth$slopeMargin <- coxsubsetbetaSlopeBoth$slopeDiff / 2
coxsubsetbetaSlopeBoth.plot <- coxsubsetbetaSlopeBoth[(coxsubsetbetaSlopeBoth$ageAtWeight < 1000)&(coxsubsetbetaSlopeBoth$ageAtWeight >50),]
coxsubsetbetaSlopeBoth.plot$withinRange <- abs(coxsubsetbetaSlopeBoth.plot$slopeChange) < coxsubsetbetaSlopeBoth.plot$slopeMargin

# add text annotations
coxsubsetbetaSlopeBoth.plot <- do.call(rbind, lapply(split(coxsubsetbetaSlopeBoth.plot, paste0(coxsubsetbetaSlopeBoth.plot$modelDiet, coxsubsetbetaSlopeBoth.plot$ageAtWeight)), function(x){
  first <- min(which(x$withinRange), na.rm= T)
  last <- max(which(x$withinRange), na.rm= T)
  x$label <- NA
  x$label[c(first, last)] <- round(x$hazardRatio[c(first, last)], 2)
  x
}))
coxsubsetbetaSlopeBoth.plot$label50 <- coxsubsetbetaSlopeBoth.plot$label
coxsubsetbetaSlopeBoth.plot$label50[(coxsubsetbetaSlopeBoth.plot$ageAtWeight %% 50) != 0] <- NA



## Now extract the diet parameters
coxsubsetbetaDietBoth <- lapply(1, function(x){
  out <- coxsubsetbetaBoth[coxsubsetbetaBoth$X1 == "DietCodeHFD",]
  out$hazardRatio <- exp(out$value)^x
  out
})
coxsubsetbetaDietBoth <- do.call(rbind, coxsubsetbetaDietBoth)
coxsubsetbetaDietBoth.plot <- coxsubsetbetaDietBoth[(coxsubsetbetaDietBoth$ageAtWeight < 1000)&(coxsubsetbetaDietBoth$ageAtWeight >50),]
coxsubsetbetaDietBoth.plot$Diet <- "HFD"

prange <- sort(-log10(range(c(coxsubsetbetaWeightBoth.plot$p, coxsubsetbetaSlopeBoth.plot$p, coxsubsetbetaDietBoth.plot$p))), decreasing = F)
prange[1] <- floor(prange[1])
prange[2] <- ceiling(prange[2])
hrange <- sort(log2(range(c(coxsubsetbetaWeightBoth.plot$hazardRatio, coxsubsetbetaSlopeBoth.plot$hazardRatio, coxsubsetbetaDietBoth.plot$hazardRatio))), decreasing = F)



weightEffectBoth <- ggplot(coxsubsetbetaWeightBoth.plot, 
                       aes(x = ageAtWeight, fill = log2(hazardRatio), y = weightChange, label = label50)) +
  geom_tile() +
  xlim(0, 1100) +
  scale_fill_gradient2(low="green", mid = "lightgrey", high = "red", midpoint = 0, limits = hrange) +
  scale_color_gradient(low = "lightgrey", high = "black", limits = prange) +
  geom_rug(data = coxsubsetbetaWeightBoth.plot[coxsubsetbetaWeightBoth.plot$weightChange == 0,] , 
           size = 2, sides ="tb",
           aes(color = -log10(p)), y = 1) +
  geom_point(data = coxsubsetbetaWeightBoth.plot[!is.na(coxsubsetbetaWeightBoth.plot$label50),] , shape = 21, color = "black") +
  geom_label_repel(color = "black") +
  geom_hline(yintercept = 0) +
  xlab("Age (days)")+
  ylab("Weight change\n(grams)")+
  guides(color = "none") 
weightEffectBoth

slopeEffectBoth <- ggplot(coxsubsetbetaSlopeBoth.plot, 
                           aes(x = ageAtWeight, fill = log2(hazardRatio), y = slopeChange, label = label50)) +
  geom_tile() +
  # geom_point() +
  # geom_hline(yintercept = exp(0.0), color = "grey") +
  xlim(0, 1100) +
  scale_fill_gradient2(low="green", mid = "lightgrey", high = "red", midpoint = 0, limits = hrange) +
  scale_color_gradient(low = "lightgrey", high = "black", limits = prange) +
  geom_rug(data = coxsubsetbetaSlopeBoth.plot[coxsubsetbetaSlopeBoth.plot$slopeChange == 0,] , 
           size = 2, sides ="tb",
           aes(color = -log10(p)), y = 1) +
  geom_point(data = coxsubsetbetaSlopeBoth.plot[!is.na(coxsubsetbetaSlopeBoth.plot$label50),] , shape = 21, color = "black") +
  geom_label_repel(color = "black") +
  geom_hline(yintercept = 0) +
  xlab("Age (days)")+
  ylab("Weight slope change\n(grams/month)")
slopeEffectBoth


coxsubsetbetaDietBoth.plot$label50 <- round(coxsubsetbetaDietBoth.plot$hazardRatio,2)
coxsubsetbetaDietBoth.plot$label50[(coxsubsetbetaDietBoth.plot$ageAtWeight %% 50) != 0] <- NA
dietEffectBoth <- ggplot(coxsubsetbetaDietBoth.plot[coxsubsetbetaDietBoth.plot$Diet == "HFD",], 
                           aes(x = ageAtWeight, fill = log2(hazardRatio), y = Diet, label = label50)) +
  geom_tile() +
  # geom_point() +
  # geom_hline(yintercept = exp(0.0), color = "grey") +
  xlim(0, 1100) +
  scale_fill_gradient2(low="green", mid = "lightgrey", high = "red", midpoint = 0, limits = hrange) +
  scale_color_gradient(low = "lightgrey", high = "black", limits = prange) +
  geom_rug(data = coxsubsetbetaDietBoth.plot[coxsubsetbetaDietBoth.plot$Diet == "HFD",] , 
           size = 2, sides ="tb",
           aes(color = -log10(p)), y = 1) +
  geom_point(data = coxsubsetbetaDietBoth.plot[!is.na(coxsubsetbetaDietBoth.plot$label50),] , shape = 21, color = "black") +
  geom_label_repel(data = coxsubsetbetaDietBoth.plot[coxsubsetbetaDietBoth.plot$Diet == "HFD",] , color = "black") +
  xlab("Age (days)")+
  ylab("Diet")+
  guides(color = "none") 
dietEffectBoth


plot_grid(dietEffectBoth, weightEffectBoth, slopeEffectBoth, nrow = 3, rel_heights = c(1,2,2), align = "hv", labels = LETTERS[1:3]) +
  ggsave("./Plots/diet.vs.bw.hazardEffect.pdf",width = 10, height= 10)
# same as before, but add the p-values separately
pvaluesplotBoth.data <- lmerWeightReshaped[(lmerWeightReshaped$X2 == "p")&
                                             (!grepl("frail", lmerWeightReshaped$X1))&
                                             (grepl("cmodel", lmerWeightReshaped$model))&
                                             (lmerWeightReshaped$modelWithSlope == 'with slope')&
                                             (lmerWeightReshaped$modelDiet == "Both")&
                                             (lmerWeightReshaped$X1 != "daysOnDiet") &
                                           (lmerWeightReshaped$ageAtWeight >50)&
                                             (lmerWeightReshaped$ageAtWeight <1000),]
pvaluesplotBoth.data$X1 <- plyr::mapvalues(pvaluesplotBoth.data$X1, from = c("DietCodeHFD", "weight", "rollingSlope"), to = c("Diet (HFD)", "Weight", "rollingSlope"))
pvaluesplotBoth <- ggplot(pvaluesplotBoth.data, 
                  aes(x = ageAtWeight, y = -log10(value), group = X1, color = X1)) +
  geom_line() +
  geom_hline(yintercept = -log10(0.05), color = 'grey') +
  geom_hline(yintercept = -log10(0.05 / length(unique(pvaluesplotBoth.data$ageAtWeight)))) +
  xlim(0, 1100) +
  xlab("Age (days)") +
  scale_color_discrete(name = "Model parameter") +
  ylab(expression(-log[10](p-value)))

plot_grid(pvaluesplotBoth,dietEffectBoth, weightEffectBoth, slopeEffectBoth, nrow = 4, rel_heights = c(2,1,2,2), align = "hv", labels = LETTERS[1:4]) +
  ggsave("./Plots/diet.vs.bw.hazardEffect.withPval.pdf",width = 10, height= 10)


# make an excel sheet with the model outputs
toOutput <- do.call(rbind, mapply(x=lmerWeight,y=names(lmerWeight), function(x,y){
  do.call(rbind, mapply(z=x, w = names(x), function(z,w){
    
    out <- reshape::melt(summary(z)$coef)
    
    out$model <- w
    out$ageAtWeight <- as.numeric(y)
    
    if(w %in% c("cmodelslope", "cmodelCDslope", "cmodelHFDslope")){
      
      out$weightMin = z$weightMin
      out$weightMax = z$weightMax
      out$rollingSlopeMin = z$rollingSlopeMin
      out$rollingSlopeMax = z$rollingSlopeMax
      out$nmice <- nrow(z$data)
      out
    }else{
      
      NULL
    }
    
    
  }, SIMPLIFY = FALSE))
}, SIMPLIFY = FALSE))

toOutput_reshaped <- reshape::cast(toOutput, ageAtWeight + model~ X1 + X2)
# remove: frailty(StrainNameCurrent_se(coef) frailty(StrainNameCurrent_se2
toOutput_reshaped$`frailty(StrainNameCurrent_se(coef)` <- NULL
toOutput_reshaped$`frailty(StrainNameCurrent_se2` <- NULL
toOutput_reshaped$model <- plyr::mapvalues(toOutput_reshaped$model,
                                           c("cmodelslope", "cmodelCDslope", "cmodelHFDslope"),
                                           c("Cox model (Both Diets)",
                                             "Cox model (CD only)",
                                             "Cox model (HFD only)"))

write.csv(toOutput_reshaped, "./Data/coxModelResults.csv", row.names = F)
