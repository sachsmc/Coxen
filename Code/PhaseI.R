
library(data.table)
library(samr)
library(lars)
library(reshape2)
library(ggplot2)

celline <- data.table(read.csv("../Data/nci60-cisplatin-lowdoserange.csv", stringsAsFactors = FALSE), 
                      key = c("CellPanelNbr", "CellLineNbr"))

expression <- data.table(read.csv("../Data/WEB_DATA_GENELOGIC_U133/WEB_DATA_GENELOGIC_U133.TXT", 
                                  stringsAsFactors = FALSE, 
                                  colClasses = c(rep("character", 4), 
                                                 "integer", "integer", "integer", 
                                        rep("character", 5), rep("numeric", 2)), na.string = "."),
                         key = c("PANELNBR", "CELLNBR"))[CHIP=="U133A"]

both <- celline[expression]

rm(expression)
### reproducing the coxen algorithm

leaveout <- both[CellPanelName %in% c("Prostate", "Renal")]
training <- both[!CellPanelName %in% c("Prostate", "Renal")]

## step one, select 12 most sensitive and 12 most resistant cell lines

lowhighcuts <- quantile(celline[!CellPanelName %in% c("Prostate", "Renal")]$logValue, c(.2, .8))
training[, high := ifelse(logValue <= lowhighcuts[1], 1, ifelse(logValue >= lowhighcuts[2], 2, NA))]

extremes <- training[!is.na(high)]
featmat <- dcast.data.table(extremes, cellname ~ FEATURE_ID, value.var = "VALUE", fun = function(x) x[1])

highlook <- unique(extremes[, list(cellname, high)])
hilook <- highlook$high
names(hilook) <- highlook$cellname


### use SAM to find subset of genes

X <- t(as.matrix(as.data.frame(featmat)[, -1]))
y <- hilook[featmat[, cellname]]

samfit <- SAM(X, y, resp.type = "Two class unpaired", fdr.output = .1)

cand_sens <- with(samfit$siggenes.table, as.numeric(c(genes.up[,2], genes.lo[,2])))
cand_sens2 <- rownames(X)[cand_sens]
plot(samfit)

### coexpression matrix

U1 <- dcast.data.table(training[FEATURE_ID %in% cand_sens2], cellname ~ FEATURE_ID, value.var = "VALUE", fun = function(x) x[1])

V1 <- dcast.data.table(leaveout[FEATURE_ID %in% cand_sens2], cellname ~ FEATURE_ID, value.var = "VALUE", fun = function(x) x[1])

rj <- diag(cor(cor(as.data.frame(U1)[,-1]),
         cor(as.data.frame(V1)[,-1])))

keepers <- names(rj)[rj > .5]


####
## forming predictions

trainpred <- dcast(training[FEATURE_ID %in% keepers], cellname + logValue ~ FEATURE_ID, value.var = "VALUE")
validate <- dcast(leaveout[FEATURE_ID %in% keepers], cellname + logValue ~ FEATURE_ID, value.var = "VALUE")

trainfit <- lm(as.formula(paste("logValue ~" , paste(paste0("`", keepers, "`"), collapse = "+"))), data = trainpred)

valid <- predict(trainfit, newdata = validate)
qplot(valid, validate$logValue) + stat_smooth(method = "lm") + geom_abline(intercept = 0, slope = 1)

mse <- sqrt(mean((valid - validate$logValue)^2))

### compare to lasso

lassfeatures <- dcast.data.table(training, cellname + logValue ~ FEATURE_ID, value.var = "VALUE", fun = mean)
validlass <- dcast.data.table(leaveout, cellname + logValue ~ FEATURE_ID, value.var = "VALUE", fun = mean)
validX <- as.data.frame(validlass)[, -c(1,2)]

ylass <- lassfeatures$logValue
Xlass <- as.data.frame(lassfeatures)[, -c(1,2)]

lassfit <- lars(as.matrix(Xlass), ylass, trace = TRUE, use.Gram = FALSE)

predlas <- predict(lassfit, newx = validX, type = "fit", s = 42)$fit
mselass <- sqrt(mean((predlas - validlass$logValue)^2))
mselass

### supervised principle components
## screen features on association with response

ks <- seq(5, 35, by = 5)
npcs <- 2:15

setkey(training, FEATURE_ID)
train.assoc <- training[, list(assoc = lm(logValue ~ VALUE)$coeff[2]), by = FEATURE_ID]

cv.dat <- function(k, npc, train.assoc){

train.super <- train.assoc[abs(assoc) > quantile(abs(assoc), 1-k/22261), FEATURE_ID]
geneset.train <- training[FEATURE_ID %in% train.super, list(FEATURE_ID, cellname, logValue, VALUE)]
geneset.valid <- leaveout[FEATURE_ID %in% train.super, list(FEATURE_ID, cellname, logValue, VALUE)]

prinX <- dcast.data.table(rbind(geneset.train, geneset.valid), cellname + logValue ~ FEATURE_ID, value.var = "VALUE", fun = mean)
components <- princomp(as.data.frame(prinX)[,-c(1:2)], scores = FALSE)

###  model based on the PCs

X.train <- dcast.data.table(geneset.train, 
                 cellname + logValue ~ FEATURE_ID, 
                 value.var = "VALUE", fun = mean)
train.pcs <- predict(components, 
                     newdata = X.train)[,1:npc]
train.fit <- lm(as.formula(paste("logValue ~ ", paste(paste("Comp", 1:npc, sep = "."), collapse = "+"))), data = cbind(X.train, train.pcs))


X.valid <- dcast.data.table(geneset.valid, cellname + logValue ~ FEATURE_ID, 
                            value.var = "VALUE", fun = mean)

valid.pcs <- predict(components, newdata = X.valid)[,1:npc]
valid.est <- predict(train.fit, newdata = cbind(X.valid, valid.pcs))
sqrt(mean(X.valid$logValue - valid.est)^2)
}

dout <- NULL
for(k in ks){ for(n in npcs){
  
  if(n > k) next
  dout <- rbind(dout, data.frame(k = k, npcs = n, rmse = cv.dat(k, n, train.assoc)))
  
}}
