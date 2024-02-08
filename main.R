# package and data loading ------------------------------------------------

#library loading
#WARNING: it is needed to have to have the correct versions of the library, in case of some old versions method pmm may cause code abort
library(multcomp)
library(mice)
library(withr)
library(Rcpp)
library(tidyverse)
library(corrplot)
library(GGally)
library(outliers)
library(EnvStats)

#loading data
z17 <- read.table("https://web.sgh.waw.pl/~akorczy/files/abappm/data/t1/z17.csv", sep=",", header = TRUE)

#checking basic patterns, statistics -----------------------------

#checking column types
str(z17)

#a pattern of missing data
md.pattern(z17)

#basic statistical values for continuous variables
summary(z17[c('y_t1', 'x1', 'y_t2')])
sd(z17$y_t2, na.rm = TRUE)
sapply(z17[c('y_t1', 'x1')], sd)

#skew and kurtosis
moments::skewness(z17[c('y_t1', 'x1', 'y_t2')], na.rm = TRUE)
moments::kurtosis(z17[c('y_t1', 'x1', 'y_t2')], na.rm = TRUE)

#values for discreet variables
table(z17$x2)
table(z17$trt)

#data correlation with heatmap visualization
cor_mtrx <- cor(z17, use = "complete.obs")
cor_mtrx
corrplot(cor_mtrx, method = "color")

#data covariance
cov(z17, use = "complete.obs")

#transforming data for visualization
z17 <- z17 %>%
  mutate(m = ifelse(is.na(y_t2), 1, 0))

z17 <- z17 %>%
  mutate(m0 = ifelse(is.na(y_t2), 0, y_t2))

z17 <- z17 %>%
  mutate(chg = y_t2 - y_t1,
         chg2 = if_else(is.na(chg), -8, chg))


#visualizing data and looking for data patterns ------------------

#checking counts of missing values in values of trt
ggplot(z17, aes(x = trt, fill = factor(m))) +
  geom_bar(position = "dodge") +
  labs(title = "Count Plot", x = "trt", y = "Count", fill = "m")

#checking counts of missing values in values of x2
ggplot(z17, aes(x = trt, fill = factor(m))) +
  geom_bar(position = "dodge") +
  labs(title = "Count Plot", x = "x2", y = "Count", fill = "m")

#histogram of x1 colored based on missing values in y_t2
ggplot(z17, aes(x = x1, fill = factor(m))) +
  geom_histogram(binwidth = 1, position = "dodge", color = "white") +
  labs(title = "Histogram Plot", x = "x1", y = "Frequency", fill = "m")

#histogram of y_t1 colored based on missing values in y_t2
ggplot(z17, aes(x = x1, fill = factor(m))) +
  geom_histogram(binwidth = 1, position = "dodge", color = "white") +
  labs(title = "Histogram Plot", x = "y_t1", y = "Frequency", fill = "m")

#comparing shapes of y_t1 and y_t2
z17 %>%
  gather(key = "variable", value = "value", y_t1, y_t2) %>%
  ggplot(aes(x = value, fill = variable)) +
  geom_histogram(alpha = 0.4, position = "identity", bins = 30) +
  labs(title = "Histogram Plot", x = "Values", y = "Frequency") +
  theme_minimal() +
  scale_fill_manual(values = c(rgb(0.2, 0.4, 0.6, alpha = 0.4), rgb(0.6, 0.2, 0.4, alpha = 0.4)))

#a visual representation of relationship between change and y_t1
plot(z17$y_t1, z17$chg2, pch=20, col = factor(z17$trt), xlab = "y_t1", ylab = "Change")
plot(z17$y_t1, z17$chg2, pch=20, col = factor(z17$m), xlab = "y_t1", ylab = "Change")

#pairplot - checking pairs of every two variables, colored based on missing values in y_t2
ggpairs(z17, columns = c("y_t1", "x1", "trt", "x2", "y_t2"),
        aes(color = as.factor(m), alpha=0.4))


#comparing differences ------------------

#continuous
t.test(subset(z17, m==0)$y_t1, subset(z17, m==1)$y_t1)
t.test(subset(z17, m==0)$x1, subset(z17,m==1)$x1)

#discrete
chisq.test(z17$x2, z17$m, correct = FALSE)
chisq.test(z17$trt, z17$m, correct = FALSE)

#testing for outliers ------------------
#Grubbs Test
grubbs.test(z17$y_t1) 
grubbs.test(z17$y_t1, opposite = TRUE) 
grubbs.test(z17$y_t2)
grubbs.test(z17$y_t2, opposite = TRUE)
grubbs.test(z17$x1)
grubbs.test(z17$x1, opposite = TRUE)

#Rosner Test
rosnerTest(z17$y_t1)
rosnerTest(z17$y_t2)
rosnerTest(z17$x1)


#histogram of y_t1
ggplot(z17) +
  geom_histogram(aes(x=y_t1), fill='cornflowerblue', col='white') +
  labs(title = "Histogram of y_t1") +
  theme_minimal()

#histogram of y_t2
ggplot(z17) +
  geom_histogram(aes(x=y_t2), fill='cornflowerblue', col='white') +
  labs(title = "Histogram of y_t2") +
  theme_minimal()

#histogram of x1
ggplot(z17) +
  geom_histogram(aes(x=x1), fill='cornflowerblue', col='white') +
  labs(title = "Histogram of x1") +
  theme_minimal()


#checking data model - MCAR, MAR -----------------------------------------

#MCAR
#indicating reference category
z17$trt <- relevel(factor(z17$trt), ref = 2)
options(contrasts=c("contr.treatment","contr.poly")) 

#ANCOVA analysis
ancova01 <- lm(chg ~ y_t1 + trt, data = z17)
anova(ancova01)
summary(ancova01)

#test post-hoc Tukey HSD
postHocs <- glht(ancova01, linfct = mcp(trt = "Tukey")) 
summary(postHocs)


#MAR
#Analysis of different imputation options

#setting seed
seed1 <- 1234
#assigning clean dataset
z17new <- z17[, c("y_t1", "trt", "x1", "x2", "y_t2")]
missing_vector <- z17[,c("y_t1", "m")]

#LINEAR REGRESSION
#multivariate imputations based on linear regression with 25 multiple imputations
imp01 <- mice(z17new, method = "norm.predict", m = 25, seed=seed1) 
#filling the dataset
matmi<-complete(imp01,"long")
matmi <- merge(matmi, missing_vector, by="y_t1")
#calculating change after imputation
matmi$chg <- matmi$y_t2-matmi$y_t1 
#analyzing the imputation technique on scatter plot
xyplot(chg~y_t1,matmi,group=m
       ,par.settings = list(superpose.symbol = list(col = c("blue","red"),pch = c(1,3),cex=1.2))
       ,main="Regression imputation"
       ,key=list(title=" "
                 ,space="right"
                 ,points=list(pch=c(1,3),col=c("blue","red")),text=list(c('Y-known','Y-imputed')),cex.title=0.5,cex=0.9)) 

xyplot(y_t2~y_t1,matmi,group=m
       ,par.settings = list(superpose.symbol = list(col = c("blue","red"),pch = c(1,3),cex=1.2))
       ,main="Regression imputation"
       ,key=list(title=" "
                 ,space="right"
                 ,points=list(pch=c(1,3),col=c("blue","red")),text=list(c('Y-known','Y-imputed')),cex.title=0.5,cex=0.9)) 


#BAYESIAN LINEAR REGRESSION
imp02 <- mice(z17new, method = "norm", m = 25, seed=seed1) 
matmi2<-complete(imp02,"long")
matmi2 <- merge(matmi2, missing_vector, by="y_t1")
matmi2$chg <- matmi2$y_t2-matmi2$y_t1

xyplot(chg~y_t1,matmi2,group=m
       ,par.settings = list(superpose.symbol = list(col = c("blue","red"),pch = c(1,3),cex=1.2))
       ,main="Bayesian linear regression imputation"
       ,key=list(title=" "
                 ,space="right"
                 ,points=list(pch=c(1,3),col=c("blue","red")),text=list(c('Y-known','Y-imputed')),cex.title=0.5,cex=0.9))

xyplot(y_t2~y_t1,matmi2,group=m
       ,par.settings = list(superpose.symbol = list(col = c("blue","red"),pch = c(1,3),cex=1.2))
       ,main="Bayesian linear regression imputation"
       ,key=list(title=" "
                 ,space="right"
                 ,points=list(pch=c(1,3),col=c("blue","red")),text=list(c('Y-known','Y-imputed')),cex.title=0.5,cex=0.9)) 


#DECISION TREE
imp03 <- mice(z17new, method = "cart", m = 25, seed=seed1) 
matmi3<-complete(imp03,"long") 
matmi3 <- merge(matmi3, missing_vector, by="y_t1")
matmi3$chg <- matmi3$y_t2-matmi3$y_t1 

xyplot(chg~y_t1,matmi3,group=m
       ,par.settings = list(superpose.symbol = list(col = c("blue","red"),pch = c(1,3),cex=1.2))
       ,main="Decision tree imputation"
       ,key=list(title=" "
                 ,space="right"
                 ,points=list(pch=c(1,3),col=c("blue","red")),text=list(c('Y-known','Y-imputed')),cex.title=0.5,cex=0.9))

xyplot(y_t2~y_t1,matmi3,group=m
       ,par.settings = list(superpose.symbol = list(col = c("blue","red"),pch = c(1,3),cex=1.2))
       ,main="Decision tree imputation"
       ,key=list(title=" "
                 ,space="right"
                 ,points=list(pch=c(1,3),col=c("blue","red")),text=list(c('Y-known','Y-imputed')),cex.title=0.5,cex=0.9))


#HOT DECK
imp04 <- mice(z17new, method = "pmm", m = 25, seed=seed1) 
matmi4<-complete(imp04,"long")
matmi4 <- merge(matmi4, missing_vector, by="y_t1")
matmi4$chg <- matmi4$y_t2-matmi4$y_t1

xyplot(chg~y_t1,matmi4,group=m
       ,par.settings = list(superpose.symbol = list(col = c("blue","red"),pch = c(1,3),cex=1.2))
       ,main="Hot deck imputation"
       ,key=list(title=" "
                 ,space="right"
                 ,points=list(pch=c(1,3),col=c("blue","red")),text=list(c('Y-known','Y-imputed')),cex.title=0.5,cex=0.9)) 

xyplot(y_t2~y_t1,matmi4,group=m
       ,par.settings = list(superpose.symbol = list(col = c("blue","red"),pch = c(1,3),cex=1.2))
       ,main="Hot deck imputation"
       ,key=list(title=" "
                 ,space="right"
                 ,points=list(pch=c(1,3),col=c("blue","red")),text=list(c('Y-known','Y-imputed')),cex.title=0.5,cex=0.9))


# covariance model ------------------------------------------------------------

#MAR covariance model estimated for chosen imputation method
fit <- with(imp04, lm(y_t2-y_t1 ~ y_t1 + trt))
tab <- summary(pool(fit))
tab


# MNAR --------------------------------------------------------------------

#creating a function for automatically including several elements
sensitivityanalysis <- function(vector) {
  for (i in vector) {
    
    imp11 <- mice(z17new,
                  method = "mnar.norm",
                  blots = list(y_t2 = list(ums = i)),
                  seed = seed1)
    matmi11<-complete(imp11,"long")
    matmi11 <- merge(matmi11, missing_vector, by="y_t1")
    matmi11$chg <- matmi11$y_t2-matmi11$y_t1
    
    print(xyplot(chg~y_t1,matmi11,group=m
                 ,par.settings = list(superpose.symbol = list(col = c("blue","red"),pch = c(1,3),cex=1.2))
                 ,main=paste0("MNAR ums=",  ums=i)
                 ,key=list(title=" "
                           ,space="right"
                           ,points=list(pch=c(1,3),col=c("blue","red")),text=list(c('Y-known','Y-imputed')),cex.title=0.5,cex=0.9))) 
    
    print(xyplot(y_t2~y_t1,matmi11,group=m
                 ,par.settings = list(superpose.symbol = list(col = c("blue","red"),pch = c(1,3),cex=1.2))
                 ,main=paste0("MNAR ums=",  ums=i)
                 ,key=list(title=" "
                           ,space="right"
                           ,points=list(pch=c(1,3),col=c("blue","red")),text=list(c('Y-known','Y-imputed')),cex.title=0.5,cex=0.9)))
    fit <- with(imp11, lm(y_t2-y_t1 ~ y_t1 + trt))
    tab <- summary(pool(fit))
    print(tab)
  }
}

#intercept shift sensitivity analysis
sensitivityanalysis(c("-5", "-4", "-3", "-2", "-1", "0", "1", "2", "3", "4", "5"))
