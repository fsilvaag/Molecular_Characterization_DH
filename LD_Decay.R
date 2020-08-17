
#How to plot LD decay!!!
#Adpated by Fernando SIlva Aguilar from code 
# giveng by Arthur Pereira da Silva
#Iowa State University
#Department of Agronomy

rm(list = ls())
library(ggplot2)
library(MASS)
library(scales)
library(data.table)
library(tidyverse)
library(gtools)

#Import the LD table - output of Tassel (File with the columns Dist_bp and	R^2 only)
path1 = "C:/Users/deer1/Google Drive/Alejandro/01. Filtrado/"
All <- as.data.frame(fread(paste0(path1,"LD_Decay.txt"), 
                           header = T))

# Remove empty cells in the R squared value
All1 <- All %>% drop_na(`R^2`)
All1_1 <- All1[,c("Locus1","Locus2","Dist_bp","R^2")]
colnames(All1_1) <- c("Locus1","Locus2","dist", "rsq")

#Specify the number of individuals in your analysis
# only change "n"

n = 502 
Cstart <- c(C=0.1)

# Subset data per chromosme
for (i in 1:10) {
  assign(paste0("Ch",i),
         subset(All1_1, All1_1$Locus1 == i & All1_1$Locus2 ==i))
}
#convert dist column to numeric variable
l <- mixedsort(ls()[grep("Ch",ls())])
for (i in 1:10) {
  d <- get(l[i])
  d$dist <- as.numeric(paste0(d$dist))
  assign(paste0(l[i]),d)
}


#Fit the model
for (i in 1:10) {
  d <- get(l[i])
  assign(paste0("modelCh",i),
         nls(rsq ~ ((10+C*dist)/((2+C*dist)*(11+C*dist)))*(1+((3+C*dist)*(12+12*C*dist+(C*dist)^2))/(n*(2+C*dist)*(11+C*dist))), 
             data=d, 
             start=Cstart, 
             control=nls.control(maxiter=100)))
}

#Extract recombination parameter
m <- mixedsort(ls()[grep("modelCh",ls())])
for (i in 1:10) {
  d <- get(m[i])
  assign(paste0("rho",i),
         summary(d)$parameters[1])
}


#Adjusted LD values
p <- mixedsort(ls()[grep("rho",ls())])
for (i in 1:10) {
  rho <- get(p[i])
  d <- get(l[i])
  assign(paste0("newrsq",i),
         ((10+rho*d$dist)/((2+rho*d$dist)*(11+rho*d$dist)))*
           (1+((3+rho*d$dist)*(12+12*rho*d$dist+(rho*d$dist)^2))/
              (n*(2+rho*d$dist)*(11+rho*d$dist))))
}

y <- mixedsort(ls()[grep("newrsq",ls())])
for (i in 1:10) {
  ne <- get(y[i])
  d <- get(l[i])
  assign(paste0("newfileCh",i),
         data.frame(d$dist, ne))
}

e <- mixedsort(ls()[grep("newfile",ls())])
for (i in 1:10) {
  d <- get(e[i])
  assign(paste0("maxld",i),
         max(d$ne))
}

halfdecay = 0.2

for (i in 1:10) {
  d <- get(e[i])
  assign(paste0("halfdecaydistCh",i),
         d$d.dist[which.min(abs(d$ne-halfdecay))])
}

#halfdecaydist <- newfile$Ch1.dist[which.min(abs(newfile$newrsq-halfdecay))]
for (i in 1:10) {
  d <- get(e[i])
  assign(paste0(e[i]),
         d[order(d$d.dist),])
}

########################  TABLE FOR LD DECAY

li <- mixedsort(ls()[grep("halfdecaydistCh",ls())])
for (i in 1:10) {
  d <- as.data.frame(get(li[i]))
  assign(paste0("halfdecaydistChr_",i),
         d)
}

lk <- mixedsort(ls()[grep("halfdecaydistChr",ls())])
list_df = lapply(lk, get)
s <- rbindlist(list_df)


s$chrom <- 1:10
s <- s[,c("chrom","get(li[i])")]
colnames(s)[2] <- "LD_decay_bp"
s$LD_decay_Kb <- s$LD_decay_bp/1000

fwrite(s,paste0(path1,"ld_decay_table.csv"),
       quote = F, row.names = T, sep = ",")



# Graph example
#g1 <- ggplot(data = newfile, aes(x = newfile$Ch10.dist, y = newfile$newrsq)) +
#  geom_line() +
#  geom_hline(aes(yintercept=0.2), color = "blue") +
#  geom_vline(aes(xintercept=halfdecaydist), color = "blue") +
#  scale_x_continuous(limits = c(1,1500000), breaks = c(1,100000,300000,600000,900000,1200000),labels = c("1","100","300","600","900","1200")) +
#  labs(x = "Distance Kb", y = expression(r^2)) +
#  theme_classic()
#ggsave("LD_Ch10.png", plot = g1, width = 10, height = 10, units = "cm")
#write.table(newfile,"LD10.txt")



