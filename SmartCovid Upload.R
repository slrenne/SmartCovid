#libraries
library(rethinking)
library(tidyverse)
library(RColorBrewer)
library(gridExtra)

#dbs
write.table(performance, file = "/Users/lollo/Dropbox/SmartCovid/performance.txt",sep="\t")
write.table(psyco, file = "/Users/lollo/Dropbox/SmartCovid/psyco.txt",sep="\t")
performance <- read.delim("~/Dropbox/SmartCovid/dbs/performance.txt")
psyco <- read.delim("~/Dropbox/SmartCovid/dbs/psyco.txt")

#data wrangling
pf2 <- performance %>% gather(key = "cases", value = "wrong", "BR1.1":"GI2.5")
pf2$PID <- as.factor(pf2$PID)
pf2$PID <- factor(pf2$PID, c("P1","P2", "P3", "P4", "P5", "P6", "P7", "P8", "P9",
 "P10","P11","P12","P13","P14","P15","P16","P17"))
pf2$cases <- as.factor(pf2$cases)
pf2$spec <- ifelse(pf2$cases %in% c("BR1.1", "BR1.2", "BR1.3", "BR1.4", 
                                    "BR1.5", "BR2.2", "BR2.3", "BR2.4", "BR2.5"), "Breast", 
                   ifelse(pf2$cases %in% c("GI1.1", "GI1.2", "GI1.3", "GI1.4", "GI1.5", "GI2.1", "GI2.2",
                                    "GI2.3", "GI2.4", "GI2.5"), "GI", 
                          ifelse(pf2$cases %in% c("URO1",  "URO2",  "URO3",  "URO4",  "URO5"), "Uro", "ERR")))
sum(pf2$spec=="ERR")  
pf2$spec <- as.factor(pf2$spec)  
pf2$category <- as.factor(pf2$category)
pf2$category <- factor(pf2$category, c("nplVSnot","benVSmal","diagnosis","histotype", "grade"))
pf2$biopsy <- ifelse(pf2$cases %in% c("BR1.2", "BR1.4", "BR1.5", "BR2.2", 
        "BR2.3", "URO4", "URO5", "GI1.1", "GI2.2", "GI2.4", "GI2.5"), 0L, 1L)

#####################################################################################################
#graphical exploration of data and summary

#Figures for R1
#figure 1
png("/Users/lollo/Dropbox/SmartCovid/paper/3 JIMR/R1/Figure1.png", units = "in", width = 6, height = 5, res = 300)
plot(NULL, xlim = c(0.5,5.5),ylim = c(0,0.15), xaxt ="n", yaxt ="n", 
     ylab = "", xlab = "", main = "Errors among Categories, Median (IQR)")
points(Tables$Median, pch=16)  
for (i in 1:5) lines(c(i,i), c(Tables$f_qu[i], Tables$t_qu[i]),
                     lwd = 2)
axis( 1 , at=1:5 , labels = labs)
axis( 2 , seq(from=0,to=0.15,by=0.03) , labels = paste(0:5*3,"%",sep=""))
dev.off()

#Figure 2
plot2b <- pf2 %>% filter(!is.na(wrong)) %>%
  filter(category=="nplVSnot") %>%
  ggplot(aes(as.factor(Plevels), fill= as.factor(wrong))) +
  geom_bar(position = "fill", stat = "count")  +
  labs (x="",
        y="Neoplastic Vs Not")+
  scale_x_discrete(breaks = c("1","2","3","4"), labels = c("Resident", "Junior","Expert", "Senior"))+
  scale_y_continuous(breaks=seq(from=0, to=1, by=0.25), 
                     labels = c("0%","25%","50%","75%","100%"))+
  scale_fill_manual(name="", breaks = c("1","0"), 
                    labels = c("Wrong", "Correct"), values=c("#de2d26","#2ca25f"))+
  theme(legend.position = "none")

plot2c <- pf2 %>% filter(!is.na(wrong)) %>%
  filter(category=="nplVSnot") %>%
  ggplot(aes(as.factor(biopsy), fill= as.factor(wrong))) +
  geom_bar(position = "fill", stat = "count") +
  labs (x="", y="")+
  scale_x_discrete(breaks = c("0","1"), labels = c("Surgery","Biopsy"))+
  scale_y_continuous(breaks=seq(from=0, to=1, by=0.25), 
                     labels = c("0%","25%","50%","75%","100%"))+
  scale_fill_manual(name="", breaks = c("1","0"), 
                    labels = c("Wrong", "Correct"), values=c("#de2d26", "#2ca25f"))

plot3b <-  pf2 %>% filter(!is.na(wrong)) %>%
  filter(category=="benVSmal") %>%
  ggplot(aes(as.factor(Plevels), fill= as.factor(wrong))) +
  geom_bar(position = "fill", stat = "count")  +
  labs (x="Pathologist's Category",
        y="Benign Vs Malignant")+
  scale_x_discrete(breaks = c("1","2","3","4"), labels = c("Resident", "Junior","Expert", "Senior"))+
  scale_y_continuous(breaks=seq(from=0, to=1, by=0.25), 
                     labels = c("0%","25%","50%","75%","100%"))+
  scale_fill_manual(name="", breaks = c("1","0"), 
                    labels = c("Wrong", "Correct"), values=c("#de2d26","#2ca25f"))+
  theme(legend.position = "none")

plot3c <- pf2 %>% filter(!is.na(wrong)) %>%
  filter(category=="benVSmal") %>%
  ggplot(aes(as.factor(biopsy), fill= as.factor(wrong))) +
  geom_bar(position = "fill", stat = "count") +
  labs (x="Specimen's Kind",
        y="")+
  scale_x_discrete(breaks = c("0","1"), labels = c("Surgery","Biopsy"))+
  scale_y_continuous(breaks=seq(from=0, to=1, by=0.25), 
                     labels = c("0%","25%","50%","75%","100%"))+
  scale_fill_manual(name="", breaks = c("1","0"), 
                    labels = c("Wrong", "Correct"), values=c("#de2d26","#2ca25f"))
png("/Users/lollo/Dropbox/SmartCovid/paper/3 JIMR/R1/Fig2.png", units = "in", width = 6, height = 4, res = 300)
grid.arrange(plot2b, plot2c, plot3b, plot3c, nrow=2, ncol=2, top = "Proportion of Errors")
dev.off()

#Table 1
pf2 %>% filter(!is.na(wrong)) %>%
  group_by(PID) %>% 
  summarise(prop = mean(wrong), wrong = sum(wrong), tot = n())  
pf2 %>% filter(!is.na(wrong)) %>%
  group_by(Plevels) %>% 
  summarise(prop = mean(wrong), wrong = sum(wrong), tot = n()) 
pf2 %>% filter(!is.na(wrong)) %>%
  group_by(category) %>% 
  summarise(prop = mean(wrong), wrong = sum(wrong), tot = n())
pf2 %>% filter(!is.na(wrong)) %>%
  group_by(biopsy) %>% 
  summarise(prop = mean(wrong), wrong = sum(wrong), tot = n()) 
pf2 %>% filter(!is.na(wrong)) %>%
  group_by(spec) %>% 
  summarise(prop = mean(wrong), wrong = sum(wrong), tot = n()) 

##Supplemntary 
# plotting the results for each pathologist

p1 <- pf2 %>% filter(!is.na(wrong)) %>%
  ggplot(aes(PID, fill= as.factor(wrong))) +
  geom_bar(position = "fill", stat = "count") +
  labs (x="",
        y="Proportion of cases")+
  scale_x_discrete(breaks = c("P1","P2", "P3", "P4", "P5", "P6", "P7", "P8", "P9",
                              "P10","P11","P12","P13","P14","P15","P16","P17"), 
                   labels = 1:17)+
  scale_y_continuous(breaks=seq(from=0, to=1, by=0.25), 
                     labels = c("0%","25%","50%","75%","100%"))+
  scale_fill_manual(name="", breaks = c("1","0"), 
                    labels = c("Wrong", "Correct"), values=c("#de2d26","#2ca25f"))+ 
  theme(legend.position = "none")

p2<-pf2 %>% filter(!is.na(wrong)) %>%
  ggplot(aes(PID, fill= as.factor(wrong))) +
  geom_bar(position = "fill", stat = "count") +
  labs (x="",
        y="")+
  scale_x_discrete(breaks=NULL)+
  scale_y_continuous(breaks=seq(from=0, to=1, by=0.25), 
                    labels = c("0%","25%","50%","75%","100%"))+
  scale_fill_manual(name="", breaks = c("1","0"), 
                    labels = c("Wrong", "Correct"), values=c("#de2d26","#2ca25f"))+
  facet_wrap(vars(category))+ 
  theme(legend.position = "none")

Slevels.labs <- c("Surgery", "Biopsy")
names(Slevels.labs) <- c(0,1)
p3<-pf2 %>% filter(!is.na(wrong)) %>%
  ggplot(aes(PID, fill= as.factor(wrong))) +
  geom_bar(position = "fill", stat = "count") +
  labs (x="Pathologists",
        y="Proportion of cases")+
  scale_x_discrete(breaks=NULL)+
  scale_y_continuous(breaks=seq(from=0, to=1, by=0.25), 
                    labels = c("0%","25%","50%","75%","100%"))+
  scale_fill_manual(name="", breaks = c("1","0"), 
                    labels = c("Wrong", "Correct"), values=c("#de2d26","#2ca25f"))+
  facet_wrap(vars(biopsy), labeller = labeller(biopsy = Slevels.labs))+ 
  theme(legend.position = "none")

p4<-pf2 %>% filter(!is.na(wrong)) %>%
  ggplot(aes(PID, fill= as.factor(wrong))) +
  geom_bar(position = "fill", stat = "count") +
  labs (x="Pathologists",
        y="")+
  scale_x_discrete(breaks=NULL)+
  scale_y_continuous(breaks=seq(from=0, to=1, by=0.25), 
                    labels = c("0%","25%","50%","75%","100%"))+
  scale_fill_manual(name="", breaks = c("1","0"), 
                    labels = c("Wrong", "Correct"), values=c("#de2d26","#2ca25f"))+
  facet_wrap(vars(spec))+ 
  theme(legend.position = "none")

png("/Users/lollo/Dropbox/SmartCovid/paper/3 JIMR/R1/Supp_PID.png", 
    units = "in", width = 6, height = 4, res = 300)
grid.arrange(p1, p2, p3, p4, nrow=2, ncol=2, top = "Proportion of Errors Among Pathologists")
dev.off()


# plotting the results for each career level

p1 <- pf2 %>% filter(!is.na(wrong)) %>%
  ggplot(aes(as.factor(Plevels), fill= as.factor(wrong))) +
  geom_bar(position = "fill", stat = "count") +
  labs (x="",
        y="Proportion of cases")+
  scale_x_discrete(breaks = c("1","2","3","4"), 
                   labels = c("Resident", "Junior","Expert", "Senior"))+
  scale_y_continuous(breaks=seq(from=0, to=1, by=0.25), 
                     labels = c("0%","25%","50%","75%","100%"))+
  scale_fill_manual(name="", breaks = c("1","0"), 
                    labels = c("Wrong", "Correct"), values=c("#de2d26","#2ca25f"))+
  theme(legend.position = "none")

p2<- pf2 %>% filter(!is.na(wrong)) %>%
  ggplot(aes(as.factor(Plevels), fill= as.factor(wrong))) +
  geom_bar(position = "fill", stat = "count") +
  labs (x="",
        y="")+
  scale_x_discrete(breaks = c("1","2","3","4"), 
                   labels = c("R", "J","S", "E"))+
  scale_y_continuous(breaks=seq(from=0, to=1, by=0.25), 
                     labels = c("0%","25%","50%","75%","100%"))+
  scale_fill_manual(name="", breaks = c("1","0"), 
                    labels = c("Wrong", "Correct"), values=c("#de2d26","#2ca25f"))+
  theme(legend.position = "none")+
  facet_wrap(vars(category))+ 
  theme(legend.position = "none")

Slevels.labs <- c("Surgery", "Biopsy")
names(Slevels.labs) <- c(0,1)

p3<- pf2 %>% filter(!is.na(wrong)) %>%
  ggplot(aes(as.factor(Plevels), fill= as.factor(wrong))) +
  geom_bar(position = "fill", stat = "count") +
  labs (x="Pathologist's Category",
        y="Proportion of cases")+
  scale_x_discrete(breaks = c("1","2","3","4"), 
                   labels = c("R", "J","E", "S"))+
  scale_y_continuous(breaks=seq(from=0, to=1, by=0.25), 
                     labels = c("0%","25%","50%","75%","100%"))+
  scale_fill_manual(name="", breaks = c("1","0"), 
                    labels = c("Wrong", "Correct"), values=c("#de2d26","#2ca25f"))+
  theme(legend.position = "none")+
  facet_wrap(vars(biopsy), labeller = labeller(biopsy = Slevels.labs)) 

p4<- pf2 %>% filter(!is.na(wrong)) %>%
  ggplot(aes(as.factor(Plevels), fill= as.factor(wrong))) +
  geom_bar(position = "fill", stat = "count") +
  labs (x="Pathologist's Category",
        y="")+
  scale_x_discrete(breaks = c("1","2","3","4"), 
                   labels = c("R", "J","E", "S"))+
  scale_y_continuous(breaks=seq(from=0, to=1, by=0.25), 
                     labels = c("0%","25%","50%","75%","100%"))+
  scale_fill_manual(name="", breaks = c("1","0"), 
                    labels = c("Wrong", "Correct"), values=c("#de2d26","#2ca25f"))+
  theme(legend.position = "none")+
  facet_wrap(vars(spec))

png("/Users/lollo/Dropbox/SmartCovid/paper/3 JIMR/R1/Supp_PLev.png", 
    units = "in", width = 6, height = 4, res = 300)
grid.arrange(p1, p2, p3, p4, nrow=2, ncol=2, top = "Proportion of Errors Among Pathologist's Category")
dev.off()

########################################################################################################
#Modeling
#list preparation for Stan
i <- !is.na(pf2$wrong)
dat <- list( 
  W = as.integer(pf2$wrong[i]),
  I = as.integer(pf2$PID[i]),
  L = as.integer(pf2$Plevels[i]),
  C = as.integer(pf2$category[i]),
  B = as.integer(pf2$biopsy[i]+1),
  S = as.integer(pf2$spec[i]))
#multilevel model including the pathologist ID, the career level, the category in the questionnaire
#the kind of specimen and the subspecialty
# pf.m1 <- ulam(
#   alist(
#     W ~ dbinom( 1 , p ) ,
#     logit(p) <- a[I] + b[L] + g[C] + d[B] + e[S], 
#     a[I] ~ dnorm( a_bar , sigma_a ),
#     b[L] ~ dnorm( 0 , sigma_b ),
#     g[C] ~ dnorm( 0 , sigma_g ),
#     d[B] ~ dnorm( 0 , sigma_d ),
#     e[S] ~ dnorm( 0 , sigma_e ), 
#     a_bar ~ dnorm( 0 , 1.5 ),
#     sigma_a ~ dexp(1), 
#     sigma_b ~ dexp(1),
#     sigma_g ~ dexp(1),
#     sigma_d ~ dexp(1),
#     sigma_e ~ dexp(1)
#     ), data=dat , chains=4 , cores=4 )

#the same model with non-centered parameter to deal with divergent transitions
pf.m1.nc <- ulam(
  alist(
    W ~ dbinom( 1 , p ) ,
    logit(p) <- a_bar + a[I]*sigma_a + b[L]*sigma_b + 
      g[C]*sigma_g + d[B]*sigma_d + e[S]*sigma_e, 
    a[I] ~ dnorm( 0 , 1 ),
    b[L] ~ dnorm( 0 , 1 ),
    g[C] ~ dnorm( 0 , 1 ),
    d[B] ~ dnorm( 0 , 1 ),
    e[S] ~ dnorm( 0 , 1 ), 
    a_bar ~ dnorm( 0 , 1.5 ),
    sigma_a ~ dexp(1), 
    sigma_b ~ dexp(1),
    sigma_g ~ dexp(1),
    sigma_d ~ dexp(1),
    sigma_e ~ dexp(1)
  ), data=dat , chains=4 , cores=4, log_lik = TRUE)

#"hairy caterpillar ocular inspection test"
traceplot(pf.m1.nc) #saved with Rstudio
#model diagnostics; supplementary table 1
precis(pf.m1.nc,2)

#overfitting control, supplementary table 2
WAIC(pf.m1.nc)
LOO(pf.m1.nc)

#prior predictive simulation
prior <- extract.prior(pf.m1.nc)

#displaying the prior predictive simulation
precis(as.data.frame(prior))
png("/Users/lollo/Dropbox/SmartCovid/paper/3 JIMR/R1/SuppPriorcoeff.png", units = "in", width = 5, height = 8, res = 300)
plot(precis(as.data.frame(prior)), labels=c("Pathologist #1","Pathologist #2", "Pathologist #3", "Pathologist #4", "Pathologist #5", 
                                           "Pathologist #6", "Pathologist #7", "Pathologist #8", "Pathologist #9",
                                           "Pathologist #10","Pathologist #11","Pathologist #12","Pathologist #13",
                                           "Pathologist #14","Pathologist #15","Pathologist #16","Pathologist #17",
                                           "Level, Resident", "Level, Junior","Level, Expert", "Level, Senior",
                                           "Category, neoplastic VS not","Category, benign VS malignant","Category, diagnosis",
                                           "Category, histotype", "Category, grade", "Specimen, surgical",
                                           "Specimen, bioptic", "Specialty, breast", "Specialty, GI", 
                                           "Specialty, uro", "alpha_bar (avg pathologist)", "sigma_a (pathologist)", "sigma_b (level)", 
                                           "sigma_g (category)", "sigma_d (specimen)", "sigma_e (specialty)"),
     main = "Prior Coefficients", xlab ="")
abline(v=0, h= c(6.5,9.5,11.5,16.5,20.5), lty = 1)
dev.off()

#prior predictive simulation  for Pathologists' levels
p_link_abar <- function( L ) {
  logodds <- with( prior , a_bar + b[,L] ) 
  return( inv_logit(logodds) )
}
p_raw <- sapply( 1:4 , function(i) p_link_abar( i ) ) 
p_mu <- apply( p_raw , 2 , mean )
p_ci <- apply( p_raw , 2 , HPDI )
wr <- by( dat$W,dat$L , mean, simplify = FALSE )

png("/Users/lollo/Dropbox/SmartCovid/paper/3 JIMR/R1/SuppPriorSimAvgPath.png", units = "in", width = 5, height = 5, res = 300)
plot( NULL , xlab="Pathologist's Category" , ylab="Proportion of error" ,
      ylim=c(0,1) , xaxt="n", yaxt="n" , xlim=c(1,4), main = "Prior Predictive Simulation of Average Pathologist" )
axis( 1 , at=1:4 , labels=c("Resident","Junior","Expert","Senior") )
axis(2, at = seq(from=0, to=1, by=0.25), 
     labels = c("0%","25%","50%","75%","100%"))
lines( 1:4 , p_mu )
shade( p_ci , 1:4 )
dev.off()

#plotting model coefficients
png("/Users/lollo/Dropbox/SmartCovid/paper/3 JIMR/R1/model1_coeff.png", units = "in", width = 5, height = 8, res = 300)
plot(precis(pf.m1.nc,3), labels = c("Pathologist #1","Pathologist #2", "Pathologist #3", "Pathologist #4", "Pathologist #5", 
                                    "Pathologist #6", "Pathologist #7", "Pathologist #8", "Pathologist #9",
                                    "Pathologist #10","Pathologist #11","Pathologist #12","Pathologist #13",
                                    "Pathologist #14","Pathologist #15","Pathologist #16","Pathologist #17",
                                    "Level, Resident", "Level, Junior","Level, Expert", "Level, Senior",
                                    "Category, neoplastic VS not","Category, benign VS malignant","Category, diagnosis",
                                    "Category, histotype", "Category, grade", "Specimen, surgical",
                                    "Specimen, bioptic", "Specialty, breast", "Specialty, GI", 
                                    "Specialty, uro", "alpha_bar (avg pathologist)", "sigma_a (pathologist)", "sigma_b (level)", 
                                    "sigma_g (category)", "sigma_d (specimen)", "sigma_e (specialty)"),
     main = "Model Coefficients", xlab ="")
abline(v=0, h= c(6.5,9.5,11.5,16.5,20.5), lty = 1)
dev.off()

#postcheck for Pathologists' levels
post <- extract.samples(pf.m1.nc)
p_link_abar <- function( L ) {
  logodds <- with( post , a_bar + b[,L] ) 
  return( inv_logit(logodds) )
}
p_raw <- sapply( 1:4 , function(i) p_link_abar( i ) ) 
p_mu <- apply( p_raw , 2 , mean )
p_ci <- apply( p_raw , 2 , HPDI )
wr <- by( dat$W,dat$L , mean, simplify = FALSE )

png("/Users/lollo/Dropbox/SmartCovid/paper/3 JIMR/R1/Fig3.png", units = "in", width = 5, height = 5, res = 300)
plot( NULL , xlab="Pathologist's Category" , ylab="Proportion of error" ,
      ylim=c(0,1) , xaxt="n", yaxt="n" , xlim=c(1,4), main = "Average Pathologist" )
axis( 1 , at=1:4 , labels=c("Resident","Junior","Expert","Senior") )
axis(2, at = seq(from=0, to=1, by=0.25), 
     labels = c("0%","25%","50%","75%","100%"))
lines( 1:4 , p_mu )
shade( p_ci , 1:4 )
lines(1:4,wr,lty=3)
dev.off()

######################################################################################################àà
#psycology part data wrangling
psy2 <- psyco %>% gather(key = "item.n", value = "score", "c.1":"s.6")
psy2$item.n <- as.factor(psy2$item.n)
psy2$item <- ifelse(psy2$item.n %in% c("a.1", "a.2", "a.3", "a.4"),"attitude", 
                    ifelse(psy2$item.n %in% c("c.1", "c.2", "c.3", "c.4", "c.5", "c.6", "c.7"), "confidence",
                           ifelse(psy2$item.n %in% c("s.1", "s.2", "s.3", "s.4", "s.5", "s.6"), "satisfaction", 
                                  "ERR")))
psy2$item <- as.factor(psy2$item)
psy2$PID <- factor(psy2$PID, c("P1","P2", "P3", "P4", "P5", "P6", "P7", "P8", "P9",
                             "P10","P11","P12","P13","P14","P15","P16","P17"))

psy2$er_rate <- psy2$W/psy2$T

#data exploration
Plevels.labs <- c("resident", "junior", "expert", "senior")
names(Plevels.labs) <- 1:4
png("/Users/lollo/Dropbox/SmartCovid/paper/3 JIMR/R1/Fig4.png", units = "in", 
     width = 5, height = 4, res = 300)
psy2 %>% 
  ggplot(aes(x=as.factor(score), y=..prop.., group=item))+
  stat_count() + 
  geom_bar()+
  facet_grid(Plevels ~ item, 
             labeller = labeller(Plevels = Plevels.labs))+
  labs(x="",y ="")+
  scale_y_continuous(labels = scales::percent)+
  scale_x_discrete(breaks = 1:5, labels = c("Strongly Disagree","Moderately Disagree",
                                            "Neutral/ I don't know","Moderately Agree","Strongly Agree"))+
  theme(axis.text.x = element_text(angle = 45,hjust = 1))+
  ggtitle("Response to survey by Domain and Career Level")
dev.off()

#alternative version required by reviewer, discarded
# psy2 %>% 
#   ggplot(aes(x=as.factor(score), y=..prop.., group= item, fill=item))+
#   geom_bar(position = "dodge" )+
#   facet_grid(rows = vars(Plevels), labeller = labeller(Plevels = Plevels.labs))+
#   labs(x="",y ="")+
#   scale_y_continuous(labels = scales::percent)+
#   scale_x_discrete(breaks = 1:5, labels = c("Strongly Disagree",
#                                             "Moderately Disagree",
#                                             "Neutral/ I don't know",
#                                             "Moderately Agree",
#                                             "Strongly Agree"))+
#   theme(axis.text.x = element_text(angle = 45,hjust = 1)) +
#   scale_fill_brewer(palette="Set1") + 
#   labs(fill = "Domain")


