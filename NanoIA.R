setwd("C:/Users/luisr/Documentos/Documents/Documentos - Copia/NanoIA")

#Simplex Design

library(ggplot2)
library(mixexp)

?MixModel

install.packages("lmtest")
library(lmtest)

install.packages("MASS")
library(MASS)

install.packages("sandwich")
library(sandwich)

install.packages("EnvStats")
library(EnvStats)

library(performance)
check_model(ucuubaandiroba_PDI_lin)
check_model(ucuubaandiroba_PDI_fullcub)

library(DHARMa)
plot(simulateResiduals(ucuubaandiroba_PDI_lin))
plot(simulateResiduals(ucuubaandiroba_PDI_fullcub))

install.packages("qqplotr")
library(qqplotr)

install.packages("car")
library(car)

simplex_design <- read.csv("proportions.csv")
head(simplex_design)


simplex_plot <- ggtern(
                  data = simplex_design,
                  aes(x = Lipid.1, y = Lipid.3, z = Lipid.2),
                ) +
                stat_density_tern(geom = "polygon",
                                  n = 400,
                                  aes(fill = ..level..,
                                      alpha = ..level..)) +
  scale_fill_gradient(low = "blue",
                      high = "red",
                      name = "Density curve",
                      breaks = c(5,10, 15, 20),
                      labels = c("Low", "", "", "High")) +
  scale_L_continuous(breaks = 0:5 / 5, labels = 0:5/ 5) +
  scale_R_continuous(breaks = 0:5 / 5, labels = 0:5/ 5) +
  scale_T_continuous(breaks = 0:5 / 5, labels = 0:5/ 5) +
  geom_point(
    shape = 21,
    size = 3,
    color = "black",
    fill = "#592693",
    stroke = 1.5
  ) +
                theme_rgbg() +
                theme_showarrows() +
                ggtitle("Simplex Design") +
                xlab("Lipid 1") +
                ylab("Surfactant") +
                zlab("Lipid 2") +
                theme(
                  plot.title = element_text(hjust = 0.5,
                                            face = "bold"),
                  tern.axis.arrow.text.T = element_text(vjust = 0),
                  tern.axis.arrow.text.L = element_text(vjust = 0),
                  tern.axis.arrow.text.R = element_text(vjust = 1)
                ) +
  guides(fill = guide_colorbar(order = 1), alpha = guide_none()) 
              

simplex_plot

ggsave(
  "Simplex_Plot.png",
  plot = simplex_plot,
  bg = "transparent"
)

?ggtitle
?plot.title
?theme
?ggsave
?MixModel

simplex_plot

#Nano murumuru-cupuaçu

murucupuacu <- read.csv("NanoMurumuruCupuacu.csv")

#PDI
murucupuacu_lin <- MixModel(frame = murucupuacu, 
                            response = "PDI",
                            mixcomps = c("Murumuru", "Cupuaçu", "Brijo10"),
                            model = 1)
summary(murucupuacu_lin)

plot(murucupuacu_lin)

anovaPE(murucupuacu_lin)

check_model(murucupuacu_lin)

plot(simulateResiduals(murucupuacu_lin))

murucupuacu_quad <- MixModel(frame = murucupuacu, 
                            response = "PDI",
                            mixcomps = c("Murumuru", "Cupuaçu", "Brijo10"),
                            model = 2)
summary(murucupuacu_quad)

check_model(murucupuacu_quad)

plot(murucupuacu_quad)

murucupuacu_fullcub <- MixModel(frame = murucupuacu, 
                            response = "PDI",
                            mixcomps = c("Murumuru", "Cupuaçu", "Brijo10"),
                            model = 3)
summary(murucupuacu_fullcub)

murucupuacu_splcub <- MixModel(frame = murucupuacu, 
                            response = "PDI",
                            mixcomps = c("Murumuru", "Cupuaçu", "Brijo10"),
                            model = 4)
summary(murucupuacu_splcub)

plot(murucupuacu_splcub)

check_model(murucupuacu_splcub)


AIC(murucupuacu_lin, murucupuacu_quad, murucupuacu_splcub, murucupuacu_fullcub)
BIC(murucupuacu_lin, murucupuacu_quad, murucupuacu_splcub, murucupuacu_fullcub)

anova(murucupuacu_lin, murucupuacu_splcub)

?bgtest
bgtest(murucupuacu_splcub)



plot(predict(murucupuacu_splcub), rstudent(murucupuacu_splcub),
     xlab=expression( paste("predicted response ", hat(y)) ) ,ylim=c(-5,5)   ,
     ylab="studentized residuals",type="n",
     main=expression(paste("(",frac(paste("y - ", hat(y)),sigma), ") ~ ", hat(y)   )) ) 

text(predict(murucupuacu_splcub), rstudent(murucupuacu_splcub),label=1:length(murucupuacu$PDI)) 
abline(h=c(-2,0,2),lty=c(2,1,2)) 
grid() 

murucupuacu1 <- murucupuacu

murucupuacu1$PDI[murucupuacu1$PDI == 0.458] <- NA

murucupuacu1_splcub <- MixModel(frame = murucupuacu1, 
                               response = "PDI",
                               mixcomps = c("Murumuru", "Cupuaçu", "Brijo10"),
                               model = 4)
summary(murucupuacu1_splcub)

plot(predict(murucupuacu1_splcub), rstudent(murucupuacu1_splcub),
     xlab=expression( paste("predicted response ", hat(y)) ) ,ylim=c(-5,5)   ,
     ylab="studentized residuals",type="n",
     main=expression(paste("(",frac(paste("y - ", hat(y)),sigma), ") ~ ", hat(y)   )) ) 

text(predict(murucupuacu1_splcub), rstudent(murucupuacu1_splcub),label=1:length(murucupuacu1$PDI)) 
abline(h=c(-2,0,2),lty=c(2,1,2)) 
grid() 

ModelPlot(model = murucupuacu1_splcub,
          res = 500,
          at = c(1,0.5,0.4,0.35,0.3,0.25,0.2,0.15,0.1),
          dimensions = list(x1 = "Murumuru", x2 = "Cupuaçu", x3 = "Brijo10"),
          lims = c(0.166, 0.5, 0.166, 0.5, 0.166, 0.5835),
          contour = TRUE, fill = TRUE, pseudo = TRUE, axislab.offset = 0.1,
          axislabs = c("Murumuru", "Cupuaçu", "Brij O10"),
          cornerlabs = c("", "", ""),
          grid = TRUE, grid.pars = list(lty = 0),
          main = "PDI Contour Plot")

#size

murucupuacu1 <- murucupuacu1[-12, ]

murucupuacu_size_lin <- MixModel(frame = murucupuacu, 
                            response = "Size",
                            mixcomps = c("Murumuru", "Cupuaçu", "Brijo10"),
                            model = 1)
summary(murucupuacu_size_lin)

plot(murucupuacu_size_lin)

check_model(murucupuacu_size_lin)

murucupuacu_size_quad <- MixModel(frame = murucupuacu, 
                             response = "Size",
                             mixcomps = c("Murumuru", "Cupuaçu", "Brijo10"),
                             model = 2)
summary(murucupuacu_size_quad)

plot(murucupuacu_size_quad)

check_model(murucupuacu_size_quad)

murucupuacu_size_fullcub <- MixModel(frame = murucupuacu, 
                                response = "Size",
                                mixcomps = c("Murumuru", "Cupuaçu", "Brijo10"),
                                model = 3)
summary(murucupuacu_size_fullcub)

murucupuacu_size_splcub <- MixModel(frame = murucupuacu, 
                               response = "Size",
                               mixcomps = c("Murumuru", "Cupuaçu", "Brijo10"),
                               model = 4)
summary(murucupuacu_size_splcub)

plot(predict(murucupuacu_size_quad), rstudent(murucupuacu_size_quad),
     xlab=expression( paste("predicted response ", hat(y)) ) ,ylim=c(-5,5)   ,
     ylab="studentized residuals",type="n",
     main=expression(paste("(",frac(paste("y - ", hat(y)),sigma), ") ~ ", hat(y)   )) ) 

text(predict(murucupuacu_size_quad), rstudent(murucupuacu_size_quad),label=1:length(murucupuacu$Size)) 
abline(h=c(-2,0,2),lty=c(2,1,2)) 
grid() 

AIC(murucupuacu_size_lin, murucupuacu_size_quad, murucupuacu_size_splcub, 
    murucupuacu_size_fullcub)
BIC(murucupuacu_size_lin, murucupuacu_size_quad, murucupuacu_size_splcub, 
    murucupuacu_size_fullcub)
anova(murucupuacu_size_lin, murucupuacu_size_quad)
anova(murucupuacu_size_quad, murucupuacu_size_splcub)

bgtest(murucupuacu_size_quad)


murucupuacu1.1 <- murucupuacu1

murucupuacu1.1 <- murucupuacu1.1[-14, ]

murucupuacu1.1_size_quad <- MixModel(frame = murucupuacu1.1, 
                                response = "Size",
                                mixcomps = c("Murumuru", "Cupuaçu", "Brijo10"),
                                model = 2)
summary(murucupuacu1.1_size_quad)

plot(predict(murucupuacu1.1_size_quad), rstudent(murucupuacu1.1_size_quad),
     xlab=expression( paste("predicted response ", hat(y)) ) ,ylim=c(-5,5)   ,
     ylab="studentized residuals",type="n",
     main=expression(paste("(",frac(paste("y - ", hat(y)),sigma), ") ~ ", hat(y)   )) ) 

text(predict(murucupuacu1.1_size_quad), rstudent(murucupuacu1.1_size_quad),label=1:length(murucupuacu1.1$Size)) 
abline(h=c(-2,0,2),lty=c(2,1,2)) 
grid() 

ModelPlot(model = murucupuacu1.1_size_quad,
          res = 500,
          dimensions = list(x1 = "Murumuru", x2 = "Cupuaçu", x3 = "Brijo10"),
          lims = c(0.166, 0.5, 0.166, 0.5, 0.166, 0.5835),
          contour = TRUE, fill = TRUE, pseudo = TRUE, axislab.offset = 0.1,
          axislabs = c("Murumuru", "Cupuaçu", "Brij O10"),
          cornerlabs = c("", "", ""),
          grid = TRUE, grid.pars = list(lty = 0),
          main = "Size Contour Plot")
#Nano bacuri açaí

bacuriacai <- read.csv("NanoBacuriAcai.csv")

#PDI
bacuriacai_PDI_lin <- MixModel(frame = bacuriacai, 
                            response = "PDI",
                            mixcomps = c("Bacuri", "Acai", "Brijo10"),
                            model = 1)
summary(bacuriacai_PDI_lin)

bacuriacai_PDI_quad <- MixModel(frame = bacuriacai, 
                             response = "PDI",
                             mixcomps = c("Bacuri", "Acai", "Brijo10"),
                             model = 2)
summary(bacuriacai_PDI_quad)

bacuriacai_PDI_fullcub <- MixModel(frame = bacuriacai, 
                                response = "PDI",
                                mixcomps = c("Bacuri", "Acai", "Brijo10"),
                                model = 3)
summary(bacuriacai_PDI_fullcub)

bacuriacai_PDI_splcub <- MixModel(frame = bacuriacai, 
                               response = "PDI",
                               mixcomps = c("Bacuri", "Acai", "Brijo10"),
                               model = 4)
summary(bacuriacai_PDI_splcub)

AIC(bacuriacai_PDI_lin, bacuriacai_PDI_quad, bacuriacai_PDI_splcub, 
    bacuriacai_PDI_fullcub)
BIC(bacuriacai_PDI_lin, bacuriacai_PDI_quad, bacuriacai_PDI_splcub, 
    bacuriacai_PDI_fullcub)
anova(bacuriacai_PDI_lin, bacuriacai_PDI_fullcub)

bgtest(bacuriacai_PDI_lin)

plot(predict(bacuriacai_PDI_lin), rstudent(bacuriacai_PDI_lin),
     xlab=expression( paste("predicted response ", hat(y)) ) ,ylim=c(-5,5)   ,
     ylab="studentized residuals",type="n",
     main=expression(paste("(",frac(paste("y - ", hat(y)),sigma), ") ~ ", hat(y)   )) ) 

text(predict(bacuriacai_PDI_lin), rstudent(bacuriacai_PDI_lin),label=1:length(bacuriacai$PDI)) 
abline(h=c(-2,0,2),lty=c(2,1,2)) 
grid() 


ModelPlot(model = bacuriacai_PDI_lin,
          res = 500,
          dimensions = list(x1 = "Bacuri", x2 = "Acai", x3 = "Brijo10"),
          lims = c(0.166, 0.5, 0.166, 0.5, 0.166, 0.5835),
          at = c(0, 0.2, 0.4, 0.5, 0.6, 0.8, 1, 2),
          contour = TRUE, fill = TRUE, pseudo = TRUE, axislab.offset = 0.1,
          axislabs = c("Bacuri", "Açaí", "Brij O10"),
          cornerlabs = c("", "", ""),
          grid = TRUE, grid.pars = list(lty = 0),
          main = "PDI Contour Plot")

#size

bacuriacai_size_lin <- MixModel(frame = bacuriacai, 
                                 response = "Size",
                                 mixcomps = c("Bacuri", "Acai", "Brijo10"),
                                 model = 1)
summary(bacuriacai_size_lin)

bacuriacai_size_quad <- MixModel(frame = bacuriacai, 
                                  response = "Size",
                                  mixcomps = c("Bacuri", "Acai", "Brijo10"),
                                  model = 2)
summary(bacuriacai_size_quad)

bacuriacai_size_fullcub <- MixModel(frame = bacuriacai, 
                                     response = "Size",
                                     mixcomps = c("Bacuri", "Acai", "Brijo10"),
                                     model = 3)
summary(bacuriacai_size_fullcub)

bacuriacai_size_splcub <- MixModel(frame = bacuriacai, 
                                    response = "Size",
                                    mixcomps = c("Bacuri", "Acai", "Brijo10"),
                                    model = 4)
summary(bacuriacai_size_splcub)

AIC(bacuriacai_size_lin, bacuriacai_size_quad, bacuriacai_size_splcub, 
    bacuriacai_size_fullcub)
?AIC
anova(bacuriacai_size_lin, bacuriacai_size_fullcub)
anova(bacuriacai_size_lin, bacuriacai_size_splcub)
anova(bacuriacai_size_lin, bacuriacai_size_quad)

bgtest(bacuriacai_size_lin)

plot(predict(bacuriacai_size_lin), rstudent(bacuriacai_size_lin),
     xlab=expression( paste("predicted response ", hat(y)) ) ,ylim=c(-5,5)   ,
     ylab="studentized residuals",type="n",
     main=expression(paste("(",frac(paste("y - ", hat(y)),sigma), ") ~ ", hat(y)   )) ) 

text(predict(bacuriacai_size_lin), rstudent(bacuriacai_size_lin),label=1:length(bacuriacai$Size)) 
abline(h=c(-2,0,2),lty=c(2,1,2)) 
grid() 

ModelPlot(model = bacuriacai_size_splcub,
          res = 500,
          dimensions = list(x1 = "Bacuri", x2 = "Acai", x3 = "Brijo10"),
          at = c(-1000, 100, 200, 300, 400, 500, 600, 700, 800, 900),
          lims = c(0.166, 0.5, 0.166, 0.5, 0.166, 0.5835),
          contour = TRUE, fill = TRUE, pseudo = TRUE, axislab.offset = 0.1,
          axislabs = c("Bacuri", "Açaí", "Brij O10"),
          cornerlabs = c("", "", ""),
          grid = TRUE, grid.pars = list(lty = 0),
          main = "Size Contour Plot")

#Nano bacuri pracaxi

bacuripracaxi <- read.csv("NanoBacuriPracaxi.csv")

#PDI
bacuripracaxi_PDI_lin <- MixModel(frame = bacuripracaxi, 
                               response = "PDI",
                               mixcomps = c("Bacuri", "Pracaxi", "Brijo10"),
                               model = 1)
summary(bacuripracaxi_PDI_lin)

bacuripracaxi_PDI_quad <- MixModel(frame = bacuripracaxi, 
                                response = "PDI",
                                mixcomps = c("Bacuri", "Pracaxi", "Brijo10"),
                                model = 2)
summary(bacuripracaxi_PDI_quad)

bacuripracaxi_PDI_fullcub <- MixModel(frame = bacuripracaxi, 
                                   response = "PDI",
                                   mixcomps = c("Bacuri", "Pracaxi", "Brijo10"),
                                   model = 3)
summary(bacuripracaxi_PDI_fullcub)

bacuripracaxi_PDI_splcub <- MixModel(frame = bacuripracaxi, 
                                  response = "PDI",
                                  mixcomps = c("Bacuri", "Pracaxi", "Brijo10"),
                                  model = 4)
summary(bacuripracaxi_PDI_splcub)

AIC(bacuripracaxi_PDI_lin, bacuripracaxi_PDI_quad, 
    bacuripracaxi_PDI_splcub, bacuripracaxi_PDI_fullcub)

anova(bacuripracaxi_PDI_lin, bacuripracaxi_PDI_quad)

bgtest(bacuripracaxi_PDI_lin)

bgtest(bacuripracaxi_PDI_quad)

plot(predict(bacuripracaxi_PDI_lin), rstudent(bacuripracaxi_PDI_lin),
     xlab=expression( paste("predicted response ", hat(y)) ) ,ylim=c(-5,5)   ,
     ylab="studentized residuals",type="n",
     main=expression(paste("(",frac(paste("y - ", hat(y)),sigma), ") ~ ", hat(y)   )) ) 

text(predict(bacuripracaxi_PDI_lin), rstudent(bacuripracaxi_PDI_lin),label=1:length(bacuripracaxi$PDI)) 
abline(h=c(-2,0,2),lty=c(2,1,2)) 
grid() 

plot(bacuripracaxi_PDI_lin$fitted.values, bacuripracaxi_PDI_lin$residuals)
abline(h=0)

?boxcox
boxcox(bacuriacai_PDI_lin)

?vcov
?vcovHC

?coeftest

coeftest(bacuripracaxi_PDI_lin, vcov. = vcovHC(bacuripracaxi_PDI_lin,
                                               type = "HC3"))

coeftest(bacuripracaxi_PDI_quad, vcov. = vcovHC(bacuripracaxi_PDI_quad,
                                               type = "HC1"))

ModelPlot(model = bacuripracaxi_PDI_lin,
          res = 500,
          dimensions = list(x1 = "Bacuri", x2 = "Pracaxi", x3 = "Brijo10"),
          lims = c(0.166, 0.5, 0.166, 0.5, 0.166, 0.5835),
          at = c(-100, 0.1, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1, 2, 3, 100),
          contour = TRUE, fill = TRUE, pseudo = TRUE, axislab.offset = 0.15,
          axislabs = c("Bacuri", "Pracaxi", "Brij O10"),
          cornerlabs = c("", "", ""),
          grid = TRUE, grid.pars = list(lty = 0),
          main = "PDI Contour Plot")

#size

bacuripracaxi_size_lin <- MixModel(frame = bacuripracaxi, 
                                response = "Size",
                                mixcomps = c("Bacuri", "Pracaxi", "Brijo10"),
                                model = 1)
summary(bacuripracaxi_size_lin)

bacuripracaxi_size_quad <- MixModel(frame = bacuripracaxi, 
                                 response = "Size",
                                 mixcomps = c("Bacuri", "Pracaxi", "Brijo10"),
                                 model = 2)
summary(bacuripracaxi_size_quad)

bacuripracaxi_size_fullcub <- MixModel(frame = bacuripracaxi, 
                                    response = "Size",
                                    mixcomps = c("Bacuri", "Pracaxi", "Brijo10"),
                                    model = 3)
summary(bacuripracaxi_size_fullcub)

bacuripracaxi_size_splcub <- MixModel(frame = bacuripracaxi, 
                                   response = "Size",
                                   mixcomps = c("Bacuri", "Pracaxi", "Brijo10"),
                                   model = 4)
summary(bacuripracaxi_size_splcub)

View(bacuripracaxi)

AIC(bacuripracaxi_size_lin, bacuripracaxi_size_quad, 
    bacuriacai_size_splcub, bacuriacai_size_fullcub)

BIC(bacuripracaxi_size_lin, bacuripracaxi_size_quad, 
    bacuriacai_size_splcub, bacuriacai_size_fullcub)

anova(bacuripracaxi_size_lin, bacuriacai_size_quad)
anova(bacuripracaxi_size_quad, bacuriacai_size_splcub)
anova(bacuripracaxi_size_lin, bacuripracaxi_size_fullcub)

bgtest(bacuriacai_size_fullcub)
bgtest(bacuripracaxi_size_quad)

plot(bacuripracaxi_size_fullcub$fitted.values, bacuripracaxi_size_fullcub$residuals)
abline(h=0)

coeftest(bacuriacai_size_fullcub, vcov. = vcovHC(bacuriacai_size_fullcub,
                                               type = "HC1"))

plot(predict(bacuripracaxi_size_quad), rstudent(bacuripracaxi_size_quad),
     xlab=expression( paste("predicted response ", hat(y)) ) ,ylim=c(-5,5)   ,
     ylab="studentized residuals",type="n",
     main=expression(paste("(",frac(paste("y - ", hat(y)),sigma), ") ~ ", hat(y)   )) ) 

text(predict(bacuripracaxi_size_quad), rstudent(bacuripracaxi_size_quad),label=1:length(bacuripracaxi$Size)) 
abline(h=c(-2,0,2),lty=c(2,1,2)) 
grid() 

ModelPlot(model = bacuripracaxi_size_quad,
          res = 500,
          dimensions = list(x1 = "Bacuri", x2 = "Pracaxi", x3 = "Brijo10"),
          lims = c(0.166, 0.5, 0.166, 0.5, 0.166, 0.5835),
          at = c(-30000, 100, 300, 500, 700, 800, 900, 1000, 30000),
          contour = TRUE, fill = TRUE, pseudo = TRUE, axislab.offset = 0.1,
          axislabs = c("Bacuri", "Pracaxi", "Brij O10"),
          cornerlabs = c("", "", ""),
          grid = TRUE, grid.pars = list(lty = 0),
          main = "Size Contour Plot")

#CacauTucumã

#PDI
setwd("C:/Users/luisr/Documentos/Documents/Documentos - Copia/NanoIA")
cacautucuma <- read.csv("NanoCacauTucuma.csv")

cacautucuma_PDI_lin <- MixModel(frame = cacautucuma, 
                                  response = "PDI",
                                  mixcomps = c("Cacau", "Tucuma", "Brijo10"),
                                  model = 1)
summary(cacautucuma_PDI_lin)

cacautucuma_PDI_quad <- MixModel(frame = cacautucuma, 
                                   response = "PDI",
                                   mixcomps = c("Cacau", "Tucuma", "Brijo10"),
                                   model = 2)
summary(cacautucuma_PDI_quad)

cacautucuma_PDI_fullcub <- MixModel(frame = cacautucuma, 
                                      response = "PDI",
                                    
                                      mixcomps = c("Cacau", "Tucuma", "Brijo10"),
                                      model = 3)
summary(cacautucuma_PDI_fullcub)

plot(cacautucuma_PDI_fullcub)

cacautucuma_PDI_splcub <- MixModel(frame = cacautucuma, 
                                     response = "PDI",
                                     mixcomps = c("Cacau", "Tucuma", "Brijo10"),
                                     model = 4)
summary(cacautucuma_PDI_splcub)

AIC(cacautucuma_PDI_lin, cacautucuma_PDI_quad, cacautucuma_PDI_splcub,
    cacautucuma_PDI_fullcub)
BIC(cacautucuma_PDI_lin, cacautucuma_PDI_quad, cacautucuma_PDI_splcub,
    cacautucuma_PDI_fullcub)

anova(cacautucuma_PDI_lin, cacautucuma_PDI_fullcub)

bgtest(cacautucuma_PDI_fullcub)
bgtest(cacautucuma_PDI_lin)

coeftest(cacautucuma_PDI_fullcub, vcov. = vcovHC(cacautucuma_PDI_fullcub,
                                                 type = "HC1"))



plot(predict(cacautucuma_PDI_fullcub), rstudent(cacautucuma_PDI_fullcub),
     xlab=expression( paste("predicted response ", hat(y)) ) ,ylim=c(-5,5)   ,
     ylab="studentized residuals",type="n",
     main=expression(paste("(",frac(paste("y - ", hat(y)),sigma), ") ~ ", hat(y)   )) ) 

text(predict(cacautucuma_PDI_fullcub), rstudent(cacautucuma_PDI_fullcub),label=1:length(cacautucuma$PDI)) 
abline(h=c(-2,0,2),lty=c(2,1,2)) 
grid() 


ModelPlot(model = cacautucuma_PDI_fullcub,
          res = 500,
          dimensions = list(x1 = "Cacau", x2 = "Tucuma", x3 = "Brijo10"),
          lims = c(0.166, 0.5, 0.166, 0.5, 0.166, 0.5835),
          at = c(-100, 0.1, 0.3, 0.5, 0.9, 100),
          contour = TRUE, fill = TRUE, pseudo = TRUE, axislab.offset = 0.15,
          color.palette = heat.colors,
          axislabs = c("Cacau", "Tucumã", "Brij O10"),
          cornerlabs = c("", "", ""),
          grid = TRUE, grid.pars = list(lty = 0),
          main = "PDI Contour Plot")

View(cacautucuma)
#Size
setwd("C:/Users/luisr/Documentos/Documents/Documentos - Copia/NanoIA")
cacautucuma <- read.csv("NanoCacauTucuma.csv")

cacautucuma_size_lin <- MixModel(frame = cacautucuma, 
                                response = "Size",
                                mixcomps = c("Cacau", "Tucuma", "Brijo10"),
                                model = 1)
summary(cacautucuma_size_lin)

cacautucuma_size_quad <- MixModel(frame = cacautucuma, 
                                 response = "Size",
                                 mixcomps = c("Cacau", "Tucuma", "Brijo10"),
                                 model = 2)
summary(cacautucuma_size_quad)

cacautucuma_size_fullcub <- MixModel(frame = cacautucuma, 
                                 response = "Size",
                                 mixcomps = c("Cacau", "Tucuma", "Brijo10"),
                                 model = 3)

summary(cacautucuma_size_fullcub)

plot(cacautucuma_size_fullcub)

cacautucuma_size_splcub <- MixModel(frame = cacautucuma, 
                                     response = "Size",
                                     mixcomps = c("Cacau", "Tucuma", "Brijo10"),
                                     model = 4)

summary(cacautucuma_size_splcub)

AIC(cacautucuma_size_lin, cacautucuma_size_quad,
    cacautucuma_size_splcub, cacautucuma_size_fullcub)

BIC(cacautucuma_size_lin, cacautucuma_size_quad,
    cacautucuma_size_splcub, cacautucuma_size_fullcub)

anova(cacautucuma_size_lin, cacautucuma_size_fullcub)

bgtest(cacautucuma_size_fullcub)

coeftest(cacautucuma_size_fullcub, vcov. = vcovHC(cacautucuma_size_fullcub,
                                                 type = "HC3"))

plot(predict(cacautucuma_size_fullcub), rstudent(cacautucuma_size_fullcub),
     xlab=expression( paste("predicted response ", hat(y)) ) ,ylim=c(-5,5)   ,
     ylab="studentized residuals",type="n",
     main=expression(paste("(",frac(paste("y - ", hat(y)),sigma), ") ~ ", hat(y)   )) ) 

text(predict(cacautucuma_size_fullcub), rstudent(cacautucuma_size_fullcub),label=1:length(cacautucuma$size)) 
abline(h=c(-2,0,2),lty=c(2,1,2)) 
grid() 

ModelPlot(model = cacautucuma_size_fullcub,
          res = 500,
          dimensions = list(x1 = "Cacau", x2 = "Tucuma", x3 = "Brijo10"),
          lims = c(0.166, 0.5, 0.166, 0.5, 0.166, 0.5835),
          at = c(-1000, 50, 100, 150, 200, 300, 400, 600, 800, 10000),
          contour = TRUE, fill = TRUE, pseudo = TRUE, axislab.offset = 0.1,
          axislabs = c("Cacau", "Tucumã", "Brij O10"),
          cornerlabs = c("", "", ""),
          grid = TRUE, grid.pars = list(lty = 0),
          main = "Size Contour Plot")

#Cupuaçu-Açaí

#PDI
setwd("C:/Users/luisr/Documentos/Documents/Documentos - Copia/NanoIA")
cupuacuacai <- read.csv("NanoCupuacuAcai.csv")

cupuacuacai_PDI_lin <- MixModel(frame = cupuacuacai, 
                                response = "PDI",
                                mixcomps = c("Cupuacu", "Acai", "Brijo10"),
                                model = 1)
summary(cupuacuacai_PDI_lin)

anova(cupuacuacai_PDI_lin)

plot(cupuacuacai_PDI_lin)

cupuacuacai_PDI_quad <- MixModel(frame = cupuacuacai, 
                                 response = "PDI",
                                 mixcomps = c("Cupuacu", "Acai", "Brijo10"),
                                 model = 2)
summary(cupuacuacai_PDI_quad)

plot(cupuacuacai_PDI_quad)


cupuacuacai_PDI_fullcub <- MixModel(frame = cupuacuacai, 
                                    response = "PDI",
                                    mixcomps = c("Cupuacu", "Acai", "Brijo10"),
                                    model = 3)
summary(cupuacuacai_PDI_fullcub)

cupuacuacai_PDI_splcub <- MixModel(frame = cupuacuacai, 
                                   response = "PDI",
                                   mixcomps = c("Cupuacu", "Acai", "Brijo10"),
                                   model = 4)
summary(cupuacuacai_PDI_splcub)

plot(cupuacuacai_PDI_splcub)

AIC(cupuacuacai_PDI_lin, cupuacuacai_PDI_quad, cupuacuacai_PDI_splcub, 
    cupuacuacai_PDI_fullcub)
BIC(cupuacuacai_PDI_lin, cupuacuacai_PDI_quad, cupuacuacai_PDI_splcub, 
    cupuacuacai_PDI_fullcub)

anova(cupuacuacai_PDI_lin, cupuacuacai_PDI_quad)

bgtest(cupuacuacai_PDI_lin)

plot(cupuacuacai_PDI_lin$fitted.values, cupuacuacai_PDI_lin$residuals)

coeftest(cupuacuacai_PDI_lin, vcov. = vcovHC(cupuacuacai_PDI_lin,
                                                  type = "HC1"))


plot(predict(cupuacuacai_PDI_lin), rstudent(cupuacuacai_PDI_lin),
     xlab=expression( paste("predicted response ", hat(y)) ) ,ylim=c(-5,5)   ,
     ylab="studentized residuals",type="n",
     main=expression(paste("(",frac(paste("y - ", hat(y)),sigma), ") ~ ", hat(y)   )) ) 

text(predict(cupuacuacai_PDI_lin), rstudent(cupuacuacai_PDI_lin),label=1:length(cupuacuacai$PDI)) 
abline(h=c(-2,0,2),lty=c(2,1,2)) 
grid() 


ModelPlot(model = cupuacuacai_PDI_lin,
          res = 500,
          dimensions = list(x1 = "Cupuacu", x2 = "Acai", x3 = "Brijo10"),
          lims = c(0.166, 0.5, 0.166, 0.5, 0.166, 0.5835),
          at = c(-100, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 100),
          contour = TRUE, fill = TRUE, pseudo = T, axislab.offset = 0.1,
          color.palette = heat.colors,
          axislabs = c("Cupuaçu", "Açaí", "Brij O10"),
          cornerlabs = c("", "", ""),
          grid = TRUE, grid.pars = list(lty = 0),
          main = "PDI Contour Plot")

#size

cupuacuacai_size_lin <- MixModel(frame = cupuacuacai, 
                                response = "Size",
                                mixcomps = c("Cupuacu", "Acai", "Brijo10"),
                                model = 1)
summary(cupuacuacai_size_lin)

plot(cupuacuacai_size_lin)

cupuacuacai_size_quad <- MixModel(frame = cupuacuacai, 
                                 response = "Size",
                                 mixcomps = c("Cupuacu", "Acai", "Brijo10"),
                                 model = 2)
summary(cupuacuacai_size_quad)

plot(cupuacuacai_size_quad)

cupuacuacai_size_fullcub <- MixModel(frame = cupuacuacai, 
                                    response = "Size",
                                    mixcomps = c("Cupuacu", "Acai", "Brijo10"),
                                    model = 3)
summary(cupuacuacai_size_fullcub)

cupuacuacai_size_splcub <- MixModel(frame = cupuacuacai, 
                                   response = "Size",
                                   mixcomps = c("Cupuacu", "Acai", "Brijo10"),
                                   model = 4)
summary(cupuacuacai_size_splcub)

AIC(cupuacuacai_size_lin, cupuacuacai_size_quad, 
    cupuacuacai_size_splcub, 
    cupuacuacai_size_fullcub)
BIC(cupuacuacai_size_lin, cupuacuacai_size_quad, 
    cupuacuacai_size_splcub, 
    cupuacuacai_size_fullcub)
anova(cupuacuacai_size_lin, cupuacuacai_size_fullcub)

bgtest(cupuacuacai_size_lin)

plot(predict(cupuacuacai_size_lin), rstudent(cupuacuacai_size_lin),
     xlab=expression( paste("predicted response ", hat(y)) ) ,ylim=c(-5,5)   ,
     ylab="studentized residuals",type="n",
     main=expression(paste("(",frac(paste("y - ", hat(y)),sigma), ") ~ ", hat(y)   )) ) 

text(predict(cupuacuacai_size_lin), rstudent(cupuacuacai_size_lin),label=1:length(cupuacuacai$Size)) 
abline(h=c(-2,0,2),lty=c(2,1,2)) 
grid() 


ModelPlot(model = cupuacuacai_size_lin,
          res = 500,
          dimensions = list(x1 = "Cupuacu", x2 = "Acai", x3 = "Brijo10"),
          lims = c(0.166, 0.5, 0.166, 0.5, 0.166, 0.5835),
          at = c(-10000, 100, 300, 500, 750, 1000, 9000),
          contour = TRUE, fill = TRUE, pseudo = T, axislab.offset = 0.1,
          color.palette = heat.colors,
          axislabs = c("Cupuaçu", "Açaí", "Brij O10"),
          cornerlabs = c("", "", ""),
          grid = TRUE, grid.pars = list(lty = 0),
          main = "Size Contour Plot")

#Cupuaçu-Cacau

cupuacucacau <- read.csv("NanoCupuacuCacau1.csv")
View(cupuacucacau)
summary(cupuacucacau)
class(cupuacucacau)

#PDI

cupuacucacau_PDI_lin <- MixModel(frame = cupuacucacau, 
                                response = "PDI",
                                mixcomps = c("Cupuacu", "Cacau", "Brijo10"),
                                model = 1)
summary(cupuacucacau_PDI_lin)

cupuacucacau_PDI_quad <- MixModel(frame = cupuacucacau, 
                                 response = "PDI",
                                 mixcomps = c("Cupuacu", "Cacau", "Brijo10"),
                                 model = 2)
summary(cupuacucacau_PDI_quad)

cupuacucacau_PDI_fullcub <- MixModel(frame = cupuacucacau, 
                                    response = "PDI",
                                    mixcomps = c("Cupuacu", "Cacau", "Brijo10"),
                                    model = 3)
summary(cupuacucacau_PDI_fullcub)

cupuacucacau_PDI_splcub <- MixModel(frame = cupuacucacau, 
                                   response = "PDI",
                                   mixcomps = c("Cupuacu", "Cacau", "Brijo10"),
                                   model = 4)
summary(cupuacucacau_PDI_splcub)

AIC(cupuacucacau_PDI_lin, cupuacucacau_PDI_quad, cupuacucacau_PDI_splcub)

anova(cupuacucacau_PDI_lin, cupuacucacau_PDI_quad)
anova(cupuacucacau_PDI_lin, cupuacucacau_PDI_splcub)

bgtest(cupuacucacau_PDI_lin)

coeftest(cupuacucacau_PDI_lin, vcov. = vcovHC(cupuacucacau_PDI_lin,
                                             type = "HC1"))

plot(predict(cupuacucacau_PDI_lin), rstudent(cupuacucacau_PDI_lin),
     xlab=expression( paste("predicted response ", hat(y)) ) ,ylim=c(-5,5)   ,
     ylab="studentized residuals",type="n",
     main=expression(paste("(",frac(paste("y - ", hat(y)),sigma), ") ~ ", hat(y)   )) ) 

text(predict(cupuacucacau_PDI_splcub), rstudent(cupuacucacau_PDI_splcub),
     label=1:length(cupuacucacau$PDI)) 
abline(h=c(-2,0,2),lty=c(2,1,2)) 
grid() 


ModelPlot(model = cupuacucacau_PDI_lin,
          res = 500,
          dimensions = list(x1 = "Cupuacu", x2 = "Cacau", x3 = "Brijo10"),
          lims = c(0.166, 0.5, 0.166, 0.5, 0.166, 0.5835),
          at = c(0.01, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 10),
          contour = TRUE, fill = TRUE, pseudo = T, axislab.offset = 0.1,
          color.palette = heat.colors,
          axislabs = c("Cupuaçu", "Cacau", "Brij O10"),
          cornerlabs = c("", "", ""),
          grid = TRUE, grid.pars = list(lty = 0),
          main = "PDI Contour Plot")

#size

cupuacucacau_size_lin <- MixModel(frame = cupuacucacau, 
                                 response = "Size",
                                 mixcomps = c("Cupuacu", "Cacau", "Brijo10"),
                                 model = 1)
summary(cupuacucacau_size_lin)

cupuacucacau_size_quad <- MixModel(frame = cupuacucacau, 
                                  response = "Size",
                                  mixcomps = c("Cupuacu", "Cacau", "Brijo10"),
                                  model = 2)
summary(cupuacucacau_size_quad)

cupuacucacau_size_fullcub <- MixModel(frame = cupuacucacau, 
                                     response = "Size",
                                     mixcomps = c("Cupuacu", "Cacau", "Brijo10"),
                                     model = 3)
summary(cupuacucacau_size_fullcub)

cupuacucacau_size_splcub <- MixModel(frame = cupuacucacau, 
                                    response = "Size",
                                    mixcomps = c("Cupuacu", "Cacau", "Brijo10"),
                                    model = 4)
summary(cupuacucacau_size_splcub)

cupuacucacau

AIC(cupuacucacau_size_lin, cupuacucacau_size_quad, cupuacucacau_size_splcub)

anova(cupuacucacau_size_fullcub, cupuacucacau_size_splcub)

bgtest(cupuacucacau_size_lin)

plot(predict(cupuacucacau_size_quad), rstudent(cupuacucacau_size_quad),
     xlab=expression( paste("predicted response ", hat(y)) ) ,ylim=c(-5,5)   ,
     ylab="studentized residuals",type="n",
     main=expression(paste("(",frac(paste("y - ", hat(y)),sigma), ") ~ ", hat(y)   )) ) 

text(predict(cupuacucacau_size_quad), rstudent(cupuacucacau_size_quad),label=1:length(cupuacucacau$Size)) 
abline(h=c(-2,0,2),lty=c(2,1,2)) 
grid() 


ModelPlot(model = cupuacucacau_size_lin,
          res = 500,
          dimensions = list(x1 = "Cupuacu", x2 = "Cacau", x3 = "Brijo10"),
          lims = c(0.166, 0.5, 0.166, 0.5, 0.166, 0.5835),
          at = c(-1000, 100, 300, 500, 700, 900, 10000),
          contour = TRUE, fill = TRUE, pseudo = T, axislab.offset = 0.1,
          color.palette = heat.colors,
          axislabs = c("Cupuaçu", "Cacau", "Brij O10"),
          cornerlabs = c("", "", ""),
          grid = TRUE, grid.pars = list(lty = 0),
          main = "Size Contour Plot")

#Cupuaçu-Pequi

cupuacupequi <- read.csv("NanoCupuacuPequi.csv")

#PDI

cupuacupequi_PDI_lin <- MixModel(frame = cupuacupequi, 
                                 response = "PDI",
                                 mixcomps = c("Cupuacu", "Pequi", "Brijo10"),
                                 model = 1)
summary(cupuacupequi_PDI_lin)

cupuacupequi_PDI_quad <- MixModel(frame = cupuacupequi, 
                                  response = "PDI",
                                  mixcomps = c("Cupuacu", "Pequi", "Brijo10"),
                                  model = 2)
summary(cupuacupequi_PDI_quad)

cupuacupequi_PDI_fullcub <- MixModel(frame = cupuacupequi, 
                                     response = "PDI",
                                     mixcomps = c("Cupuacu", "Pequi", "Brijo10"),
                                     model = 3)
summary(cupuacupequi_PDI_fullcub)

cupuacupequi_PDI_splcub <- MixModel(frame = cupuacupequi, 
                                    response = "PDI",
                                    mixcomps = c("Cupuacu", "Pequi", "Brijo10"),
                                    model = 4)
summary(cupuacupequi_PDI_splcub)

AIC(cupuacupequi_PDI_lin, cupuacupequi_PDI_quad, cupuacupequi_PDI_splcub,
    cupuacupequi_PDI_fullcub)

anova(cupuacupequi_PDI_lin, cupuacupequi_PDI_fullcub)
anova(cupuacupequi_PDI_lin, cupuacupequi_PDI_quad)

bgtest(cupuacupequi_PDI_quad)

eplot(predict(cupuacupequi_PDI_quad), rstudent(cupuacupequi_PDI_quad),
     xlab=expression( paste("predicted response ", hat(y)) ) ,ylim=c(-5,5)   ,
     ylab="studentized residuals",type="n",
     main=expression(paste("(",frac(paste("y - ", hat(y)),sigma), ") ~ ", hat(y)   )) ) 

text(predict(cupuacupequi_PDI_splcub), rstudent(cupuacupequi_PDI_splcub),label=1:length(cupuacupequi$PDI)) 
abline(h=c(-2,0,2),lty=c(2,1,2)) 
grid() 


ModelPlot(model = cupuacupequi_PDI_quad,
          res = 500,
          dimensions = list(x1 = "Cupuacu", x2 = "Pequi", x3 = "Brijo10"),
          lims = c(0.166, 0.5, 0.166, 0.5, 0.166, 0.5835),
          at = c(-10, 0.1, 0.3, 0.4, 0.5, 0.6, 0.8, 0.9, 10),
          contour = TRUE, fill = TRUE, pseudo = T, axislab.offset = 0.1,
          color.palette = heat.colors,
          axislabs = c("Cupuaçu", "Pequi", "Brij O10"),
          cornerlabs = c("", "", ""),
          grid = TRUE, grid.pars = list(lty = 0),
          main = "PDI Contour Plot")

#size

cupuacupequi$Size <- as.numeric(cupuacupequi$Size)


cupuacupequi_size_lin <- MixModel(frame = cupuacupequi, 
                                 response = "Size",
                                 mixcomps = c("Cupuacu", "Pequi", "Brijo10"),
                                 model = 1)
summary(cupuacupequi_size_lin)

cupuacupequi_size_quad <- MixModel(frame = cupuacupequi, 
                                  response = "Size",
                                  mixcomps = c("Cupuacu", "Pequi", "Brijo10"),
                                  model = 2)
summary(cupuacupequi_size_quad)

cupuacupequi_size_fullcub <- MixModel(frame = cupuacupequi, 
                                     response = "Size",
                                     mixcomps = c("Cupuacu", "Pequi", "Brijo10"),
                                     model = 3)
summary(cupuacupequi_size_fullcub)

cupuacupequi_size_splcub <- MixModel(frame = cupuacupequi, 
                                    response = "Size",
                                    mixcomps = c("Cupuacu", "Pequi", "Brijo10"),
                                    model = 4)
summary(cupuacupequi_size_splcub)

AIC(cupuacupequi_size_lin, cupuacupequi_size_quad, cupuacupequi_size_splcub,
    cupuacupequi_size_fullcub)

anova(cupuacupequi_size_lin, cupuacupequi_size_fullcub)
anova(cupuacupequi_size_lin, cupuacupequi_size_quad)
anova(cupuacupequi_size_quad, cupuacupequi_size_splcub)

plot(predict(cupuacupequi_size_quad), rstudent(cupuacupequi_size_quad),
     xlab=expression( paste("predicted response ", hat(y)) ) ,ylim=c(-5,5)   ,
     ylab="studentized residuals",type="n",
     main=expression(paste("(",frac(paste("y - ", hat(y)),sigma), ") ~ ", hat(y)   )) ) 

text(predict(cupuacupequi_size_quad), rstudent(cupuacupequi_size_quad),label=1:length(cupuacupequi$Size)) 
abline(h=c(-2,0,2),lty=c(2,1,2)) 
grid() 


ModelPlot(model = cupuacupequi_size_quad,
          res = 500,
          dimensions = list(x1 = "Cupuacu", x2 = "Pequi", x3 = "Brijo10"),
          lims = c(0.166, 0.5, 0.166, 0.5, 0.166, 0.5835),
          at = c(-100, 25, 50, 100, 200, 300, 500, 700, 900, 10000),
          contour = TRUE, fill = TRUE, pseudo = T, axislab.offset = 0.1,
          color.palette = heat.colors,
          axislabs = c("Cupuaçu", "Pequi", "Brij O10"),
          cornerlabs = c("", "", ""),
          grid = TRUE, grid.pars = list(lty = 0),
          main = "Size Contour Plot")

View(cupuacupequi)
summary(cupuacupequi)

#Cupuaçu-Tucumã

cupuacutucuma <- read.csv("NanoCupuacuTucuma.csv")

#PDI

cupuacutucuma_PDI_lin <- MixModel(frame = cupuacutucuma, 
                                 response = "PDI",
                                 mixcomps = c("Cupuacu", "Tucuma", "Brijo10"),
                                 model = 1)
summary(cupuacutucuma_PDI_lin)

cupuacutucuma_PDI_quad <- MixModel(frame = cupuacutucuma, 
                                  response = "PDI",
                                  mixcomps = c("Cupuacu", "Tucuma", "Brijo10"),
                                  model = 2)
summary(cupuacutucuma_PDI_quad)

cupuacutucuma_PDI_fullcub <- MixModel(frame = cupuacutucuma, 
                                     response = "PDI",
                                     mixcomps = c("Cupuacu", "Tucuma", "Brijo10"),
                                     model = 3)
summary(cupuacutucuma_PDI_fullcub)

cupuacutucuma_PDI_splcub <- MixModel(frame = cupuacutucuma, 
                                    response = "PDI",
                                    mixcomps = c("Cupuacu", "Tucuma", "Brijo10"),
                                    model = 4)
summary(cupuacutucuma_PDI_splcub)

cupuacutucuma

AIC(cupuacutucuma_PDI_lin, cupuacutucuma_PDI_quad, cupuacutucuma_PDI_splcub)

bgtest(cupuacutucuma_PDI_lin)

plot(predict(cupuacutucuma_PDI_quad), rstudent(cupuacutucuma_PDI_quad),
     xlab=expression( paste("predicted response ", hat(y)) ) ,ylim=c(-5,5)   ,
     ylab="studentized residuals",type="n",
     main=expression(paste("(",frac(paste("y - ", hat(y)),sigma), ") ~ ", hat(y)   )) ) 

text(predict(cupuacutucuma_PDI_quad), rstudent(cupuatucuma_PDI_splcub),label=1:length(cupuatucuma$PDI)) 
abline(h=c(-2,0,2),lty=c(2,1,2)) 
grid() 


ModelPlot(model = cupuacutucuma_PDI_lin,
          res = 500,
          dimensions = list(x1 = "Cupuacu", x2 = "Tucuma", x3 = "Brijo10"),
          lims = c(0.166, 0.5, 0.166, 0.5, 0.166, 0.5835),
          at = c(-10, 0.2, 0.25, 0.3, 0.35, 10),
          contour = TRUE, fill = TRUE, pseudo = T, axislab.offset = 0.1,
          color.palette = heat.colors,
          axislabs = c("Cupuaçu", "Tucumã", "Brij O10"),
          cornerlabs = c("", "", ""),
          grid = TRUE, grid.pars = list(lty = 0),
          main = "PDI Contour Plot")

cupuacutucuma

#size

cupuacutucuma_size_lin <- MixModel(frame = cupuacutucuma, 
                                  response = "Size",
                                  mixcomps = c("Cupuacu", "Tucuma", "Brijo10"),
                                  model = 1)
summary(cupuacutucuma_size_lin)

cupuacutucuma_size_quad <- MixModel(frame = cupuacutucuma, 
                                   response = "Size",
                                   mixcomps = c("Cupuacu", "Tucuma", "Brijo10"),
                                   model = 2)
summary(cupuacutucuma_size_quad)

cupuacutucuma_size_fullcub <- MixModel(frame = cupuacutucuma, 
                                      response = "Size",
                                      mixcomps = c("Cupuacu", "Tucuma", "Brijo10"),
                                      model = 3)
summary(cupuacutucuma_size_fullcub)

cupuacutucuma_size_splcub <- MixModel(frame = cupuacutucuma, 
                                     response = "Size",
                                     mixcomps = c("Cupuacu", "Tucuma", "Brijo10"),
                                     model = 4)
summary(cupuacutucuma_size_splcub)

AIC(cupuacutucuma_size_lin, cupuacutucuma_size_quad, cupuacutucuma_size_splcub)

anova(cupuacutucuma_size_lin, cupuacutucuma_size_splcub)
anova(cupuacutucuma_size_lin, cupuacutucuma_size_quad)

bgtest(cupuacutucuma_size_lin)

plot(predict(cupuacutucuma_size_splcub), rstudent(cupuacutucuma_size_splcub),
     xlab=expression( paste("predicted response ", hat(y)) ) ,ylim=c(-5,5)   ,
     ylab="studentized residuals",type="n",
     main=expression(paste("(",frac(paste("y - ", hat(y)),sigma), ") ~ ", hat(y)   )) ) 

text(predict(cupuacutucuma_size_splcub), rstudent(cupuacutucuma_size_splcub),label=1:length(cupuacutucuma$size)) 
abline(h=c(-2,0,2),lty=c(2,1,2)) 
grid() 


ModelPlot(model = cupuacutucuma_size_lin,
          res = 500,
          dimensions = list(x1 = "Cupuacu", x2 = "Tucuma", x3 = "Brijo10"),
          lims = c(0.166, 0.5, 0.166, 0.5, 0.166, 0.5835),
          at = c(-1000, 100, 125, 150, 175, 200, 250, 10000),
          contour = TRUE, fill = TRUE, pseudo = T, axislab.offset = 0.1,
          color.palette = heat.colors,
          axislabs = c("Cupuaçu", "Tucumã", "Brij O10"),
          cornerlabs = c("", "", ""),
          grid = TRUE, grid.pars = list(lty = 0),
          main = "Size Contour Plot")

#Manga-Andiroba

mangaandiroba <- read.csv("NanoMangaAndiroba.csv")

View(mangaandiroba)

#PDI

mangaandiroba_PDI_lin <- MixModel(frame = mangaandiroba, 
                                  response = "PDI",
                                  mixcomps = c("Manga", "Andiroba", "Brijo10"),
                                  model = 1)
summary(mangaandiroba_PDI_lin)

mangaandiroba_PDI_quad <- MixModel(frame = mangaandiroba, 
                                   response = "PDI",
                                   mixcomps = c("Manga", "Andiroba", "Brijo10"),
                                   model = 2)
summary(mangaandiroba_PDI_quad)

mangaandiroba_PDI_fullcub <- MixModel(frame = mangaandiroba, 
                                      response = "PDI",
                                      mixcomps = c("Manga", "Andiroba", "Brijo10"),
                                      model = 3)
summary(mangaandiroba_PDI_fullcub)
AIC(mangaandiroba_PDI_fullcub)

mangaandiroba_PDI_fullcub$terms

mangaandiroba_PDI_splcub <- MixModel(frame = mangaandiroba, 
                                     response = "PDI",
                                     mixcomps = c("Manga", "Andiroba", "Brijo10"),
                                     model = 4)
summary(mangaandiroba_PDI_splcub)

AIC(mangaandiroba_PDI_lin, mangaandiroba_PDI_quad,
    mangaandiroba_PDI_splcub, mangaandiroba_PDI_fullcub)

anova(mangaandiroba_PDI_lin, mangaandiroba_PDI_quad)
anova(mangaandiroba_PDI_lin, mangaandiroba_PDI_splcub)
anova(mangaandiroba_PDI_lin, mangaandiroba_PDI_fullcub)


anova(mangaandiroba_PDI_lin,mangaandiroba_PDI_fullcub)

bgtest(mangaandiroba_PDI_lin)

plot(mangaandiroba_PDI_lin)

coeftest(mangaandiroba_PDI_lin, vcov. = vcovHC(mangaandiroba_PDI_lin,
                                              type = "HC1"))

plot(predict(mangaandiroba_PDI_fullcub), rstudent(mangaandiroba_PDI_fullcub),
     xlab=expression( paste("predicted response ", hat(y)) ) ,ylim=c(-5,5)   ,
     ylab="studentized residuals",type="n",
     main=expression(paste("(",frac(paste("y - ", hat(y)),sigma), ") ~ ", hat(y)   )) ) 

text(predict(mangaandiroba_PDI_fullcub), rstudent(mangaandiroba_PDI_fullcub),label=1:length(mangaandiroba$PDI)) 
abline(h=c(-2,0,2),lty=c(2,1,2)) 
grid() 


ModelPlot(model = mangaandiroba_PDI_lin,
          res = 500,
          dimensions = list(x1 = "Manga", x2 = "Andiroba", x3 = "Brijo10"),
          lims = c(0.166, 0.5, 0.166, 0.5, 0.166, 0.5835),
          at = c(-10, 0.4, 0.45, 0.5, 0.55, 0.6, 10),
          contour = TRUE, fill = TRUE, pseudo = T, axislab.offset = 0.1,
          color.palette = heat.colors,
          axislabs = c("Manga", "Andiroba", "Brij O10"),
          cornerlabs = c("", "", ""),
          grid = TRUE, grid.pars = list(lty = 0),
          main = "PDI Contour Plot")

#size

mangaandiroba_size_lin <- MixModel(frame = mangaandiroba, 
                                  response = "Size",
                                  mixcomps = c("Manga", "Andiroba", "Brijo10"),
                                  model = 1)
summary(mangaandiroba_PDI_lin)

mangaandiroba_size_quad <- MixModel(frame = mangaandiroba, 
                                   response = "Size",
                                   mixcomps = c("Manga", "Andiroba", "Brijo10"),
                                   model = 2)
summary(mangaandiroba_PDI_quad)

mangaandiroba_size_fullcub <- MixModel(frame = mangaandiroba, 
                                      response = "Size",
                                      mixcomps = c("Manga", "Andiroba", "Brijo10"),
                                      model = 3)
summary(mangaandiroba_PDI_fullcub)
AIC(mangaandiroba_PDI_fullcub)

mangaandiroba_PDI_fullcub$terms

mangaandiroba_size_splcub <- MixModel(frame = mangaandiroba, 
                                     response = "Size",
                                     mixcomps = c("Manga", "Andiroba", "Brijo10"),
                                     model = 4)
summary(mangaandiroba_PDI_splcub)

AIC(mangaandiroba_size_lin, mangaandiroba_size_quad,
    mangaandiroba_size_splcub, mangaandiroba_size_fullcub)

anova(mangaandiroba_PDI_lin, mangaandiroba_PDI_quad)
anova(mangaandiroba_PDI_lin, mangaandiroba_PDI_splcub)
anova(mangaandiroba_PDI_lin, mangaandiroba_PDI_fullcub)


anova(mangaandiroba_PDI_lin,mangaandiroba_PDI_fullcub)

bgtest(mangaandiroba_size_lin)

plot(mangaandiroba_PDI_lin)

coeftest(mangaandiroba_size_lin, vcov. = vcovHC(mangaandiroba_size_lin,
                                               type = "HC1"))

plot(predict(mangaandiroba_PDI_fullcub), rstudent(mangaandiroba_PDI_fullcub),
     xlab=expression( paste("predicted response ", hat(y)) ) ,ylim=c(-5,5)   ,
     ylab="studentized residuals",type="n",
     main=expression(paste("(",frac(paste("y - ", hat(y)),sigma), ") ~ ", hat(y)   )) ) 

text(predict(mangaandiroba_PDI_fullcub), rstudent(mangaandiroba_PDI_fullcub),label=1:length(mangaandiroba$PDI)) 
abline(h=c(-2,0,2),lty=c(2,1,2)) 
grid() 

mangaandiroba

ModelPlot(model = mangaandiroba_size_lin,
          res = 500,
          dimensions = list(x1 = "Manga", x2 = "Andiroba", x3 = "Brijo10"),
          lims = c(0.166, 0.5, 0.166, 0.5, 0.166, 0.5835),
          #at = c(-10, 0.1, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 10),
          contour = TRUE, fill = TRUE, pseudo = T, axislab.offset = 0.1,
          color.palette = heat.colors,
          axislabs = c("Manga", "Andiroba", "Brij O10"),
          cornerlabs = c("", "", ""),
          grid = TRUE, grid.pars = list(lty = 0),
          main = "Size Contour Plot")

#Manga-Pracaxi

mangapracaxi <- read.csv("NanoMangaPracaxi.csv")

View(mangaandiroba)

#PDI

mangapracaxi_PDI_lin <- MixModel(frame = mangapracaxi, 
                                  response = "PDI",
                                  mixcomps = c("Manga", "Pracaxi", "Brijo10"),
                                  model = 1)
summary(mangaandiroba_PDI_lin)

mangapracaxi_PDI_quad <- MixModel(frame = mangapracaxi, 
                                   response = "PDI",
                                   mixcomps = c("Manga", "Pracaxi", "Brijo10"),
                                   model = 2)
summary(mangaandiroba_PDI_quad)

mangapracaxi_PDI_fullcub <- MixModel(frame = mangapracaxi, 
                                      response = "PDI",
                                      mixcomps = c("Manga", "Pracaxi", "Brijo10"),
                                      model = 3)
summary(mangaandiroba_PDI_fullcub)
AIC(mangaandiroba_PDI_fullcub)

mangaandiroba_PDI_fullcub$terms

mangapracaxi_PDI_splcub <- MixModel(frame = mangapracaxi, 
                                     response = "PDI",
                                     mixcomps = c("Manga", "Pracaxi", "Brijo10"),
                                     model = 4)
summary(mangaandiroba_PDI_splcub)

AIC(mangapracaxi_PDI_lin, mangapracaxi_PDI_quad,
    mangapracaxi_PDI_splcub, mangapracaxi_PDI_fullcub)

anova(mangaandiroba_PDI_lin, mangaandiroba_PDI_quad)
anova(mangaandiroba_PDI_lin, mangaandiroba_PDI_splcub)
anova(mangapracaxi_PDI_lin, mangapracaxi_PDI_fullcub)


anova(mangaandiroba_PDI_lin,mangaandiroba_PDI_fullcub)

bgtest(mangapracaxi_PDI_fullcub)
bgtest(mangapracaxi_PDI_lin)

plot(mangaandiroba_PDI_lin)

coeftest(mangapracaxi_PDI_fullcub, vcov. = vcovHC(mangapracaxi_PDI_fullcub,
                                               type = "HC1"))

plot(predict(mangaandiroba_PDI_fullcub), rstudent(mangaandiroba_PDI_fullcub),
     xlab=expression( paste("predicted response ", hat(y)) ) ,ylim=c(-5,5)   ,
     ylab="studentized residuals",type="n",
     main=expression(paste("(",frac(paste("y - ", hat(y)),sigma), ") ~ ", hat(y)   )) ) 

text(predict(mangaandiroba_PDI_fullcub), rstudent(mangaandiroba_PDI_fullcub),label=1:length(mangaandiroba$PDI)) 
abline(h=c(-2,0,2),lty=c(2,1,2)) 
grid() 

mangapracaxi

ModelPlot(model = mangapracaxi_PDI_lin,
          res = 500,
          dimensions = list(x1 = "Manga", x2 = "Pracaxi", x3 = "Brijo10"),
          lims = c(0.166, 0.5, 0.166, 0.5, 0.166, 0.5835),
          at = c(-10, 0.4, 0.45, 0.5, 0.55, 0.6, 10),
          contour = TRUE, fill = TRUE, pseudo = T, axislab.offset = 0.1,
          color.palette = heat.colors,
          axislabs = c("Manga", "Pracaxi", "Brij O10"),
          cornerlabs = c("", "", ""),
          grid = TRUE, grid.pars = list(lty = 0),
          main = "PDI Contour Plot")

#size

mangapracaxi$Size <- as.numeric()

mangapracaxi_size_lin <- MixModel(frame = mangapracaxi, 
                                   response = "Size",
                                   mixcomps = c("Manga", "Pracaxi", "Brijo10"),
                                   model = 1)
summary(mangaandiroba_PDI_lin)

mangapracaxi_size_quad <- MixModel(frame = mangapracaxi, 
                                    response = "Size",
                                    mixcomps = c("Manga", "Pracaxi", "Brijo10"),
                                    model = 2)
summary(mangaandiroba_PDI_quad)

mangapracaxi_size_fullcub <- MixModel(frame = mangapracaxi, 
                                       response = "Size",
                                       mixcomps = c("Manga", "Pracaxi", "Brijo10"),
                                       model = 3)
summary(mangaandiroba_PDI_fullcub)
AIC(mangaandiroba_PDI_fullcub)

mangaandiroba_PDI_fullcub$terms

mangapracaxi_size_splcub <- MixModel(frame = mangapracaxi, 
                                      response = "Size",
                                      mixcomps = c("Manga", "Pracaxi", "Brijo10"),
                                      model = 4)
summary(mangaandiroba_PDI_splcub)

AIC(mangapracaxi_size_lin, mangapracaxi_size_quad,
    mangapracaxi_size_splcub, mangapracaxi_size_fullcub)

anova(mangaandiroba_PDI_lin, mangaandiroba_PDI_quad)
anova(mangaandiroba_PDI_lin, mangaandiroba_PDI_splcub)
anova(mangapracaxi_PDI_lin, mangapracaxi_PDI_fullcub)


anova(mangaandiroba_PDI_lin,mangaandiroba_PDI_fullcub)

bgtest(mangapracaxi_size_lin)

plot(mangaandiroba_PDI_lin)

coeftest(mangaandiroba_size_lin, vcov. = vcovHC(mangaandiroba_size_lin,
                                                type = "HC1"))

plot(predict(mangaandiroba_PDI_fullcub), rstudent(mangaandiroba_PDI_fullcub),
     xlab=expression( paste("predicted response ", hat(y)) ) ,ylim=c(-5,5)   ,
     ylab="studentized residuals",type="n",
     main=expression(paste("(",frac(paste("y - ", hat(y)),sigma), ") ~ ", hat(y)   )) ) 

text(predict(mangaandiroba_PDI_fullcub), rstudent(mangaandiroba_PDI_fullcub),label=1:length(mangaandiroba$PDI)) 
abline(h=c(-2,0,2),lty=c(2,1,2)) 
grid() 

mangaandiroba

ModelPlot(model = mangapracaxi_size_lin,
          res = 500,
          dimensions = list(x1 = "Manga", x2 = "Pracaxi", x3 = "Brijo10"),
          lims = c(0.166, 0.5, 0.166, 0.5, 0.166, 0.5835),
          at = c(-10, 500, 800, 1100, 1400, 1700, 2100, 3000),
          contour = TRUE, fill = TRUE, pseudo = T, axislab.offset = 0.1,
          color.palette = heat.colors,
          axislabs = c("Manga", "Pracaxi", "Brij O10"),
          cornerlabs = c("", "", ""),
          grid = TRUE, grid.pars = list(lty = 0),
          main = "Size Contour Plot")

#Mururu-Açaí

murumuruacai <- read.csv("NanoMurumuruAcai.csv")

View(mangaandiroba)

#PDI

murumuruacai_PDI_lin <- MixModel(frame = murumuruacai, 
                                 response = "PDI",
                                 mixcomps = c("Murumuru", "Acai", "Brijo10"),
                                 model = 1)
summary(mangaandiroba_PDI_lin)

murumuruacai_PDI_quad <- MixModel(frame = murumuruacai, 
                                  response = "PDI",
                                  mixcomps = c("Murumuru", "Acai", "Brijo10"),
                                  model = 2)
summary(mangaandiroba_PDI_quad)

murumuruacai_PDI_fullcub <- MixModel(frame = murumuruacai, 
                                     response = "PDI",
                                     mixcomps = c("Murumuru", "Acai", "Brijo10"),
                                     model = 3)
summary(mangaandiroba_PDI_fullcub)
AIC(mangaandiroba_PDI_fullcub)

mangaandiroba_PDI_fullcub$terms

murumuruacai_PDI_splcub <- MixModel(frame = murumuruacai, 
                                    response = "PDI",
                                    mixcomps = c("Murumuru", "Acai", "Brijo10"),
                                    model = 4)
summary(mangaandiroba_PDI_splcub)

AIC(murumuruacai_PDI_lin, murumuruacai_PDI_quad,
    murumuruacai_PDI_splcub, murumuruacai_PDI_fullcub)

anova(mangaandiroba_PDI_lin, mangaandiroba_PDI_quad)
anova(mangaandiroba_PDI_lin, mangaandiroba_PDI_splcub)
anova(murumuruacai_PDI_lin, murumuruacai_PDI_fullcub)


anova(mangaandiroba_PDI_lin,mangaandiroba_PDI_fullcub)

bgtest(mangapracaxi_PDI_fullcub)
bgtest(murumuruacai_PDI_lin)

plot(mangaandiroba_PDI_lin)

coeftest(mangapracaxi_PDI_fullcub, vcov. = vcovHC(mangapracaxi_PDI_fullcub,
                                                  type = "HC1"))

plot(predict(mangaandiroba_PDI_fullcub), rstudent(mangaandiroba_PDI_fullcub),
     xlab=expression( paste("predicted response ", hat(y)) ) ,ylim=c(-5,5)   ,
     ylab="studentized residuals",type="n",
     main=expression(paste("(",frac(paste("y - ", hat(y)),sigma), ") ~ ", hat(y)   )) ) 

text(predict(mangaandiroba_PDI_fullcub), rstudent(mangaandiroba_PDI_fullcub),label=1:length(mangaandiroba$PDI)) 
abline(h=c(-2,0,2),lty=c(2,1,2)) 
grid() 

mangapracaxi

ModelPlot(model = murumuruacai_PDI_lin,
          res = 500,
          dimensions = list(x1 = "Murumuru", x2 = "Acai", x3 = "Brijo10"),
          lims = c(0.166, 0.5, 0.166, 0.5, 0.166, 0.5835),
          at = c(-10, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 10),
          contour = TRUE, fill = TRUE, pseudo = T, axislab.offset = 0.1,
          color.palette = heat.colors,
          axislabs = c("Murumuru", "Açaí", "Brijo10"),
          cornerlabs = c("", "", ""),
          grid = TRUE, grid.pars = list(lty = 0),
          main = "PDI Contour Plot")

#size

mangapracaxi$Size <- as.numeric()

murumuruacai_size_lin <- MixModel(frame = murumuruacai, 
                                  response = "Size",
                                  mixcomps = c("Murumuru", "Acai", "Brijo10"),
                                  model = 1)
summary(mangaandiroba_PDI_lin)

murumuruacai_size_quad <- MixModel(frame = murumuruacai, 
                                   response = "Size",
                                   mixcomps = c("Murumuru", "Acai", "Brijo10"),
                                   model = 2)
summary(mangaandiroba_PDI_quad)

murumuruacai_size_fullcub <- MixModel(frame = murumuruacai, 
                                      response = "Size",
                                      mixcomps = c("Murumuru", "Acai", "Brijo10"),
                                      model = 3)
summary(mangaandiroba_PDI_fullcub)
AIC(mangaandiroba_PDI_fullcub)

mangaandiroba_PDI_fullcub$terms

murumuruacai_size_splcub <- MixModel(frame = murumuruacai, 
                                     response = "Size",
                                     mixcomps = c("Murumuru", "Acai", "Brijo10"),
                                     model = 4)
summary(mangaandiroba_PDI_splcub)

AIC(murumuruacai_size_lin, murumuruacai_size_quad,
    murumuruacai_size_splcub, murumuruacai_size_fullcub)

anova(mangaandiroba_PDI_lin, mangaandiroba_PDI_quad)
anova(mangaandiroba_PDI_lin, mangaandiroba_PDI_splcub)
anova(murumuruacai_PDI_lin, murumuruacai_PDI_fullcub)


anova(mangaandiroba_PDI_lin,mangaandiroba_PDI_fullcub)

bgtest(murumuruacai_size_lin)

plot(mangaandiroba_PDI_lin)

coeftest(mangaandiroba_size_lin, vcov. = vcovHC(mangaandiroba_size_lin,
                                                type = "HC1"))

plot(predict(mangaandiroba_PDI_fullcub), rstudent(mangaandiroba_PDI_fullcub),
     xlab=expression( paste("predicted response ", hat(y)) ) ,ylim=c(-5,5)   ,
     ylab="studentized residuals",type="n",
     main=expression(paste("(",frac(paste("y - ", hat(y)),sigma), ") ~ ", hat(y)   )) ) 

text(predict(mangaandiroba_PDI_fullcub), rstudent(mangaandiroba_PDI_fullcub),label=1:length(mangaandiroba$PDI)) 
abline(h=c(-2,0,2),lty=c(2,1,2)) 
grid() 

murumuruacai

ModelPlot(model = murumuruacai_size_lin,
          res = 500,
          dimensions = list(x1 = "Murumuru", x2 = "Acai", x3 = "Brijo10"),
          lims = c(0.166, 0.5, 0.166, 0.5, 0.166, 0.5835),
          at = c(-1000, 10, 100, 200, 300, 400, 1000),
          contour = TRUE, fill = TRUE, pseudo = T, axislab.offset = 0.1,
          color.palette = heat.colors,
          axislabs = c("Murumuru", "Açaí", "Brijo10"),
          cornerlabs = c("", "", ""),
          grid = TRUE, grid.pars = list(lty = 0),
          main = "Size Contour Plot")

#Mururu-Buriti

murumuruburiti <- read.csv("NanoMurumuruBuriti.csv")

View(mangaandiroba)

#PDI

murumuruburiti_PDI_lin <- MixModel(frame = murumuruburiti, 
                                 response = "PDI",
                                 mixcomps = c("Murumuru", "Buriti", "Brijo10"),
                                 model = 1)
summary(mangaandiroba_PDI_lin)

murumuruburiti_PDI_quad <- MixModel(frame = murumuruburiti, 
                                  response = "PDI",
                                  mixcomps = c("Murumuru", "Buriti", "Brijo10"),
                                  model = 2)
summary(mangaandiroba_PDI_quad)

murumuruburiti_PDI_fullcub <- MixModel(frame = murumuruburiti, 
                                     response = "PDI",
                                     mixcomps = c("Murumuru", "Buriti", "Brijo10"),
                                     model = 3)
summary(mangaandiroba_PDI_fullcub)
AIC(mangaandiroba_PDI_fullcub)

mangaandiroba_PDI_fullcub$terms

murumuruburiti_PDI_splcub <- MixModel(frame = murumuruburiti, 
                                    response = "PDI",
                                    mixcomps = c("Murumuru", "Buriti", "Brijo10"),
                                    model = 4)
summary(mangaandiroba_PDI_splcub)

AIC(murumuruburiti_PDI_lin, murumuruburiti_PDI_quad,
    murumuruburiti_PDI_splcub, murumuruburiti_PDI_fullcub)

anova(murumuruburiti_PDI_lin, murumuruburiti_PDI_quad)
anova(mangaandiroba_PDI_lin, mangaandiroba_PDI_splcub)
anova(murumuruburiti_PDI_lin, murumuruburiti_PDI_fullcub)


anova(mangaandiroba_PDI_lin,mangaandiroba_PDI_fullcub)

bgtest(mangapracaxi_PDI_fullcub)
bgtest(murumuruburiti_PDI_lin)

plot(mangaandiroba_PDI_lin)

coeftest(mangapracaxi_PDI_fullcub, vcov. = vcovHC(mangapracaxi_PDI_fullcub,
                                                  type = "HC1"))

plot(predict(mangaandiroba_PDI_fullcub), rstudent(mangaandiroba_PDI_fullcub),
     xlab=expression( paste("predicted response ", hat(y)) ) ,ylim=c(-5,5)   ,
     ylab="studentized residuals",type="n",
     main=expression(paste("(",frac(paste("y - ", hat(y)),sigma), ") ~ ", hat(y)   )) ) 

text(predict(mangaandiroba_PDI_fullcub), rstudent(mangaandiroba_PDI_fullcub),label=1:length(mangaandiroba$PDI)) 
abline(h=c(-2,0,2),lty=c(2,1,2)) 
grid() 

murumuruburiti

ModelPlot(model = murumuruburiti_PDI_lin,
          res = 500,
          dimensions = list(x1 = "Murumuru", x2 = "Buriti", x3 = "Brijo10"),
          lims = c(0.166, 0.5, 0.166, 0.5, 0.166, 0.5835),
          #at = c(-10, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 10),
          contour = TRUE, fill = TRUE, pseudo = T, axislab.offset = 0.1,
          color.palette = heat.colors,
          axislabs = c("Murumuru", "Buriti", "Brijo10"),
          cornerlabs = c("", "", ""),
          grid = TRUE, grid.pars = list(lty = 0),
          main = "PDI Contour Plot")

#size

mangapracaxi$Size <- as.numeric()

murumuruburiti_size_lin <- MixModel(frame = murumuruburiti, 
                                  response = "Size",
                                  mixcomps = c("Murumuru", "Buriti", "Brijo10"),
                                  model = 1)
summary(mangaandiroba_PDI_lin)

murumuruburiti_size_quad <- MixModel(frame = murumuruburiti, 
                                   response = "Size",
                                   mixcomps = c("Murumuru", "Buriti", "Brijo10"),
                                   model = 2)
summary(mangaandiroba_PDI_quad)

murumuruburiti_size_fullcub <- MixModel(frame = murumuruburiti, 
                                      response = "Size",
                                      mixcomps = c("Murumuru", "Buriti", "Brijo10"),
                                      model = 3)
summary(mangaandiroba_PDI_fullcub)
AIC(mangaandiroba_PDI_fullcub)

mangaandiroba_PDI_fullcub$terms

murumuruburiti_size_splcub <- MixModel(frame = murumuruburiti, 
                                     response = "Size",
                                     mixcomps = c("Murumuru", "Buriti", "Brijo10"),
                                     model = 4)
summary(mangaandiroba_PDI_splcub)

AIC(murumuruburiti_size_lin, murumuruburiti_size_quad,
    murumuruburiti_size_splcub, murumuruburiti_size_fullcub)

anova(murumuruburiti_size_lin, murumuruburiti_size_quad)
anova(mangaandiroba_PDI_lin, mangaandiroba_PDI_splcub)
anova(murumuruburiti_size_lin, murumuruburiti_size_fullcub)


anova(mangaandiroba_PDI_lin,mangaandiroba_PDI_fullcub)

bgtest(murumuruburiti_size_lin)
bgtest(murumuruburiti_size_quad)

plot(mangaandiroba_PDI_lin)

coeftest(mangaandiroba_size_lin, vcov. = vcovHC(mangaandiroba_size_lin,
                                                type = "HC1"))

plot(predict(mangaandiroba_PDI_fullcub), rstudent(mangaandiroba_PDI_fullcub),
     xlab=expression( paste("predicted response ", hat(y)) ) ,ylim=c(-5,5)   ,
     ylab="studentized residuals",type="n",
     main=expression(paste("(",frac(paste("y - ", hat(y)),sigma), ") ~ ", hat(y)   )) ) 

text(predict(mangaandiroba_PDI_fullcub), rstudent(mangaandiroba_PDI_fullcub),label=1:length(mangaandiroba$PDI)) 
abline(h=c(-2,0,2),lty=c(2,1,2)) 
grid() 

murumuruacai

ModelPlot(model = murumuruburiti_size_quad,
          res = 500,
          dimensions = list(x1 = "Murumuru", x2 = "Buriti", x3 = "Brijo10"),
          lims = c(0.166, 0.5, 0.166, 0.5, 0.166, 0.5835),
          at = c(-1000, 10, 50, 100, 150, 200, 300, 400, 1000),
          contour = TRUE, fill = TRUE, pseudo = T, axislab.offset = 0.1,
          color.palette = heat.colors,
          axislabs = c("Murumuru", "Buriti", "Brijo10"),
          cornerlabs = c("", "", ""),
          grid = TRUE, grid.pars = list(lty = 0),
          main = "Size Contour Plot")

#Mururu-Jambu

murumurujambu <- read.csv("NanoMurumuruJambu.csv")

View(mangaandiroba)

murumurujambu

#PDI

murumurujambu_PDI_lin <- MixModel(frame = murumurujambu, 
                                   response = "PDI",
                                   mixcomps = c("Murumuru", "Jambu", "Brijo10"),
                                   model = 1)
summary(mangaandiroba_PDI_lin)

murumurujambu_PDI_quad <- MixModel(frame = murumurujambu, 
                                    response = "PDI",
                                    mixcomps = c("Murumuru", "Jambu", "Brijo10"),
                                    model = 2)
summary(mangaandiroba_PDI_quad)

murumurujambu_PDI_fullcub <- MixModel(frame = murumurujambu, 
                                       response = "PDI",
                                       mixcomps = c("Murumuru", "Jambu", "Brijo10"),
                                       model = 3)
summary(mangaandiroba_PDI_fullcub)
AIC(mangaandiroba_PDI_fullcub)

mangaandiroba_PDI_fullcub$terms

murumurujambu_PDI_splcub <- MixModel(frame = murumurujambu, 
                                      response = "PDI",
                                      mixcomps = c("Murumuru", "Jambu", "Brijo10"),
                                      model = 4)
summary(mangaandiroba_PDI_splcub)

AIC(murumurujambu_PDI_lin, murumurujambu_PDI_quad,
    murumurujambu_PDI_splcub, murumurujambu_PDI_fullcub)

anova(murumuruburiti_PDI_lin, murumuruburiti_PDI_quad)
anova(murumurujambu_PDI_lin, murumurujambu_PDI_splcub)
anova(murumuruburiti_PDI_lin, murumuruburiti_PDI_fullcub)


anova(mangaandiroba_PDI_lin,mangaandiroba_PDI_fullcub)

bgtest(mangapracaxi_PDI_fullcub)
bgtest(murumurujambu_PDI_lin)

plot(mangaandiroba_PDI_lin)

coeftest(mangapracaxi_PDI_fullcub, vcov. = vcovHC(mangapracaxi_PDI_fullcub,
                                                  type = "HC1"))

plot(predict(murumurujambu_PDI_lin), rstudent(mangaandiroba_PDI_fullcub),
     xlab=expression( paste("predicted response ", hat(y)) ) ,ylim=c(-5,5)   ,
     ylab="studentized residuals",type="n",
     main=expression(paste("(",frac(paste("y - ", hat(y)),sigma), ") ~ ", hat(y)   )) ) 

text(predict(mangaandiroba_PDI_fullcub), rstudent(mangaandiroba_PDI_fullcub),label=1:length(mangaandiroba$PDI)) 
abline(h=c(-2,0,2),lty=c(2,1,2)) 
grid() 

murumuruburiti

ModelPlot(model = murumurujambu_PDI_lin,
          res = 500,
          dimensions = list(x1 = "Murumuru", x2 = "Jambu", x3 = "Brijo10"),
          lims = c(0.166, 0.5, 0.166, 0.5, 0.166, 0.5835),
          #at = c(-10, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 10),
          contour = TRUE, fill = TRUE, pseudo = T, axislab.offset = 0.1,
          color.palette = heat.colors,
          axislabs = c("Murumuru", "Jambu", "Brijo10"),
          cornerlabs = c("", "", ""),
          grid = TRUE, grid.pars = list(lty = 0),
          main = "PDI Contour Plot")

#size

mangapracaxi$Size <- as.numeric()

murumurujambu_size_lin <- MixModel(frame = murumurujambu, 
                                    response = "Size",
                                    mixcomps = c("Murumuru", "Jambu", "Brijo10"),
                                    model = 1)
summary(mangaandiroba_PDI_lin)

murumurujambu_size_quad <- MixModel(frame = murumurujambu, 
                                     response = "Size",
                                     mixcomps = c("Murumuru", "Jambu", "Brijo10"),
                                     model = 2)
summary(mangaandiroba_PDI_quad)

murumurujambu_size_fullcub <- MixModel(frame = murumurujambu, 
                                        response = "Size",
                                        mixcomps = c("Murumuru", "Jambu", "Brijo10"),
                                        model = 3)
summary(mangaandiroba_PDI_fullcub)
AIC(mangaandiroba_PDI_fullcub)

mangaandiroba_PDI_fullcub$terms

murumurujambu_size_splcub <- MixModel(frame = murumurujambu, 
                                       response = "Size",
                                       mixcomps = c("Murumuru", "Jambu", "Brijo10"),
                                       model = 4)
summary(mangaandiroba_PDI_splcub)

AIC(murumurujambu_size_lin, murumurujambu_size_quad,
    murumurujambu_size_splcub, murumurujambu_size_fullcub)

anova(murumurujambu_size_lin, murumurujambu_size_quad)
anova(murumurujambu_PDI_lin, murumurujambu_PDI_splcub)
anova(murumurujambu_size_lin, murumurujambu_size_fullcub)

anova(murumurujambu_size_quad, murumurujambu_size_fullcub)


anova(mangaandiroba_PDI_lin,mangaandiroba_PDI_fullcub)

bgtest(murumuruburiti_size_lin)
bgtest(murumurujambu_size_quad)

plot(mangaandiroba_PDI_lin)

coeftest(murumurujambu_size_quad, vcov. = vcovHC(murumurujambu_size_quad,
                                                type = "HC1"))

plot(predict(mangaandiroba_PDI_fullcub), rstudent(mangaandiroba_PDI_fullcub),
     xlab=expression( paste("predicted response ", hat(y)) ) ,ylim=c(-5,5)   ,
     ylab="studentized residuals",type="n",
     main=expression(paste("(",frac(paste("y - ", hat(y)),sigma), ") ~ ", hat(y)   )) ) 

text(predict(mangaandiroba_PDI_fullcub), rstudent(mangaandiroba_PDI_fullcub),label=1:length(mangaandiroba$PDI)) 
abline(h=c(-2,0,2),lty=c(2,1,2)) 
grid() 

murumuruacai

ModelPlot(model = murumurujambu_size_quad,
          res = 500,
          dimensions = list(x1 = "Murumuru", x2 = "Jambu", x3 = "Brijo10"),
          lims = c(0.166, 0.5, 0.166, 0.5, 0.166, 0.5835),
          at = c(-1000, 100, 500, 1000, 1500, 2000, 10000),
          contour = TRUE, fill = TRUE, pseudo = T, axislab.offset = 0.1,
          color.palette = heat.colors,
          axislabs = c("Murumuru", "Jambu", "Brijo10"),
          cornerlabs = c("", "", ""),
          grid = TRUE, grid.pars = list(lty = 0),
          main = "Size Contour Plot")

#Tucumã-Tucumã

tucumatucuma <- read.csv("NanoTucumaTucuma.csv")

View(mangaandiroba)

tucumatucuma

#PDI

tucumatucuma_PDI_lin <- MixModel(frame = tucumatucuma, 
                                  response = "PDI",
                                  mixcomps = c("Tucuma", "Tucuma.1", "Brijo10"),
                                  model = 1)
summary(mangaandiroba_PDI_lin)

tucumatucuma_PDI_quad <- MixModel(frame = tucumatucuma, 
                                   response = "PDI",
                                   mixcomps = c("Tucuma", "Tucuma.1", "Brijo10"),
                                   model = 2)
summary(mangaandiroba_PDI_quad)

tucumatucuma_PDI_fullcub <- MixModel(frame = tucumatucuma, 
                                      response = "PDI",
                                      mixcomps = c("Tucuma", "Tucuma.1", "Brijo10"),
                                      model = 3)
summary(mangaandiroba_PDI_fullcub)
AIC(mangaandiroba_PDI_fullcub)

mangaandiroba_PDI_fullcub$terms

tucumatucuma_PDI_splcub <- MixModel(frame = tucumatucuma, 
                                     response = "PDI",
                                     mixcomps = c("Tucuma", "Tucuma.1", "Brijo10"),
                                     model = 4)
summary(mangaandiroba_PDI_splcub)

AIC(tucumatucuma_PDI_lin, tucumatucuma_PDI_quad,
    tucumatucuma_PDI_splcub, tucumatucuma_PDI_fullcub)

anova(tucumatucuma_PDI_lin, tucumatucuma_PDI_quad)
anova(murumurujambu_PDI_lin, murumurujambu_PDI_splcub)
anova(tucumatucuma_PDI_lin, tucumatucuma_PDI_fullcub)


anova(tucumatucuma_PDI_lin,tucumatucuma_PDI_quad)

bgtest(tucumatucuma_PDI_quad)
bgtest(murumurujambu_PDI_lin)

plot(mangaandiroba_PDI_lin)

coeftest(mangapracaxi_PDI_fullcub, vcov. = vcovHC(mangapracaxi_PDI_fullcub,
                                                  type = "HC1"))

plot(predict(murumurujambu_PDI_lin), rstudent(mangaandiroba_PDI_fullcub),
     xlab=expression( paste("predicted response ", hat(y)) ) ,ylim=c(-5,5)   ,
     ylab="studentized residuals",type="n",
     main=expression(paste("(",frac(paste("y - ", hat(y)),sigma), ") ~ ", hat(y)   )) ) 

text(predict(mangaandiroba_PDI_fullcub), rstudent(mangaandiroba_PDI_fullcub),label=1:length(mangaandiroba$PDI)) 
abline(h=c(-2,0,2),lty=c(2,1,2)) 
grid() 

murumuruburiti

ModelPlot(model = tucumatucuma_PDI_quad,
          res = 500,
          dimensions = list(x1 = "Tucuma", x2 = "Tucuma.1", x3 = "Brijo10"),
          lims = c(0.166, 0.5, 0.166, 0.5, 0.166, 0.5835),
          at = c(-10, 0.1, 0.2, 0.3, 0.4, 0.5, 10),
          contour = TRUE, fill = TRUE, pseudo = T, axislab.offset = 0.1,
          color.palette = heat.colors,
          axislabs = c("Tucuma", "Tucuma", "Brijo10"),
          cornerlabs = c("", "", ""),
          grid = TRUE, grid.pars = list(lty = 0),
          main = "PDI Contour Plot")

#size

mangapracaxi$Size <- as.numeric()

tucumatucuma_size_lin <- MixModel(frame = tucumatucuma, 
                                   response = "Size",
                                   mixcomps = c("Tucuma", "Tucuma.1", "Brijo10"),
                                   model = 1)
summary(mangaandiroba_PDI_lin)

tucumatucuma_size_quad <- MixModel(frame = tucumatucuma, 
                                    response = "Size",
                                    mixcomps = c("Tucuma", "Tucuma.1", "Brijo10"),
                                    model = 2)
summary(mangaandiroba_PDI_quad)

tucumatucuma_size_fullcub <- MixModel(frame = tucumatucuma, 
                                       response = "Size",
                                       mixcomps = c("Tucuma", "Tucuma.1", "Brijo10"),
                                       model = 3)
summary(mangaandiroba_PDI_fullcub)
AIC(mangaandiroba_PDI_fullcub)

mangaandiroba_PDI_fullcub$terms

tucumatucuma_size_splcub <- MixModel(frame = tucumatucuma, 
                                      response = "Size",
                                      mixcomps = c("Tucuma", "Tucuma.1", "Brijo10"),
                                      model = 4)
summary(mangaandiroba_PDI_splcub)

AIC(tucumatucuma_size_lin, tucumatucuma_size_quad,
    tucumatucuma_size_splcub, tucumatucuma_size_fullcub)

anova(tucumatucuma_size_lin, tucumatucuma_size_quad)
anova(murumurujambu_PDI_lin, murumurujambu_PDI_splcub)
anova(tucumatucuma_size_lin, tucumatucuma_size_fullcub)

anova(tucumatucuma_size_quad, tucumatucuma_size_fullcub)


anova(mangaandiroba_PDI_lin,mangaandiroba_PDI_fullcub)

bgtest(murumuruburiti_size_lin)
bgtest(tucumatucuma_size_quad)

plot(mangaandiroba_PDI_lin)

hatvalues(tucumatucuma_size_quad)
n <- nrow(tucumatucuma)
p <- length(coef(tucumatucuma_size_quad))
cutoff <- 2 * p / n
cutoff
cutoffxt <- 3 * p / n
cutoffxt

coeftest(tucumatucuma_size_quad, vcov. = vcovHC(tucumatucuma_size_quad,
                                                 type = "HC1"))

plot(predict(mangaandiroba_PDI_fullcub), rstudent(mangaandiroba_PDI_fullcub),
     xlab=expression( paste("predicted response ", hat(y)) ) ,ylim=c(-5,5)   ,
     ylab="studentized residuals",type="n",
     main=expression(paste("(",frac(paste("y - ", hat(y)),sigma), ") ~ ", hat(y)   )) ) 

text(predict(mangaandiroba_PDI_fullcub), rstudent(mangaandiroba_PDI_fullcub),label=1:length(mangaandiroba$PDI)) 
abline(h=c(-2,0,2),lty=c(2,1,2)) 
grid() 

murumuruacai

ModelPlot(model = tucumatucuma_size_quad,
          res = 500,
          dimensions = list(x1 = "Tucuma", x2 = "Tucuma.1", x3 = "Brijo10"),
          lims = c(0.166, 0.5, 0.166, 0.5, 0.166, 0.5835),
          at = c(-1000, 25, 50, 100, 150, 200, 250, 10000),
          contour = TRUE, fill = TRUE, pseudo = T, axislab.offset = 0.1,
          color.palette = heat.colors,
          axislabs = c("Tucuma", "Tucuma", "Brijo10"),
          cornerlabs = c("", "", ""),
          grid = TRUE, grid.pars = list(lty = 0),
          main = "Size Contour Plot")

#Tucumã-Urucum

tucumaurucum <- read.csv("NanoTucumaUrucum.csv")

View(mangaandiroba)

tucumaurucum

#PDI

tucumaurucum_PDI_lin <- MixModel(frame = tucumaurucum, 
                                 response = "PDI",
                                 mixcomps = c("Tucuma", "Urucum", "Brijo10"),
                                 model = 1)
summary(mangaandiroba_PDI_lin)

tucumaurucum_PDI_quad <- MixModel(frame = tucumaurucum, 
                                  response = "PDI",
                                  mixcomps = c("Tucuma", "Urucum", "Brijo10"),
                                  model = 2)
summary(mangaandiroba_PDI_quad)

tucumaurucum_PDI_fullcub <- MixModel(frame = tucumaurucum, 
                                     response = "PDI",
                                     mixcomps = c("Tucuma", "Urucum", "Brijo10"),
                                     model = 3)
summary(mangaandiroba_PDI_fullcub)
AIC(mangaandiroba_PDI_fullcub)

mangaandiroba_PDI_fullcub$terms

tucumaurucum_PDI_splcub <- MixModel(frame = tucumaurucum, 
                                    response = "PDI",
                                    mixcomps = c("Tucuma", "Urucum", "Brijo10"),
                                    model = 4)
summary(mangaandiroba_PDI_splcub)

AIC(tucumaurucum_PDI_lin, tucumaurucum_PDI_quad,
    tucumaurucum_PDI_splcub, tucumaurucum_PDI_fullcub)

anova(tucumatucuma_PDI_lin, tucumatucuma_PDI_quad)
anova(murumurujambu_PDI_lin, murumurujambu_PDI_splcub)
anova(tucumatucuma_PDI_lin, tucumatucuma_PDI_fullcub)


anova(tucumatucuma_PDI_lin,tucumatucuma_PDI_quad)

bgtest(tucumatucuma_PDI_quad)
bgtest(tucumaurucum_PDI_lin)

plot(mangaandiroba_PDI_lin)

coeftest(mangapracaxi_PDI_fullcub, vcov. = vcovHC(mangapracaxi_PDI_fullcub,
                                                  type = "HC1"))

plot(predict(murumurujambu_PDI_lin), rstudent(mangaandiroba_PDI_fullcub),
     xlab=expression( paste("predicted response ", hat(y)) ) ,ylim=c(-5,5)   ,
     ylab="studentized residuals",type="n",
     main=expression(paste("(",frac(paste("y - ", hat(y)),sigma), ") ~ ", 
                           hat(y)   )) ) 

text(predict(mangaandiroba_PDI_fullcub), rstudent(mangaandiroba_PDI_fullcub),
     label=1:length(mangaandiroba$PDI)) 
abline(h=c(-2,0,2),lty=c(2,1,2)) 
grid() 

murumuruburiti

ModelPlot(model = tucumaurucum_PDI_lin,
          res = 500,
          dimensions = list(x1 = "Tucuma", x2 = "Urucum", x3 = "Brijo10"),
          lims = c(0.166, 0.5, 0.166, 0.5, 0.166, 0.5835),
          #at = c(-10, 0.1, 0.2, 0.3, 0.4, 0.5, 10),
          contour = TRUE, fill = TRUE, pseudo = T, axislab.offset = 0.1,
          color.palette = heat.colors,
          axislabs = c("Tucumã", "Urucum", "Brijo10"),
          cornerlabs = c("", "", ""),
          grid = TRUE, grid.pars = list(lty = 0),
          main = "PDI Contour Plot")

#size

mangapracaxi$Size <- as.numeric()

tucumaurucum_size_lin <- MixModel(frame = tucumaurucum, 
                                  response = "Size",
                                  mixcomps = c("Tucuma", "Urucum", "Brijo10"),
                                  model = 1)
summary(mangaandiroba_PDI_lin)

tucumaurucum_size_quad <- MixModel(frame = tucumaurucum, 
                                   response = "Size",
                                   mixcomps = c("Tucuma", "Urucum", "Brijo10"),
                                   model = 2)
summary(mangaandiroba_PDI_quad)

tucumaurucum_size_fullcub <- MixModel(frame = tucumaurucum, 
                                      response = "Size",
                                      mixcomps = c("Tucuma", "Urucum", "Brijo10"),
                                      model = 3)
summary(mangaandiroba_PDI_fullcub)
AIC(mangaandiroba_PDI_fullcub)

mangaandiroba_PDI_fullcub$terms

tucumaurucum_size_splcub <- MixModel(frame = tucumaurucum, 
                                     response = "Size",
                                     mixcomps = c("Tucuma", "Urucum", "Brijo10"),
                                     model = 4)
summary(mangaandiroba_PDI_splcub)

AIC(tucumaurucum_size_lin, tucumaurucum_size_quad,
    tucumaurucum_size_splcub, tucumaurucum_size_fullcub)

anova(tucumaurucum_size_lin, tucumaurucum_size_quad)
anova(murumurujambu_PDI_lin, murumurujambu_PDI_splcub)
anova(tucumaurucum_size_lin, tucumaurucum_size_fullcub)

anova(tucumatucuma_size_quad, tucumatucuma_size_fullcub)


anova(mangaandiroba_PDI_lin,mangaandiroba_PDI_fullcub)

bgtest(tucumaurucum_size_lin)
bgtest(tucumatucuma_size_quad)

plot(mangaandiroba_PDI_lin)

hatvalues(tucumatucuma_size_quad)
n <- nrow(tucumatucuma)
p <- length(coef(tucumatucuma_size_quad))
cutoff <- 2 * p / n
cutoff
cutoffxt <- 3 * p / n
cutoffxt

coeftest(tucumatucuma_size_quad, vcov. = vcovHC(tucumatucuma_size_quad,
                                                type = "HC1"))

plot(predict(mangaandiroba_PDI_fullcub), rstudent(mangaandiroba_PDI_fullcub),
     xlab=expression( paste("predicted response ", hat(y)) ) ,ylim=c(-5,5)   ,
     ylab="studentized residuals",type="n",
     main=expression(paste("(",frac(paste("y - ", hat(y)),sigma), ") ~ ", hat(y)   )) ) 

text(predict(mangaandiroba_PDI_fullcub), rstudent(mangaandiroba_PDI_fullcub),label=1:length(mangaandiroba$PDI)) 
abline(h=c(-2,0,2),lty=c(2,1,2)) 
grid() 

murumuruacai

ModelPlot(model = tucumaurucum_size_lin,
          res = 500,
          dimensions = list(x1 = "Tucuma", x2 = "Urucum", x3 = "Brijo10"),
          lims = c(0.166, 0.5, 0.166, 0.5, 0.166, 0.5835),
          at = c(-1000, 200, 400, 600, 800, 1000, 1200, 1400, 1600, 10000),
          contour = TRUE, fill = TRUE, pseudo = T, axislab.offset = 0.1,
          color.palette = heat.colors,
          axislabs = c("Tucumã", "Urucum", "Brijo10"),
          cornerlabs = c("", "", ""),
          grid = TRUE, grid.pars = list(lty = 0),
          main = "Size Contour Plot")

#Ucuuba-Andiroba

ucuubaandiroba <- read.csv("NanoUcuubaAndiroba.csv")

View(mangaandiroba)

ucuubaandiroba

colnames(ucuubaandiroba)[2] <- "Andiroba"

#PDI

ucuubaandiroba_PDI_lin <- MixModel(frame = ucuubaandiroba, 
                                 response = "PDI",
                                 mixcomps = c("Ucuuba", "Andiroba", "Brijo10"),
                                 model = 1)
summary(ucuubaandiroba_PDI_lin)

anova(ucuubaandiroba_PDI_lin)

plot(ucuubaandiroba_PDI_lin)

check_model(ucuubaandiroba_PDI_lin)

plot(simulateResiduals(ucuubaandiroba_PDI_lin))


ucuubaandiroba_PDI_quad <- MixModel(frame = ucuubaandiroba, 
                                  response = "PDI",
                                  mixcomps = c("Ucuuba", "Andiroba", "Brijo10"),
                                  model = 2)
summary(mangaandiroba_PDI_quad)

anova(ucuubaandiroba_PDI_quad)

plot(ucuubaandiroba_PDI_quad)

check_model(ucuubaandiroba_PDI_quad)

plot(simulateResiduals(ucuubaandiroba_PDI_quad))


ucuubaandiroba$PDI_sqrt<-log(ucuubaandiroba$PDI)


ucuubaandiroba_PDI_fullcub <- MixModel(frame = ucuubaandiroba, 
                                     response = "PDI",
                                     mixcomps = c("Ucuuba", "Andiroba", "Brijo10"),
                                     model = 3)
summary(ucuubaandiroba_PDI_fullcub)

plot(ucuubaandiroba_PDI_lin)

AIC(mangaandiroba_PDI_fullcub)

mangaandiroba_PDI_fullcub$terms

ucuubaandiroba_PDI_splcub <- MixModel(frame = ucuubaandiroba, 
                                    response = "PDI",
                                    mixcomps = c("Ucuuba", "Andiroba", "Brijo10"),
                                    model = 4)
summary(ucuubaandiroba_PDI_splcub)

plot(ucuubaandiroba_PDI_splcub)

check_model(ucuubaandiroba_PDI_splcub)

AIC(ucuubaandiroba_PDI_lin, ucuubaandiroba_PDI_quad,
    ucuubaandiroba_PDI_splcub, ucuubaandiroba_PDI_fullcub)

BIC(ucuubaandiroba_PDI_lin, ucuubaandiroba_PDI_quad,
    ucuubaandiroba_PDI_splcub, ucuubaandiroba_PDI_fullcub)

anova(ucuubaandiroba_PDI_lin, ucuubaandiroba_PDI_quad)
anova(ucuubaandiroba_PDI_lin, ucuubaandiroba_PDI_splcub)
anova(ucuubaandiroba_PDI_lin, ucuubaandiroba_PDI_fullcub)


anova(murumurujambu_PDI_lin, murumurujambu_PDI_splcub)
anova(ucuubaandiroba_PDI_lin, ucuubaandiroba_PDI_fullcub)

library(performance)
check_model(ucuubaandiroba_PDI_lin)
check_model(ucuubaandiroba_PDI_fullcub)

library(DHARMa)
plot(simulateResiduals(ucuubaandiroba_PDI_lin))
plot(simulateResiduals(ucuubaandiroba_PDI_fullcub))


anova(tucumatucuma_PDI_lin,tucumatucuma_PDI_quad)

bgtest(ucuubaandiroba_PDI_lin)
bgtest(ucuubaandiroba_PDI_quad)
bgtest(ucuubaandiroba_PDI_fullcub)

hatvalues(ucuubaandiroba_PDI_fullcub)
n <- nrow(ucuubaandiroba)
p <- length(coef(ucuubaandiroba_PDI_fullcub))
cutoff <- 2 * p / n
cutoff
cutoffxt <- 3 * p / n
cutoffxt



plot(mangaandiroba_PDI_lin)

coeftest(mangapracaxi_PDI_fullcub, vcov. = vcovHC(mangapracaxi_PDI_fullcub,
                                                  type = "HC1"))

plot(predict(ucuubaandiroba_PDI_lin), rstudent(ucuubaandiroba_PDI_lin),
     xlab=expression( paste("predicted response ", hat(y)) ) ,ylim=c(-5,5)   ,
     ylab="studentized residuals",type="n",
     main=expression(paste("(",frac(paste("y - ", hat(y)),sigma), ") ~ ", 
                           hat(y)   )) ) 

text(predict(ucuubaandiroba_PDI_lin), rstudent(ucuubaandiroba_PDI_lin),
     label=1:length(ucuubaandiroba$PDI)) 
abline(h=c(-2,0,2),lty=c(2,1,2)) 
grid() 

murumuruburiti

ModelPlot(model = ucuubaandiroba_PDI_lin,
          res = 500,
          dimensions = list(x1 = "Ucuuba", x2 = "Andiroba", x3 = "Brijo10"),
          lims = c(0.166, 0.5, 0.166, 0.5, 0.166, 0.5835),
          #at = c(-10, 0.1, 0.2, 0.3, 0.4, 0.5, 10),
          contour = TRUE, fill = TRUE, pseudo = T, axislab.offset = 0.1,
          color.palette = heat.colors,
          axislabs = c("Ucuuba", "Andiroba", "Brijo10"),
          cornerlabs = c("", "", ""),
          grid = TRUE, grid.pars = list(lty = 0),
          main = "PDI Contour Plot")

#size

ucuubaandiroba

ucuubaandiroba[c(12, 13, 14, 15), 5] <- NA

ucuubaandiroba$Size <- as.numeric(gsub(",", "", ucuubaandiroba$Size))

ucuubaandiroba_size_lin <- MixModel(frame = ucuubaandiroba, 
                                  response = "Size",
                                  mixcomps = c("Ucuuba", "Andiroba", "Brijo10"),
                                  model = 1)
summary(ucuubaandiroba_size_lin)

anova(ucuubaandiroba_size_lin)

plot(ucuubaandiroba_size_lin)

check_model(ucuubaandiroba_size_lin)

plot(simulateResiduals(ucuubaandiroba_size_lin))

ucuubaandiroba_size_quad <- MixModel(frame = ucuubaandiroba, 
                                   response = "Size",
                                   mixcomps = c("Ucuuba", "Andiroba", "Brijo10"),
                                   model = 2)
summary(ucuubaandiroba_size_quad)

anova(ucuubaandiroba_size_quad)

plot(ucuubaandiroba_size_quad)

check_model(ucuubaandiroba_size_quad)

plot(simulateResiduals(ucuubaandiroba_size_quad))

ucuubaandiroba_size_fullcub <- MixModel(frame = ucuubaandiroba, 
                                      response = "Size",
                                      mixcomps = c("Ucuuba", "Andiroba", "Brijo10"),
                                      model = 3)
summary(ucuubaandiroba_size_fullcub)

anova(ucuubaandiroba_size_fullcub)

plot(ucuubaandiroba_size_fullcub)

check_model(ucuubaandiroba_size_quad)

plot(simulateResiduals(ucuubaandiroba_size_fullcub))

AIC(mangaandiroba_PDI_fullcub)

mangaandiroba_PDI_fullcub$terms

ucuubaandiroba_size_splcub <- MixModel(frame = ucuubaandiroba, 
                                     response = "Size",
                                     mixcomps = c("Ucuuba", "Andiroba", "Brijo10"),
                                     model = 4)
summary(ucuubaandiroba_size_splcub)

plot(ucuubaandiroba_size_splcub)

check_model(ucuubaandiroba_size_splcub)

plot(simulateResiduals(ucuubaandiroba_size_splcub))

AIC(ucuubaandiroba_size_lin, ucuubaandiroba_size_quad,
    ucuubaandiroba_size_splcub, ucuubaandiroba_size_fullcub)

BIC(ucuubaandiroba_size_lin, ucuubaandiroba_size_quad,
    ucuubaandiroba_size_splcub, ucuubaandiroba_size_fullcub)

anova(ucuubaandiroba_size_lin, ucuubaandiroba_size_quad)
anova(ucuubaandiroba_size_quad, ucuubaandiroba_size_splcub)
anova(tucumaurucum_size_lin, tucumaurucum_size_fullcub)

anova(ucuubaandiroba_size_lin, ucuubaandiroba_size_fullcub)


anova(mangaandiroba_PDI_lin,mangaandiroba_PDI_fullcub)

bgtest(tucumaurucum_size_lin)
bgtest(ucuubaandiroba_size_quad)

plot(mangaandiroba_PDI_lin)

hatvalues(tucumatucuma_size_quad)
n <- nrow(tucumatucuma)
p <- length(coef(tucumatucuma_size_quad))
cutoff <- 2 * p / n
cutoff
cutoffxt <- 3 * p / n
cutoffxt

coeftest(tucumatucuma_size_quad, vcov. = vcovHC(tucumatucuma_size_quad,
                                                type = "HC3"))

plot(predict(ucuubaandiroba_size_quad), rstudent(ucuubaandiroba_size_quad),
     xlab=expression( paste("predicted response ", hat(y)) ) ,ylim=c(-5,5)   ,
     ylab="studentized residuals",type="n",
     main=expression(paste("(",frac(paste("y - ", hat(y)),sigma), ") ~ ", hat(y)   )) ) 

text(predict(ucuubaandiroba_size_quad), rstudent(ucuubaandiroba_size_quad),
     label=1:length(ucuubaandiroba$size)) 
abline(h=c(-2,0,2),lty=c(2,1,2)) 
grid() 

murumuruacai

ModelPlot(model = ucuubaandiroba_size_quad,
          res = 500,
          dimensions = list(x1 = "Ucuuba", x2 = "Andiroba", x3 = "Brijo10"),
          lims = c(0.166, 0.5, 0.166, 0.5, 0.166, 0.5835),
          at = c(-1000, 200, 400, 600, 800, 1000, 1200, 1400, 1600, 10000),
          contour = TRUE, fill = TRUE, pseudo = T, axislab.offset = 0.1,
          color.palette = heat.colors,
          axislabs = c("Ucuuba", "Andiroba", "Brijo10"),
          cornerlabs = c("", "", ""),
          grid = TRUE, grid.pars = list(lty = 0),
          main = "Size Contour Plot")
