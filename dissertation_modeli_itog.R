#------------------------------------Пакеты-------------------------------------
#установка пакетов
install.packages("plm") #для панельных данных
install.packages("stargazer") #для таблиц сравнения регрессий
install.packages("MASS") #для порядкового логита/пробита
install.packages("corrplot") #для корреляционных матриц
install.packages("Hmisc") #для корреляционных матриц
install.packages("car") #для бутстрапа, для графиков
install.packages("spdep") #моделирование пространственных взаимосвязей
install.packages("lattice") #отрисовка графиков
install.packages("splm") #пространственные модели на панельных данных
install.packages("spatialreg") #пространственные модели
install.packages("clusterSim")#кластерный анализ
install.packages("psych") #Для описательных статистик по группам
install.packages("factoextra") #для визуализации кластерного анализа

#активация пакетов
library("plm") 
library("stargazer") 
library("MASS")
library("corrplot")
library("Hmisc") 
library("car") 
library("spdep")
library("lattice")
library("splm")
library("spatialreg")
library("clusterSim")
library("psych")
library("factoextra")
#------------------------------------Полезные функции-------------------------------------

# function to calculate corrected SEs for OLS regression 
cse = function(reg) {
  rob = sqrt(diag(vcovHC(reg, type = "HC1")))
  return(rob)
}

# clustered SEs, clustered on "group"... could also cluster on "time" 
# compute Stata-like degrees of freedom adjustment for number of groups
# See http://www.richard-bluhm.com/clustered-ses-in-r-and-stata-2/

clse = function(reg) { 
  # index(reg, "id") returns the id or entity variable vector 
  G = length(unique(index(reg,"id")))
  N = length(index(reg,"id"))
  dfa = (G/(G - 1))   # note Bluhm multiplies this by finite-sample df adjustment
  rob = sqrt(diag(dfa*vcovHC(reg, method="arellano", type = "HC1", 
                             cluster = "group")))
  return(rob)
}

#Для создания сбалансированной панели

balanced <- function(data, ID, TIME, VARS, required=c("all", "shared")){
  if(is.character(ID)) {
    ID <- match(ID, names(data))
  }
  if(is.character(TIME)) {
    TIME <- match(TIME, names(data))
    if(missing(VARS)) { 
      VARS <- setdiff(1:ncol(data), c(ID,TIME))
    } else if (is.character(VARS)) {
      VARS <- match(VARS, names(data))
    }
    required <- match.arg(required)
    idf <- do.call(interaction, c(data[, ID, drop=FALSE], drop=TRUE))
    timef <- do.call(interaction, c(data[, TIME, drop=FALSE], drop=TRUE))
    complete <- complete.cases(data[, VARS])
    tbl <- table(idf[complete], timef[complete])
    if (required == "all") {
      keep <- which(rowSums(tbl == 1) == ncol(tbl))
      idx <- as.numeric(idf) %in% keep
    } else if (required == "shared") {
      keep <- which(colSums(tbl == 1) == nrow(tbl))
      idx <- as.numeric(timef) %in% keep
    }
    data[idx, ]
  }
}

# AIC and BIC function for splm object дописала для pooled и RE
# Взято отсюда: https://github.com/rfsaldanha/ecoespacialunicamp/blob/master/AICsplm.R

AICsplm = function(object, k=2, criterion=c("AIC", "BIC")){ 
  sp = summary(object)
  l = sp$logLik #логарифм правдоподобия
  np = length(coef(sp)) #кол-во коэффициентов
  N = nrow(sp$model) #кол-во наблюдений
  if (is.null(sp$effects)==TRUE) {
    T = length(sp$res.eff[[1]]$res.tfe)
    np = np
  }
  else if (sp$effects=="sptpfe") {
    n = length(sp$res.eff[[1]]$res.sfe) 
    T = length(sp$res.eff[[1]]$res.tfe) 
    np = np+n+T
  }
  else if (sp$effects=="spfe") {
    n = length(sp$res.eff[[1]]$res.sfe)
    np = np+n+1 
  }
  else if (sp$effects=="tpfe") {
    T = length(sp$res.eff[[1]]$res.tfe)
    np = np+T+1
  }
  if (criterion=="AIC"){
    aic = -2*l+k*np
    names(aic) <- "AIC"
    return(aic)
  }
  if (criterion=="BIC"){
    bic = -2*l+log(N)*np
    names(bic) <- "BIC"
    return(bic)
  } 
}
#------------------------Корреляции для обоснования зависимой переменной-----------------

USA_data <- read.csv("USA_data.csv", sep=";", dec=",", header=TRUE, encoding = "UTF-8")

cor.test(USA_data$ART_children, USA_data$Multiple_birth)

#------------------------Корреляции для обоснования зависимой переменной-----------------

data_full <- read.csv("data_full.csv", sep=";", dec=",", header=TRUE, encoding = "UTF-8")

cor.test(data_full$Clinics, data_full$total_mult_birth_Rosstat)

#-----------------------------------------------Корреляция по клиникам-----------------------------------

data_clinics_corr <- read.csv("data_clinics_corr.csv", sep=";", dec=",", header=TRUE)

df_all_info <- data.frame(data_clinics_corr[2:4], data_clinics_corr[6:7])
df_eff <- data.frame(data_clinics_corr[2:7])
df_price <- data.frame(data_clinics_corr[2:4], data_clinics_corr[6:8])
df_all_in <- data.frame(data_clinics_corr[2:8])

nc_1 <- complete.cases(df_all_info)
nc_2 <- complete.cases(df_eff)
nc_3 <- complete.cases(df_price)
nc_4 <- complete.cases(df_all_in)

corr_matr1 <- cor(df_all_info[nc_1,])
corr_matr2 <- cor(df_eff[nc_2,])
corr_matr3 <- cor(df_price[nc_3,])
corr_matr4 <- cor(df_all_in[nc_4,])

pc1 <- rcorr(as.matrix(df_all_info)) 
corrplot(corr_matr1, method = 'color', tl.col = 'black', type = 'lower', addCoef.col = 'black',
         tl.srt = 90, diag = TRUE, is.corr = FALSE, tl.cex = 1.5, cl.cex = 1.2, number.cex = 1.5,
         p.mat = pc1$P, sig.level = 0.01, pch.cex = 10)

pc2 <- rcorr(as.matrix(df_eff)) 
corrplot(corr_matr2, method = 'color', tl.col = 'black', type = 'lower', addCoef.col = 'black',
         tl.srt = 90, diag = TRUE, is.corr = FALSE, tl.cex = 1.5, cl.cex = 1.2, number.cex = 1.5,
         p.mat = pc2$P, sig.level = 0.01, pch.cex = 10)

pc3 <- rcorr(as.matrix(df_price)) 
corrplot(corr_matr3, method = 'color', tl.col = 'black', type = 'lower', addCoef.col = 'black',
         tl.srt = 90, diag = TRUE, is.corr = FALSE, tl.cex = 1.5, cl.cex = 1.2, number.cex = 1.5,
         p.mat = pc3$P, sig.level = 0.01, pch.cex = 10)

pc4 <- rcorr(as.matrix(df_all_in)) 
corrplot(corr_matr4, method = 'color', tl.col = 'black', type = 'lower', addCoef.col = 'black',
         tl.srt = 90, diag = TRUE, is.corr = FALSE, tl.cex = 1.5, cl.cex = 1.2, number.cex = 1.5,
         p.mat = pc4$P, sig.level = 0.01, pch.cex = 10)

#--------------------------------На лаге------------------------------------

data_lag <- read.csv("data_lag.csv", sep=";", dec=",", header=TRUE, encoding = "UTF-8")

#------------------------------------Предмодельный анализ-------------------------------------

#Диаграммы разброса и распределения
df_pl <- data.frame(data_lag[4],data_lag[15], data_lag[18], data_lag[20:22], data_lag[27:28], 
                    data_lag[31:33])

scatterplotMatrix(df_pl[1:85,], smooth=FALSE, regLine=FALSE, col='darkgrey', 
                  diagonal=list(method ="histogram", breaks="FD"))

#Создание корреляционной матрицы
df <- data.frame(data_lag[4:7], data_lag[9:11],data_lag[15], data_lag[18], data_lag[20:22],
                 data_lag[23], data_lag[25:29], data_lag[31:41])
net_cases <- complete.cases(df)
corr_matr_data_lag <- cor(df[net_cases, ])

corrplot(corr_matr_data_lag, method = 'color', tl.col = 'black', type = 'lower', addCoef.col = 'black',
         tl.srt = 90, diag = TRUE, is.corr = FALSE, tl.cex = 0.8, cl.cex = 0.8, number.cex = 0.5)

#Добавим уровень значимости по Пирсону
pearson_lag <- rcorr(as.matrix(df)) 
corrplot(corr_matr_data_lag, method = 'color', tl.col = 'black', type = 'lower', addCoef.col = 'black',
         tl.srt = 90, diag = TRUE, is.corr = FALSE, tl.cex = 0.8, cl.cex = 0.8, number.cex = 0.5,
         p.mat = pearson_lag$P, sig.level = 0.01)

#Описательная статистика
stargazer(data_lag, type="text", digits=2, summary.stat = c("n", "min", "max", "mean", "sd"))
#Статистика по годам
describeBy(data_lag, group = data_lag$Year, digits=2)
#------------------------------------Модели-------------------------------------

# Пулд МНК

Pooled_OLS_ml1 = plm(pct_IVF_Children ~ Ins + Ins_Inc + Ins_Ineq +
                      l1_Clinics + l1_pct_Spec_Clinics + l1_pct_Gov_Clinics + Year1Clinic_NA +  
                      l1_WM_Infertility + l1_pct_WM_older35 + l1_pct_WM_older35_2 + l1_Koef_Marriage +
                       log(l1_Aver_Real_Income) + l1_Inequal + Recession + 
                       l1_pct_ARTsities + l1_pct_High_Educ + l1_pct_Internet + 
                      l1_Doctors + l1_Beds + l1_WM_Consult + 
                      pct_Christianity + pct_Islam + pct_Atheist, 
                    data = data_lag, index = c("Region","Year"), model="pooling")

Pooled_OLS_ml2 = plm(pct_IVF_Children ~ Ins + Ins_Inc + Ins_Ineq +
                       l1_Clinics + l1_pct_Spec_Clinics + l1_pct_Gov_Clinics + Year1Clinic +  
                       l1_WM_Infertility + l1_pct_WM_older35 + l1_pct_WM_older35_2 + l1_Koef_Marriage + 
                       log(l1_Aver_Real_Income) + l1_Inequal + Recession + 
                       l1_pct_ARTsities + l1_pct_High_Educ + l1_pct_Internet + 
                       l1_Doctors + l1_Beds + l1_WM_Consult + 
                       pct_Christianity + pct_Islam + pct_Atheist, 
                     data = data_lag, index = c("Region","Year"), model="pooling")

Pooled_OLS_ml3 = plm(pct_IVF_Children ~ Ins + Ins_Inc + Ins_Ineq +
                       l1_Clinics + l1_pct_Spec_Clinics + l1_pct_Gov_Clinics +
                       l1_WM_Infertility + l1_pct_WM_older35 + l1_pct_WM_older35_2 + l1_Koef_Marriage + 
                       log(l1_Aver_Real_Income) + l1_Inequal + Recession + 
                       l1_pct_ARTsities + l1_pct_High_Educ + l1_pct_Internet + 
                       l1_Doctors + l1_Beds + l1_WM_Consult, 
                    data = data_lag, index = c("Region","Year"), model="pooling")

# Модель с фиксированными эффектами

Fixed_effects_ml1 = plm(pct_IVF_Children ~ Ins + Ins_Inc + Ins_Ineq + 
                          l1_Clinics + l1_pct_Spec_Clinics + l1_pct_Gov_Clinics + 
                          l1_WM_Infertility + l1_pct_WM_older35 + l1_pct_WM_older35_2 + l1_Koef_Marriage + 
                          log(l1_Aver_Real_Income) + l1_Inequal + Recession + 
                          l1_pct_ARTsities + l1_pct_High_Educ + l1_pct_Internet + 
                          l1_Doctors + l1_Beds + l1_WM_Consult, 
                       data = data_lag, index = c("Region","Year"), model="within", 
                       effect="time")

Fixed_effects_ml2 = plm(pct_IVF_Children ~ Ins + Ins_Inc + Ins_Ineq +
                          l1_Clinics + l1_pct_Spec_Clinics + l1_pct_Gov_Clinics + 
                          l1_WM_Infertility + l1_pct_WM_older35 + l1_Koef_Marriage + 
                          log(l1_Aver_Real_Income) + l1_Inequal + Recession + 
                          l1_pct_ARTsities + l1_pct_High_Educ + l1_pct_Internet + 
                          l1_Doctors + l1_Beds + l1_WM_Consult, 
                       data = data_lag, index = c("Region","Year"), model="within", 
                       effect="individual")

#модель со случайными эффектами

Random_effects_ml1 = plm(pct_IVF_Children ~ Ins + Ins_Inc + Ins_Ineq +
                           l1_Clinics + l1_pct_Spec_Clinics + l1_pct_Gov_Clinics +  
                           l1_WM_Infertility + l1_pct_WM_older35 + l1_pct_WM_older35_2 + l1_Koef_Marriage + 
                           log(l1_Aver_Real_Income) + l1_Inequal + Recession + 
                           l1_pct_ARTsities + l1_pct_High_Educ + l1_pct_Internet + 
                           l1_Doctors + l1_Beds + l1_WM_Consult, 
                        data = data_lag, index = c("Region","Year"), model="random")

#----------------------------------------все модели------------------------------------

stargazer(Pooled_OLS_ml1, Pooled_OLS_ml2, Pooled_OLS_ml3, 
          Fixed_effects_ml1, Fixed_effects_ml2, Random_effects_ml1,
          se=list(clse(Pooled_OLS_ml1), clse(Pooled_OLS_ml2), clse(Pooled_OLS_ml3), clse(Fixed_effects_ml1), 
                  clse(Fixed_effects_ml2), clse(Random_effects_ml1)),
          title="Panel models, 1 year lag, full sample", type="text", 
          column.labels=c("PO all full", "PO all restr", "PO stad", "FE time eff", "FE ind eff", "RE"), 
          df=FALSE, digits=2)

#------------------------------тесты на сравнение моделей-----------------------

#Тест на линейное ограничение
pFtest(Fixed_effects_ml1, Pooled_OLS_ml3)
#H0 - пулд, здесь p-value маленькое, следовательно, H1 -> FE
pFtest(Fixed_effects_ml2, Pooled_OLS_ml3)
#H0 - пулд, здесь p-value маленькое, следовательно, H1 -> FE

#тест Бреуша - Пагана для сравнивнения пулд и RE 
plmtest(Pooled_OLS_ml3,effect="time",type="bp")
#H0 отсутствие индивидуальных эффектов. p-value маленькое значит H1 лучше RE
#тест Бреуша - Пагана для сравнивнения пулд и RE 
plmtest(Pooled_OLS_ml3,effect="individual",type="bp")
#H0 отсутствие индивидуальных эффектов. p-value маленькое значит H1 лучше RE

#Тест Хаусмана для сравнения FE и RE моделей
phtest(Random_effects_ml1, Fixed_effects_ml1)
#H0 есть индивидуальные случайные эффекты, p-value маленькое, значит лучше FE (5% уровень значимости)
phtest(Random_effects_ml1, Fixed_effects_ml2)
#H0 нет индивидуальных случайных эффектов, p-value маленькое, значит лучше FE

#В общем лучше FE

#------------------------------альтернативные метрики-----------------------

#с индивидуальными эффектами
Fixed_effects_ml3 = plm(pct_IVF_Children ~ Ins + Ins_Inc + Ins_Ineq + 
                          l1_Clinics_RAHR + l1_pct_Spec_Clinics + l1_pct_Gov_Clinics +  
                          l1_WM_Infertility + l1_pct_WM_older35 + l1_Koef_Marriage + 
                          log(l1_Aver_Real_Income) + l1_Inequal + Recession + 
                          l1_pct_ARTsities_RAHR + l1_pct_High_Educ + l1_pct_Internet + 
                          l1_Doctors + l1_Beds + l1_WM_Consult, 
                       data = data_lag, index = c("Region","Year"), model="within", 
                       effect="individual")

Fixed_effects_ml4 = plm(pct_IVF_Children ~ Ins + Ins_Inc + Ins_Ineq + 
                          l1_Clinics + l1_pct_Spec_Clinics + l1_pct_Gov_Clinics +  
                          l1_Total_WM_Infert + l1_pct_WM_older35 + l1_Koef_Marriage + 
                          log(l1_Aver_Real_Income) + l1_Inequal + Recession + 
                          l1_pct_ARTsities + l1_pct_High_Educ + l1_pct_Internet + 
                          l1_Doctors + l1_Beds + l1_WM_Consult, 
                       data = data_lag, index = c("Region","Year"), model="within", 
                       effect="individual")

Fixed_effects_ml5 = plm(pct_IVF_Children ~ Ins + Ins_Inc + Ins_Ineq + 
                          l1_Clinics + l1_pct_Spec_Clinics + l1_pct_Gov_Clinics +  
                          l1_WM_Infertility + l1_pct_WM_older35 + l1_Koef_Marriage + 
                          log(l1_Aver_Real_Income) + l1_Inequal + Recession + 
                          l1_pct_ARTsities + l1_pct_High_Educ + l1_pct_Internet + 
                          l1_Doctors + l1_Beds + l1_WM_Consult + l1_Reiting_Mother_Health, 
                       data = data_lag, index = c("Region","Year"), model="within", 
                       effect="individual")

#Модель 5 ограничивывает выборку снизу: количество лет до включения

Fixed_effects_ml6 = plm(pct_IVF_Children ~ Ins + Ins_Inc + Ins_Ineq + 
                          l1_Clinics + l1_pct_Spec_Clinics + l1_pct_Gov_Clinics +  
                          l1_Air_Polution	+ l1_Water_Polution	+ l1_Industry_Ind	 + 
                          l1_pct_WM_older35 + l1_Koef_Marriage + 
                          log(l1_Aver_Real_Income) + l1_Inequal + Recession + 
                          l1_pct_ARTsities + l1_pct_High_Educ + l1_pct_Internet + 
                          l1_Doctors + l1_Beds + l1_WM_Consult, 
                       data = data_lag, index = c("Region","Year"), model="within", 
                       effect="individual")

#Модель 6 включает перетягивающие в разные стороны эффект переменные

Fixed_effects_ml7 = plm(pct_IVF_Children ~ Ins + Ins_Inc + Ins_Ineq + 
                          l1_Clinics + l1_pct_Spec_Clinics + l1_pct_Gov_Clinics +  
                          Ecology	 + 
                          l1_pct_WM_older35 + l1_Koef_Marriage + 
                          log(l1_Aver_Real_Income) + l1_Inequal + Recession + 
                          l1_pct_ARTsities + l1_pct_High_Educ + l1_pct_Internet + 
                          l1_Doctors + l1_Beds + l1_WM_Consult, 
                        data = data_lag, index = c("Region","Year"), model="within", 
                        effect="individual")

Fixed_effects_ml8 = plm(pct_IVF_Children ~ Ins + Ins_Inc + Ins_Ineq + 
                          l1_Clinics + l1_pct_Spec_Clinics + l1_pct_Gov_Clinics +  
                          l1_WM_Infertility + l1_pct_WM_older35 + l1_Koef_Marriage + 
                          log(l1_Aver_Real_Income) + l1_Inequal + Recession + 
                          l1_pct_Urban + l1_pct_High_Educ + l1_pct_Internet + 
                          l1_Doctors + l1_Beds + l1_WM_Consult, 
                        data = data_lag, index = c("Region","Year"), model="within", 
                        effect="individual")

stargazer(Fixed_effects_ml2, Fixed_effects_ml3, Fixed_effects_ml4, Fixed_effects_ml5, Fixed_effects_ml6,
          Fixed_effects_ml7, Fixed_effects_ml8,
          se=list(clse(Fixed_effects_ml2), clse(Fixed_effects_ml3), clse(Fixed_effects_ml4), 
                  clse(Fixed_effects_ml5), clse(Fixed_effects_ml6), clse(Fixed_effects_ml7), clse(Fixed_effects_ml8)),
          title="Fixed individual effects: different variables, 1 year lag", type="text", 
          column.labels=c("st", "RAHR", "Total WM infert", "Reiting", "Ecology 3", "Ecology", "Urban"), 
          df=FALSE, digits=2)

#------------------------------До и после включения в ОМС-----------------------

Fixed_effects_ml9 = plm(pct_IVF_Children ~ 
                           l1_Clinics + l1_pct_Spec_Clinics + l1_pct_Gov_Clinics + 
                           l1_WM_Infertility + l1_pct_WM_older35 + l1_Koef_Marriage + 
                           log(l1_Aver_Real_Income) + l1_Inequal + Recession + 
                           l1_pct_ARTsities + l1_pct_High_Educ + l1_pct_Internet + 
                           l1_Doctors + l1_Beds + l1_WM_Consult, 
                         data = data_lag[254:585,], index = c("Region","Year"), model="within", 
                         effect="individual")

Fixed_effects_ml10 = plm(pct_IVF_Children ~ 
                           l1_Clinics + l1_pct_Spec_Clinics + l1_pct_Gov_Clinics + 
                           l1_WM_Infertility + l1_pct_WM_older35 + l1_Koef_Marriage + 
                           log(l1_Aver_Real_Income) + l1_Inequal + Recession + 
                           l1_pct_ARTsities + l1_pct_High_Educ + l1_pct_Internet + 
                           l1_Doctors + l1_Beds + l1_WM_Consult, 
                         data = data_lag[1:253,], index = c("Region","Year"), model="within", 
                         effect="individual")


stargazer(Fixed_effects_ml2, Fixed_effects_ml9, Fixed_effects_ml10,
          se=list(clse(Fixed_effects_ml2), clse(Fixed_effects_ml9),clse(Fixed_effects_ml10)),
          title="Panel models, 1 year lag, different time period", type="text", 
          column.labels=c("FE, ind, stand", "FE, ind, bef2014", "FE, ind, aft2014"), 
          df=FALSE, digits=2)

#--------------------------------пространственные матрицы------------------------------------

matr_dist_simpl1 <- read.csv("matrix_simpl_distance_v1.csv", sep=";", dec=",", header=TRUE, encoding = "UTF-8")
matr_dist_simpl_1 <- as.matrix(matr_dist_simpl1[2:86])
matr_dist_simpl__1 = mat2listw(matr_dist_simpl_1, row.names = matr_dist_simpl1[,1], style="W")

matr_dist_railw1 <- read.csv("matrix_railway_v1.csv", sep=";", dec=",", header=TRUE, encoding = "UTF-8")
matr_dist_railw_1 <- as.matrix(matr_dist_railw1[2:86])
matr_dist_railw__1 = mat2listw(matr_dist_railw_1, row.names = matr_dist_railw1[,1], style="W")

matr_dist_railw2 <- read.csv("matrix_railway_v2.csv", sep=";", dec=",", header=TRUE, encoding = "UTF-8")
matr_dist_railw_2 <- as.matrix(matr_dist_railw2[2:86])
matr_dist_railw__2 = mat2listw(matr_dist_railw_2, row.names = matr_dist_railw2[,1], style="W")

matr_dist_railw3 <- read.csv("matrix_railway_v3.csv", sep=";", dec=",", header=TRUE, encoding = "UTF-8")
matr_dist_railw_3 <- as.matrix(matr_dist_railw3[2:86])
matr_dist_railw__3 = mat2listw(matr_dist_railw_3, row.names = matr_dist_railw3[,1], style="W")

matr_dist_railw4 <- read.csv("matrix_railway_v4.csv", sep=";", dec=",", header=TRUE, encoding = "UTF-8")
matr_dist_railw_4 <- as.matrix(matr_dist_railw4[2:86])
matr_dist_railw__4 = mat2listw(matr_dist_railw_4, row.names = matr_dist_railw4[,1], style="W")

matr_dist_railw5 <- read.csv("matrix_railway_v5.csv", sep=";", dec=",", header=TRUE, encoding = "UTF-8")
matr_dist_railw_5 <- as.matrix(matr_dist_railw5[2:86])
matr_dist_railw__5 = mat2listw(matr_dist_railw_5, row.names = matr_dist_railw5[,1], style="W")

matr_clinics_d <- read.csv("matrix_clinics_v2.csv", sep=";", dec=",", header=TRUE, encoding = "UTF-8")
matr_clinics <- as.matrix(matr_clinics_d[2:86])

matr_dist_clinics1 <- matr_clinics*matr_dist_simpl_1
matr_dist_clinics__1 = mat2listw(matr_dist_clinics1, row.names = matr_dist_simpl1[,1], style="W")
matr_dist_clinics2 <- matr_clinics*matr_dist_railw_1
matr_dist_clinics__2 = mat2listw(matr_dist_clinics2, row.names = matr_dist_railw1[,1], style="W")
matr_dist_clinics3 <- matr_clinics*matr_dist_railw_2
matr_dist_clinics__3 = mat2listw(matr_dist_clinics3, row.names = matr_dist_railw2[,1], style="W")
matr_dist_clinics4 <- matr_clinics*matr_dist_railw_3
matr_dist_clinics__4 = mat2listw(matr_dist_clinics4, row.names = matr_dist_railw3[,1], style="W")
matr_dist_clinics5 <- matr_clinics*matr_dist_railw_4
matr_dist_clinics__5 = mat2listw(matr_dist_clinics5, row.names = matr_dist_railw4[,1], style="W")
matr_dist_clinics6 <- matr_clinics*matr_dist_railw_5
matr_dist_clinics__6 = mat2listw(matr_dist_clinics6, row.names = matr_dist_railw5[,1], style="W")

ramp2 = colorRampPalette(c("white","blue"))
levels2 = 1 / 1:50 # шкала 1, 0.5, 0.33, 0.25 ... 0.1

levelplot(listw2mat(matr_dist_clinics__1), main = "Матрица весов1",at = levels2, col.regions = ramp2(10))
levelplot(listw2mat(matr_dist_clinics__2), main = "Матрица весов2",at = levels2, col.regions = ramp2(10))
levelplot(listw2mat(matr_dist_clinics__3), main = "Матрица весов3",at = levels2, col.regions = ramp2(10))
levelplot(listw2mat(matr_dist_clinics__4), main = "Матрица весов4",at = levels2, col.regions = ramp2(10))
levelplot(listw2mat(matr_dist_clinics__5), main = "Матрица весов5",at = levels2, col.regions = ramp2(10))
levelplot(listw2mat(matr_dist_clinics__6), main = "Матрица весов6",at = levels2, col.regions = ramp2(10))

#Вспомогательная матрица без Крыма и Севастополя
matr_dist_simpl1_ks <- read.csv("matrix_simpl_distance_v1-ks.csv", sep=";", dec=",", header=TRUE, encoding = "UTF-8")
matr_dist_simpl_1_ks <- as.matrix(matr_dist_simpl1_ks[2:84])
matr_dist_simpl__1_ks = mat2listw(matr_dist_simpl_1_ks, row.names = matr_dist_simpl1_ks[,1], style="W")
matr_dist_railw1_ks <- read.csv("matrix_railway_v1-ks.csv", sep=";", dec=",", header=TRUE, encoding = "UTF-8")
matr_dist_railw_1_ks <- as.matrix(matr_dist_railw1_ks[2:84])
matr_dist_railw__1_ks = mat2listw(matr_dist_railw_1_ks, row.names = matr_dist_railw1_ks[,1], style="W")
matr_dist_railw4_ks <- read.csv("matrix_railway_v4-ks.csv", sep=";", dec=",", header=TRUE, encoding = "UTF-8")
matr_dist_railw_4_ks <- as.matrix(matr_dist_railw4_ks[2:84])
matr_dist_railw__4_ks = mat2listw(matr_dist_railw_4_ks, row.names = matr_dist_railw4_ks[,1], style="W")
matr_dist_railw5_ks <- read.csv("matrix_railway_v5-ks.csv", sep=";", dec=",", header=TRUE, encoding = "UTF-8")
matr_dist_railw_5_ks <- as.matrix(matr_dist_railw5_ks[2:84])
matr_dist_railw__5_ks = mat2listw(matr_dist_railw_5_ks, row.names = matr_dist_railw5_ks[,1], style="W")
matr_clinics_d_ks <- read.csv("matrix_clinics_v2-ks.csv", sep=";", dec=",", header=TRUE, encoding = "UTF-8")
matr_clinics_ks <- as.matrix(matr_clinics_d_ks[2:84])
matr_dist_clinics1_ks <- matr_clinics_ks*matr_dist_simpl_1_ks
matr_dist_clinics__1_ks = mat2listw(matr_dist_clinics1_ks, row.names = matr_dist_simpl1_ks[,1], style="W")
matr_dist_clinics2_ks <- matr_clinics_ks*matr_dist_railw_1_ks
matr_dist_clinics__2_ks = mat2listw(matr_dist_clinics2_ks, row.names = matr_dist_railw1_ks[,1], style="W")
matr_dist_clinics5_ks <- matr_clinics_ks*matr_dist_railw_4_ks
matr_dist_clinics__5_ks = mat2listw(matr_dist_clinics5_ks, row.names = matr_dist_railw4_ks[,1], style="W")
matr_dist_clinics6_ks <- matr_clinics_ks*matr_dist_railw_5_ks
matr_dist_clinics__6_ks = mat2listw(matr_dist_clinics6_ks, row.names = matr_dist_railw5_ks[,1], style="W")

#--------------------------------Индексы Морана------------------------------------

moran.test(data_lag$pct_IVF_Children[1:85], matr_dist_clinics__1)
moran.test(data_lag$pct_IVF_Children[86:170], matr_dist_clinics__1)
moran.test(data_lag$pct_IVF_Children[171:253], matr_dist_clinics__1_ks)
moran.test(data_lag$pct_IVF_Children[254:336], matr_dist_clinics__1_ks)
moran.test(data_lag$pct_IVF_Children[337:419], matr_dist_clinics__1_ks)
moran.test(data_lag$pct_IVF_Children[420:502], matr_dist_clinics__1_ks)
moran.test(data_lag$pct_IVF_Children[503:585], matr_dist_clinics__1_ks)

moran.test(data_lag$pct_IVF_Children[1:85], matr_dist_clinics__2)
moran.test(data_lag$pct_IVF_Children[86:170], matr_dist_clinics__2)
moran.test(data_lag$pct_IVF_Children[171:253], matr_dist_clinics__2_ks)
moran.test(data_lag$pct_IVF_Children[254:336], matr_dist_clinics__2_ks)
moran.test(data_lag$pct_IVF_Children[337:419], matr_dist_clinics__2_ks)
moran.test(data_lag$pct_IVF_Children[420:502], matr_dist_clinics__2_ks)
moran.test(data_lag$pct_IVF_Children[503:585], matr_dist_clinics__2_ks)

moran.test(data_lag$pct_IVF_Children[1:85], matr_dist_clinics__5)
moran.test(data_lag$pct_IVF_Children[86:170], matr_dist_clinics__5)
moran.test(data_lag$pct_IVF_Children[171:253], matr_dist_clinics__5_ks)
moran.test(data_lag$pct_IVF_Children[254:336], matr_dist_clinics__5_ks)
moran.test(data_lag$pct_IVF_Children[337:419], matr_dist_clinics__5_ks)
moran.test(data_lag$pct_IVF_Children[420:502], matr_dist_clinics__5_ks)
moran.test(data_lag$pct_IVF_Children[503:585], matr_dist_clinics__5_ks)

moran.test(data_lag$pct_IVF_Children[1:85], matr_dist_clinics__6)
moran.test(data_lag$pct_IVF_Children[86:170], matr_dist_clinics__6)
moran.test(data_lag$pct_IVF_Children[171:253], matr_dist_clinics__6_ks)
moran.test(data_lag$pct_IVF_Children[254:336], matr_dist_clinics__6_ks)
moran.test(data_lag$pct_IVF_Children[337:419], matr_dist_clinics__6_ks)
moran.test(data_lag$pct_IVF_Children[420:502], matr_dist_clinics__6_ks)
moran.test(data_lag$pct_IVF_Children[503:585], matr_dist_clinics__6_ks)

moran.plot(data_lag$pct_IVF_Children[1:85], matr_dist_clinics__2)
moran.plot(data_lag$pct_IVF_Children[86:170], matr_dist_clinics__2)
moran.plot(data_lag$pct_IVF_Children[171:253], matr_dist_clinics__2_ks)
moran.plot(data_lag$pct_IVF_Children[254:336], matr_dist_clinics__2_ks)
moran.plot(data_lag$pct_IVF_Children[337:419], matr_dist_clinics__2_ks)
moran.plot(data_lag$pct_IVF_Children[420:502], matr_dist_clinics__2_ks)
moran.plot(data_lag$pct_IVF_Children[503:585], matr_dist_clinics__2_ks)

moran.plot(data_lag$pct_IVF_Children[1:85], matr_dist_clinics__5)
moran.plot(data_lag$pct_IVF_Children[86:170], matr_dist_clinics__5)
moran.plot(data_lag$pct_IVF_Children[171:253], matr_dist_clinics__5_ks)
moran.plot(data_lag$pct_IVF_Children[254:336], matr_dist_clinics__5_ks)
moran.plot(data_lag$pct_IVF_Children[337:419], matr_dist_clinics__5_ks)
moran.plot(data_lag$pct_IVF_Children[420:502], matr_dist_clinics__5_ks)
moran.plot(data_lag$pct_IVF_Children[503:585], matr_dist_clinics__5_ks)

moran.plot(data_lag$pct_IVF_Children[1:85], matr_dist_clinics__6)
moran.plot(data_lag$pct_IVF_Children[86:170], matr_dist_clinics__6)
moran.plot(data_lag$pct_IVF_Children[171:253], matr_dist_clinics__6_ks)
moran.plot(data_lag$pct_IVF_Children[254:336], matr_dist_clinics__6_ks)
moran.plot(data_lag$pct_IVF_Children[337:419], matr_dist_clinics__6_ks)
moran.plot(data_lag$pct_IVF_Children[420:502], matr_dist_clinics__6_ks)
moran.plot(data_lag$pct_IVF_Children[503:585], matr_dist_clinics__6_ks)
#--------------------------------пространственные модели------------------------------------
 
data_lag_restricted <- data.frame(data_lag[2:7], data_lag[9:11], data_lag[15], data_lag[18:23],
                                  data_lag[25], data_lag[27:28], data_lag[31:33])

Balanced_data_lag <- balanced(data_lag_restricted, "Region", "Year")

#Потеряли: Ивановская область (14) (женское бесплодие за 2016 г)
# Крым (54) и Севастополь (69) (огр по времени)
# Республика Ингушения (50) (женское бесплодие за 2014 г)
# Чеченская Республика (82) (неравенство 2011-2010, доступ к интернету 2010-13)

#----------------------Обрезанные матрицы для сбалансированной панели-------------------------

matr_dist_railwb1 <- read.csv("matrix_railway_v1_b.csv", sep=";", dec=",", header=TRUE, encoding = "UTF-8")
matr_dist_railw_b1 <- as.matrix(matr_dist_railwb1[2:81])
matr_dist_railw__b1 = mat2listw(matr_dist_railw_b1, row.names = matr_dist_railwb1[,1], style="W")

matr_dist_railwb2 <- read.csv("matrix_railway_v5_b.csv", sep=";", dec=",", header=TRUE, encoding = "UTF-8")
matr_dist_railw_b2 <- as.matrix(matr_dist_railwb2[2:81])
matr_dist_railw__b2 = mat2listw(matr_dist_railw_b2, row.names = matr_dist_railwb2[,1], style="W")

matr_clinicsb <- read.csv("matrix_clinics_v2_b.csv", sep=";", dec=",", header=TRUE, encoding = "UTF-8")
matr_clinicsb <- as.matrix(matr_clinicsb[2:81])

matr_dist_clinicsb1 <- matr_clinicsb*matr_dist_railw_b1
matr_dist_clinics__b1 = mat2listw(matr_dist_clinicsb1, row.names = matr_dist_railwb1[,1], style="W")
matr_dist_clinicsb2 <- matr_clinicsb*matr_dist_railw_b2
matr_dist_clinics__b2 = mat2listw(matr_dist_clinicsb2, row.names = matr_dist_railwb2[,1], style="W")

#--------------------------------пространственные модели------------------------------------
SL_c1 <- matr_dist_clinicsb1%*%Balanced_data_lag$l1_Clinics[1:80]
SL_c2 <- matr_dist_clinicsb1%*%Balanced_data_lag$l1_Clinics[81:160]
SL_c3 <- matr_dist_clinicsb1%*%Balanced_data_lag$l1_Clinics[161:240]
SL_c4 <- matr_dist_clinicsb1%*%Balanced_data_lag$l1_Clinics[241:320]
SL_c5 <- matr_dist_clinicsb1%*%Balanced_data_lag$l1_Clinics[321:400]
SL_c6 <- matr_dist_clinicsb1%*%Balanced_data_lag$l1_Clinics[401:480]
SL_c7 <- matr_dist_clinicsb1%*%Balanced_data_lag$l1_Clinics[481:560]

SL_c1 = c(SL_c1, SL_c2, SL_c3, SL_c4, SL_c5, SL_c6, SL_c7)

Balanced_data_lag$SLClinics1 <- SL_c1

SL_c11 <- matr_dist_clinicsb2%*%Balanced_data_lag$l1_Clinics[1:80]
SL_c22 <- matr_dist_clinicsb2%*%Balanced_data_lag$l1_Clinics[81:160]
SL_c33 <- matr_dist_clinicsb2%*%Balanced_data_lag$l1_Clinics[161:240]
SL_c44 <- matr_dist_clinicsb2%*%Balanced_data_lag$l1_Clinics[241:320]
SL_c55 <- matr_dist_clinicsb2%*%Balanced_data_lag$l1_Clinics[321:400]
SL_c66 <- matr_dist_clinicsb2%*%Balanced_data_lag$l1_Clinics[401:480]
SL_c77 <- matr_dist_clinicsb2%*%Balanced_data_lag$l1_Clinics[481:560]

SL_c2 = c(SL_c11, SL_c22, SL_c33, SL_c44, SL_c55, SL_c66, SL_c77)

Balanced_data_lag$SLClinics2 <- SL_c2


Spatial_Fixed_effects_m1 <- spml(pct_IVF_Children ~ Ins + Ins_Inc + Ins_Ineq +
                                   l1_Clinics + l1_pct_Spec_Clinics + l1_pct_Gov_Clinics +  
                                   l1_WM_Infertility + l1_pct_WM_older35 + l1_Koef_Marriage + 
                                   log(l1_Aver_Real_Income) + l1_Inequal + Recession + 
                                   l1_pct_ARTsities + l1_pct_High_Educ + l1_pct_Internet + 
                                   l1_Doctors + l1_Beds + l1_WM_Consult, 
                                 data = Balanced_data_lag, listw = matr_dist_clinics__b1, index = c("Region","Year"),
                                 model="within", effect = "individual", spatial.error="none", 
                                 lag=TRUE, quiet=FALSE)

summary(Spatial_Fixed_effects_m1)

Spatial_Fixed_effects_m2 <- spml(pct_IVF_Children ~ Ins + Ins_Inc + Ins_Ineq +
                                   l1_Clinics + l1_pct_Spec_Clinics + l1_pct_Gov_Clinics +  
                                   l1_WM_Infertility + l1_pct_WM_older35 + l1_Koef_Marriage + 
                                   log(l1_Aver_Real_Income) + l1_Inequal + Recession + 
                                   l1_pct_ARTsities + l1_pct_High_Educ + l1_pct_Internet + 
                                   l1_Doctors + l1_Beds + l1_WM_Consult, 
                                 data = Balanced_data_lag, listw = matr_dist_clinics__b1, index = c("Region","Year"),
                                 model="within", effect = "individual", spatial.error="b", 
                                 lag=FALSE, quiet=FALSE)

summary(Spatial_Fixed_effects_m2)

Spatial_Fixed_effects_m3 <- plm(pct_IVF_Children ~ Ins + Ins_Inc + Ins_Ineq +
                                  l1_Clinics + l1_pct_Spec_Clinics + l1_pct_Gov_Clinics +  
                                  l1_WM_Infertility + l1_pct_WM_older35 + l1_Koef_Marriage + 
                                  log(l1_Aver_Real_Income) + l1_Inequal + Recession + 
                                  l1_pct_ARTsities + l1_pct_High_Educ + l1_pct_Internet + 
                                  l1_Doctors + l1_Beds + l1_WM_Consult +
                                  SLClinics1, 
                                data = Balanced_data_lag, index = c("Region","Year"),
                                model="within", effect = "individual")

summary(Spatial_Fixed_effects_m3)

Spatial_Fixed_effects_m4 <- spml(pct_IVF_Children ~ Ins + Ins_Inc + Ins_Ineq +
                                   l1_Clinics + l1_pct_Spec_Clinics + l1_pct_Gov_Clinics +  
                                   l1_WM_Infertility + l1_pct_WM_older35 + l1_Koef_Marriage + 
                                   log(l1_Aver_Real_Income) + l1_Inequal + Recession + 
                                   l1_pct_ARTsities + l1_pct_High_Educ + l1_pct_Internet + 
                                   l1_Doctors + l1_Beds + l1_WM_Consult +
                                   SLClinics1, 
                                 data = Balanced_data_lag, listw = matr_dist_clinics__b1, index = c("Region","Year"),
                                 model="within", effect = "individual", spatial.error="none", 
                                 lag=TRUE, quiet=FALSE)

summary(Spatial_Fixed_effects_m4)

Spatial_Fixed_effects_m5 <- spml(pct_IVF_Children ~ Ins + Ins_Inc + Ins_Ineq +
                                   l1_Clinics + l1_pct_Spec_Clinics + l1_pct_Gov_Clinics +  
                                   l1_WM_Infertility + l1_pct_WM_older35 + l1_Koef_Marriage + 
                                   log(l1_Aver_Real_Income) + l1_Inequal + Recession + 
                                   l1_pct_ARTsities + l1_pct_High_Educ + l1_pct_Internet + 
                                   l1_Doctors + l1_Beds + l1_WM_Consult, 
                                 data = Balanced_data_lag, listw = matr_dist_clinics__b2, index = c("Region","Year"),
                                 model="within", effect = "individual", spatial.error="none", lag=TRUE,
                                 quiet=FALSE)

summary(Spatial_Fixed_effects_m5)

Spatial_Fixed_effects_m6 <- spml(pct_IVF_Children ~ Ins + Ins_Inc + Ins_Ineq +
                                   l1_Clinics + l1_pct_Spec_Clinics + l1_pct_Gov_Clinics +  
                                   l1_WM_Infertility + l1_pct_WM_older35 + l1_Koef_Marriage + 
                                   log(l1_Aver_Real_Income) + l1_Inequal + Recession + 
                                   l1_pct_ARTsities + l1_pct_High_Educ + l1_pct_Internet + 
                                   l1_Doctors + l1_Beds + l1_WM_Consult, 
                                 data = Balanced_data_lag, listw = matr_dist_clinics__b2, index = c("Region","Year"),
                                 model="within", effect = "individual", spatial.error="b", lag=FALSE,
                                 quiet=FALSE)

summary(Spatial_Fixed_effects_m6)

Spatial_Fixed_effects_m7 <- plm(pct_IVF_Children ~ Ins + Ins_Inc + Ins_Ineq +
                                  l1_Clinics + l1_pct_Spec_Clinics + l1_pct_Gov_Clinics +  
                                  l1_WM_Infertility + l1_pct_WM_older35 + l1_Koef_Marriage + 
                                  log(l1_Aver_Real_Income) + l1_Inequal + Recession + 
                                  l1_pct_ARTsities + l1_pct_High_Educ + l1_pct_Internet + 
                                  l1_Doctors + l1_Beds + l1_WM_Consult +
                                  SLClinics2, 
                                data = Balanced_data_lag, index = c("Region","Year"),
                                model="within", effect = "individual")

summary(Spatial_Fixed_effects_m7)

Spatial_Fixed_effects_m8<- spml(pct_IVF_Children ~ Ins + Ins_Inc + Ins_Ineq +
                                   l1_Clinics + l1_pct_Spec_Clinics + l1_pct_Gov_Clinics +  
                                   l1_WM_Infertility + l1_pct_WM_older35 + l1_Koef_Marriage + 
                                   log(l1_Aver_Real_Income) + l1_Inequal + Recession + 
                                   l1_pct_ARTsities + l1_pct_High_Educ + l1_pct_Internet + 
                                   l1_Doctors + l1_Beds + l1_WM_Consult +
                                  SLClinics2, 
                                 data = Balanced_data_lag, listw = matr_dist_clinics__b2, index = c("Region","Year"),
                                 model="within", effect = "individual", spatial.error="none", lag=TRUE,
                                 quiet=FALSE)

summary(Spatial_Fixed_effects_m8)

#Чем больше логарифм функции правдоподобия, тем лучше
summary(Spatial_Fixed_effects_m1)$logLik
summary(Spatial_Fixed_effects_m2)$logLik
summary(Spatial_Fixed_effects_m3)$logLik
summary(Spatial_Fixed_effects_m4)$logLik
summary(Spatial_Fixed_effects_m5)$logLik
summary(Spatial_Fixed_effects_m6)$logLik
summary(Spatial_Fixed_effects_m7)$logLik
summary(Spatial_Fixed_effects_m8)$logLik

#Чем меньше критерий Акаике, тем лучше
AICsplm(Spatial_Fixed_effects_m1, criterion="AIC")
AICsplm(Spatial_Fixed_effects_m2, criterion="AIC")
AICsplm(Spatial_Fixed_effects_m3, criterion="AIC")
AICsplm(Spatial_Fixed_effects_m4, criterion="AIC")
AICsplm(Spatial_Fixed_effects_m5, criterion="AIC")
AICsplm(Spatial_Fixed_effects_m6, criterion="AIC")
AICsplm(Spatial_Fixed_effects_m7, criterion="AIC")
AICsplm(Spatial_Fixed_effects_m8, criterion="AIC")

#------------------------------------данные по клиникам-------------------------------------

data_clinics <- read.csv("data_clinics.csv", sep=";", dec=",", header=TRUE, encoding = "UTF-8")
View(data_clinics)

#Описательная статистика
stargazer(data_clinics, type="text", digits=2, summary.stat = c("n", "min", "max", "mean", "sd"))

#------------------------------------Корреляционные матрицы-------------------------------------

corr_matr_clinics <- cor(data_clinics[3:23])
corrplot(corr_matr_clinics, method = 'color', tl.col = 'black', type = 'lower', addCoef.col = 'black',
         tl.srt = 90, diag = TRUE, is.corr = FALSE, tl.cex = 0.8, cl.cex = 0.8, number.cex = 0.6)

#Добавим уровень значимости по Пирсону
pearson_clinics <- rcorr(as.matrix(data_clinics[3:23])) 
corrplot(corr_matr_clinics, method = 'color', tl.col = 'black', type = 'lower', addCoef.col = 'black',
         tl.srt = 90, diag = TRUE, is.corr = FALSE, tl.cex = 0.8, cl.cex = 0.8, number.cex = 0.6,
         p.mat = pearson_clinics$P, sig.level = 0.1)

#вообще модели бессмысленны

#------------------------------------Модели-------------------------------------

Pooled_OLS_clinics_mf1 = plm(pct_Wait_IVF_Age_1 ~ Year1Clinic	+ Aver_eff	+ pct_MHI + 
                               l1_Inequal + 
                               l1_WM_infertility +
                               l1_Koef_Marriage + pct_WM_older35 + pct_WM_older35_2 +
                               l1_pct_Internet	+ l1_pct_High_Educ,
                            data = data_clinics, index = c("X.U.FEFF.Region","Year"), model="pooling")

Pooled_OLS_clinics_mf2 = plm(pct_Wait_IVF_Age_2 ~ Year1Clinic	+ Aver_eff	+ pct_MHI + 
                               l1_Inequal + 
                               l1_WM_infertility +
                               l1_Koef_Marriage + pct_WM_older35 + pct_WM_older35_2 +
                               l1_pct_Internet	+ l1_pct_High_Educ,
                             data = data_clinics, index = c("X.U.FEFF.Region","Year"), model="pooling")

Pooled_OLS_clinics_mf3 = plm(pct_Wait_IVF_l1_WM_infertility ~ Year1Clinic	+ Aver_eff	+ pct_MHI + 
                               l1_Inequal + 
                               l1_WM_infertility +
                               l1_Koef_Marriage + pct_WM_older35 + pct_WM_older35_2 +
                               l1_pct_Internet	+ l1_pct_High_Educ, 
                             data = data_clinics, index = c("X.U.FEFF.Region","Year"), model="pooling")

Pooled_OLS_clinics_mf4 = plm(pct_Wait_IVF_l1_infertility ~ Year1Clinic	+ Aver_eff	+ pct_MHI + 
                               l1_Inequal + 
                               l1_WM_infertility +
                               l1_Koef_Marriage + pct_WM_older35 + pct_WM_older35_2 +
                               l1_pct_Internet	+ l1_pct_High_Educ, 
                            data = data_clinics, index = c("X.U.FEFF.Region","Year"), model="pooling")

#Значимы неравенство (модели 1-2), год первой клиники (модели 3-4) и бесплодие

#------------------------------------бутстрап-------------------------------------

#Есть два метода делать выборку по наблюдениям и по остаткам. Мы используем первый (correlation mod)

Pooled_OLS_clinics_mf2_lm = lm(pct_Wait_IVF_Age_2 ~ Year1Clinic	+ Aver_eff	+ pct_MHI + 
                                 l1_Inequal + 
                                 l1_WM_infertility +
                                 l1_Koef_Marriage + pct_WM_older35 + pct_WM_older35_2 +
                                 l1_pct_Internet	+ l1_pct_High_Educ, 
                             data = data_clinics)

betahat.boot <- Boot(Pooled_OLS_clinics_mf2_lm, R=25000, method="case")
summary(betahat.boot)
confint(betahat.boot)

beta_int = coefficients(Pooled_OLS_clinics_mf2)[1]
beta_year = coefficients(Pooled_OLS_clinics_mf2)[2]
beta_eff = coefficients(Pooled_OLS_clinics_mf2)[3]
beta_MHI = coefficients(Pooled_OLS_clinics_mf2)[4]
beta_Internet = coefficients(Pooled_OLS_clinics_mf2)[10]
beta_Educ = coefficients(Pooled_OLS_clinics_mf2)[11]
beta_Ineq = coefficients(Pooled_OLS_clinics_mf2)[5]
beta_Marr = coefficients(Pooled_OLS_clinics_mf2)[7]
beta_Infert = coefficients(Pooled_OLS_clinics_mf2)[6]
beta_35 = coefficients(Pooled_OLS_clinics_mf2)[8]
beta_352 = coefficients(Pooled_OLS_clinics_mf2)[9]

bias = c(mean(betahat.boot$t[,1]-beta_int), mean(betahat.boot$t[,2]-beta_year), 
         mean(betahat.boot$t[,3]-beta_eff), mean(betahat.boot$t[,4]-beta_MHI),
         mean(betahat.boot$t[,5]-beta_Ineq), mean(betahat.boot$t[,6]-beta_Infert), 
         mean(betahat.boot$t[,7]-beta_Marr), mean(betahat.boot$t[,8]-beta_35),
         mean(betahat.boot$t[,9]-beta_352), mean(betahat.boot$t[,10]-beta_Internet), 
         mean(betahat.boot$t[,11]-beta_Educ))
beta = c(beta_int, beta_year, beta_eff, beta_MHI, beta_Ineq, beta_Infert, beta_Marr, beta_35, beta_352, 
         beta_Internet, beta_Educ)
new_coef2 <- beta - bias
Bootstrap_OLS <- Pooled_OLS_clinics_mf2
Bootstrap_OLS$coefficients <- new_coef2

stargazer(Pooled_OLS_clinics_mf1, Pooled_OLS_clinics_mf3, Pooled_OLS_clinics_mf4,
          Pooled_OLS_clinics_mf2, Bootstrap_OLS,
          se=list(clse(Pooled_OLS_clinics_mf1), clse(Pooled_OLS_clinics_mf3),
                  clse(Pooled_OLS_clinics_mf4), clse(Pooled_OLS_clinics_mf2), clse(Bootstrap_OLS)), 
          title="Pooled OLS", type="text", 
          column.labels=c("Reproductive Age", "Older35", "WM Infertility", "Infertility", "Bootstrap"), 
          df=FALSE, digits=2)

#==============================================================================
#КЛАСТЕРНЫЙ АНАЛИЗ
#==============================================================================

#Обязательная нормализация данных
#Возьмём значимые факторы из модели с экологией

data_1517 <- read.csv("data_1517.csv", sep=";", dec=",", header=TRUE, encoding = "UTF-8")
names <- read.csv("names_sh.csv", sep=";", dec=",", header=TRUE, encoding = "UTF-8")
data_1517_sig <- data.frame(data_1517[2], data_1517[4], data_1517[18], data_1517[20:22], 
                            data_1517[27:28], data_1517[32], data_1517[42])

for_cluster=data.Normalization(data_1517_sig[2:10],type="n1",normalization="column")

#Метод к-средних

wss <- sapply(1:15, function(k){
  kmeans(for_cluster, k, nstart=100)$tot.withinss
})
fviz_nbclust(for_cluster, kmeans, method="wss")

gap_stat <- clusGap(for_cluster, FUN = kmeans, nstart=500, K.max=15, B=100)
fviz_gap_stat(gap_stat)

k=kmeans(for_cluster, centers=6, nstart=1000)
k
rownames(for_cluster) <- names$Region
fviz_cluster(k, data=for_cluster, labelsize = 8)
