#======================================================================================================#
#                                  Extreme Partial Least Squares
#======================================================================================================#


#======================================================================================================#
#                                      0. Required packages  
#======================================================================================================#
library("dplyr");library("data.table");library("stats");library("base")
library("graphics");library("grDevices");library("copula");library("VineCopula")

#======================================================================================================#
#                                      I. Simulations process
#======================================================================================================#

#-------------------------------------------------------------------------------------------------------
#                                       1. Basic functions 
#-------------------------------------------------------------------------------------------------------

# Euclidean norm of a vector 
norm.vec <- function(x) sqrt(sum(x^2))

# Indicator function 
ind <- function(Yi, y)
  ifelse(Yi >= y, 1, 0)

# Normalization function
Fn <- function(X, x) {
  sum <- 0
  for (i in c(1:length(X))) {
    sum <- sum + ind(X[i], x)
  }
  return(sum / (length(X) + 1))
}
normalize<-function(variable){
  variable <- lapply(variable, function(x) {Fn(variable, x)}) %>% unlist
}

# Power link function
g<-function(x,c) x^(c) 

# Quantile function
yn<-function(Y,alpha){
  return(unname(quantile(Y, alpha)))
}

# Simulation of Pareto distribution
pareto.simu<-function(pareto.param,a,n){
  u<-runif(n,0,1) 
  y<- a*u^(-1/pareto.param)  
  return(y)}

# Kendall's tau from copula parameter (archimedean)
tau_arch_copula <- function(eta, copula){
  
  if (tolower(copula) == "copula2"){
    alpha <- eta[1]
    kappa <- eta[2]
    output <- 1 - 2 * alpha * kappa/(1 + 2 * kappa) }
  
  if (tolower(copula) != "copula2") {
    output <- tau(archmCopula(tolower(copula), param = eta, dim = 2)) }
  return(output)}

# Kendall's tau from copula parameter (Normal)
tau_normal_copula<-function(r){return(2/pi*(asin(r)))}

#-------------------------------------------------------------------------------------------------------
#                                       2. Extreme-PLS function
#-------------------------------------------------------------------------------------------------------

# Description: this function returns the beta estimate by Extreme-PLS model

# Inputs: X: matrix of covariates/ Y: vector of response variable/ y: numeric threshold

# Output: Beta estimate vector

# Extreme-PLS function
epls<-function(X,Y,y){ 
  bool<-Y>y
  X.bool<-X[bool,]
  Y.bool<-Y[bool]  
  cy<-length(Y.bool)/length(Y)    
  uy<-sum(Y.bool)/length(Y)
  if(is.matrix(X.bool)){
    mxy<-as.vector(t(X.bool)%*%Y.bool)/length(Y)
    mx <-colSums(X.bool)/length(Y)}
  if(is.vector(X.bool)){
    mxy<-as.vector(X.bool*Y.bool)/length(Y)
    mx <-X.bool/length(Y)
  }
  w<-cy*mxy-mx*uy
  return(w/norm.vec(w))
}

#-------------------------------------------------------------------------------------------------------
#                                       3. Simulation functions
#-------------------------------------------------------------------------------------------------------

#__________________________ 3.1 Simulation of X and Y function ________________________________________________________#

# Description: this function returns a list containing the vector of simulated Y and simulated X according to model M_1

# Inputs:
#n: integer, sample size
#p: integer, dimension
#beta: numeric, vector of dimension direction
#power: numeric, vector of link function power
#sigma: numeric, matrix of standard deviaton of noise
#r: numeric, signal noise
#dist: charachter, name of Y distribution  (pareto or student)
#dist.param: numeric, value or vector of Y distribution parameters
#copula.fam: integer, defining the bivariate copula family: 1 = Gaussian copula, 5=Frank copula
#copula.param: numeric, copula parameter

# Output: the first index of the list returns the simulated Y, then the other indices return the simulated X according to the link function. For example if the size of the power vector is 4,
#that is to say we will have 4 matrices X. See example below

#X.Y <- X_Y_simu(n,p,beta, power...) where power=c(3/2,1,1/2,1/4)
#X.Y[[1]] return the Y sampling according the chosen distribution Pareto or Student
#X.Y[[2]] return the matrix X1 corresponding to link function g:x->x^{3/2}
#X.Y[[3]] return the matrix X2 corresponding to link function g:x->x^1
#X.Y[[3]] return the matrix X3 corresponding to link function g:x->sqrt(x)
#X.Y[[4]] return the matrix X4 corresponding to link function g:x->x^{1/4}


# Simulation of X and Y function
X_Y_simu<-function(n,p,beta,power,sigma,r,dist,dist.param,copula.fam,copula.param){
  #Simulation of Y
  if(dist=="pareto"){Y<-pareto.simu(dist.param[1],dist.param[2],n)} #Pareto 
  if(dist=="student"){Y<-abs(rt(n,param.dist[1]))}                  #Student 
  List_X<-list()
  List_X[[1]]<-Y
  #Simulation of Epsilon conditionally on Y using Frank or Normal copulas 
  #Normalization of Y
  U1<-normalize(Y) 
  #Copula family
  if(copula.fam==0){obj <- BiCop(family = copula.fam) }              #Independence copula 0
  if(copula.fam==5 | copula.fam==1){obj <- BiCop(family = copula.fam, copula.param)} #Frank or Normal copula
  obj <- BiCop(family = copula.fam, copula.param)                   
  
  #V1 Simulation : conditionnal copula. BiCopCondSim simulates a bivariate parametric copula, where one of the variables is fixed
  V1<-matrix(0,n,p)
  for(i in 1:p){
    V1[,i]<-BiCopCondSim(n, cond.val = U1, cond.var = 1, obj) 
  }
  
  #Simulation of X from Y and Epsilon for Pareto distribution
  if(dist=="pareto"){
    sd<-vector()
    for(l in 1:length(power)){ 
      sd[l]<- g(dist.param[2]*n^{1/dist.param[1]},power[l])/(2*r)  
      List_X[[l+1]]<-g(Y,c[l])%*%t(beta)+scale(qnorm(V1,mean=0,sd=sd[l]),scale = FALSE)
    }
  } 
  
  #Simulation of X from Y and Epsilon for Student distribution
  if(dist=="student"){
    sd<-vector()
    for(l in 1:length(power)){
      sd[l]<- g(qt(1-1/n,dist.param[1]),power[l])/(2*r) 
      List_X[[l+1]]<-g(Y,power[l])%*%t(beta)+scale(qnorm(V1,mean=0,sd=sd[l]),scale = FALSE)
    }
  }
  return(List_X)
}

#__________________________ 3.2 Simulation process function ________________________________________________________#

# Description: this function returns a table containing the mean proximity criterion for each number of exceedances and link function power

# Inputs:
#n: integer, sample size
#p: integer, dimension
#beta: numeric, vector of dimension direction
#power: numeric, vector of link function power
#sigma: numeric, matrix of standard deviaton of noise
#r: numeric, signal noise
#dist: charachter, name of Y distribution  (pareto or student)
#dist.param: numeric, value or vector of Y distribution parameters
#copula.fam: integer, defining the bivariate copula family: 1 = Gaussian copula, 5=Frank copula
#copula.param: numeric, copula parameter
#N: integer, number of replications
#k.threshold<- integer, sample size coresponding to the minimum threshold


# Output: Dataframe containing the mean proximity criterion for each number of exceedances and link function power 

# Simulation process function
simu_process<-function(n,p,beta,c,sigma,r,dist,pareto.params,copula.fam,copula.param,N,k.threshold){
  #Cos: results storage table containing: simu = number of replications, k=number of exceedances, cos^2 for each link function power
  cos<-data.table(simu=numeric(),k=numeric(),cos2.x.1.5=numeric(),cos2.x=numeric(),cos2.x.0.5=numeric(),cos2.x.0.25=numeric())
  #Begining of the algorithm
  for(j in 1:N){ 
    print(paste("Replication",j))
    #Simulation of X and Y
    X.Y<-X_Y_simu(n,p,beta,c,sigma,r,dist,pareto.params,copula.fam,copula.param)
    #Estimation of beta 
    Y<-sort(X.Y[[1]]) #Simulated Y   
    for(i in k.threshold:length(Y) ){
      print(i)
      #Threshold y
      y<-Y[i]
      #Number of exceedances
      k<-(1-(i/1000))*n
      #Estimation de beta 
      w<-vector()
      cos2<-vector()
      for(l in 1:length(c)){ 
        w<-epls(X.Y[[l+1]],X.Y[[1]],y)
        cos2[l]<-(t(w)%*%beta)^2
      }
      
      cos<-as.data.table(rbind(cos,data.table(simu=j,k=k,cos2.x.1.5=cos2[1],cos2.x=cos2[2],cos2.x.0.5=cos2[3],
                                              cos2.x.0.25=cos2[4]),fill=TRUE))
    }
  }
  #Calculation of the mean proximity criterion for each number of exceedances
  cos[,m_cos2.x.1.5.V1:=mean(cos2.x.1.5),by=k]
  cos[,m_cos2.x.V1:=mean(cos2.x),by=k]
  cos[,m_cos2.x.0.5.V1:=mean(cos2.x.0.5),by=k]
  cos[,m_cos2.x.0.25.V1:=mean(cos2.x.0.25),by=k]
  cos.rep <-unique(cos[,c("k","m_cos2.x.1.5.V1","m_cos2.x.V1","m_cos2.x.0.5.V1","m_cos2.x.0.25.V1")])
  
  return(as.data.frame(cos.rep))
}

#-------------------------------------------------------------------------------------------------------
#                                      2. Simulation process
#-------------------------------------------------------------------------------------------------------

#__________________________2.1 Inputs parameters ________________________________________________________#

######### Dimension & sample size
#Two dimension
p1<-3;p2<-30
#Selected dimension : p1 or p2
p<-p1 
#Sample size
n<-1000 

######### Beta vectors
#For dimension 3
beta1<-c(1,1,0)/norm.vec(c(1,1,0))  
#For dimension 30
beta2<-c(rep(1,15),rep(0,15))/norm.vec(c(rep(1,15),rep(0,15))) 
#Selected beta
beta<-beta1
######### Link function power
c<-c(3/2,1,0.5,0.25)

######### Noise
#Standard deviation of epsilon
sigma<-0.9*diag(p)
#Noise signal ratio
r<-5 

######### Pareto and student parameters for Y distribution
#pareto.params[1] : Pareto index and also equal to 1/gamma where gamma is the tail index
#pareto.params[2] :  minimum possible value of Y 
pareto.params<-c(5,2) 

#student.param : Student degrees of freedom and also equal to 1/gamma where gamma is the tail index
student.param<-5

#Selected distribution : "pareto" or "student"
dist<-"pareto" 

#######Frank copula param for the dependence between Y and epsilon 
#copula. fam   : 0= Independent copula, 1 = Gaussian copula, 5=Frank copula
#copula. param : Frank parametrs in {0,10,20} corresponding to Kendall tau in c(0,0.67,0.82). 
#                Gaussian parameters in {0,0.856,0.96} corresponding to Kendall tau in c(0,0.67,0.82).

copula.fam<-5          
copula.param<-10

####### Number of replications
N<-100

###### The 200 largest observations 
k.threshold<-800

#__________________________2.2 Simulation Results  ________________________________________________________#


cos.rep<-simu_process(n,p,beta,c,sigma,r,dist,pareto.params,copula.fam,copula.param,N,k.threshold)

#__________________________2.3 Graphical representation  ________________________________________________________#

plot(cos.rep[,1],cos.rep[,2],type='l',ylim=c(0,1),xlab="",ylab="",col="red",
     font.lab=4, font.sub=4,font.lab=4, font.sub=4,lty=1,lwd=3,font.axis=2)
#Add horizontal grid  
axis(2, at = seq(0,1,0.2), tck = 1, lty = "dotted", lwd = par("lwd"), col = "grey", labels = NA)
#Add vertical grid
axis(1, at = seq(0,200,25), tck = 1, lty = "dotted", lwd = par("lwd"), col = "grey", labels = NA)
#Add box around plot
box()
lines(cos.rep[,1],cos.rep[,3],type='l',xlab="k*n",ylab="cos2",col="forestgreen",
      font.lab=4, font.sub=4,font.lab=4, font.sub=4,lty=1,lwd=3,font.axis=2)

lines(cos.rep[,1],cos.rep[,4],type='l',xlab="k*n",ylab="cos2",col="orange",
      font.lab=4, font.sub=4,font.lab=4, font.sub=4,lty=1,lwd=3,font.axis=2)
lines(cos.rep[,1],cos.rep[,5],type='l',xlab="k*n",ylab="cos2",col="black",
      font.lab=4, font.sub=4,font.lab=4, font.sub=4,lty=1,lwd=3,font.axis=2)



#======================================================================================================#
#                                      II. Application on real data: FADN
#======================================================================================================#

#-------------------------------------------------------------------------------------------------------
#                                      1. Data processing 
#-------------------------------------------------------------------------------------------------------

#_____________________________________ 1.1 Data loading _________________________________________________________#

List_data <- NULL   
List_temp<-NULL
List_precip<-NULL
List_sun<-NULL
#Data RICA from 2000 to 2016
i<-2000
for (file in list.files("../Data/Data Rica",full.names = TRUE)){
  List_data[[i]] <- fread(input = file,dec = ".",sep = ";",header = T) 
  i <- i+1 
}
List_data[[2014]]
#Data temperature (C) from 2014 to 2016   
i<-2014
for (file in list.files("../Data/Temperature",full.names = TRUE)){
  List_temp[[i]] <- fread(input = file) 
  i <- i+1 
}
#Data precipitation (mm) from 2014 to 2016
i<-2014
for (file in list.files("../Data/Precipitation",full.names = TRUE)){
  List_precip[[i]] <- fread(input = file) 
  i <- i+1 
}
#Data sunshine (number of hours) from 2014 to 2016
i<-2014
for (file in list.files("../Data/Sunshine",full.names = TRUE)){
  List_sun[[i]] <- fread(input = file) 
  i <- i+1 
}

#_____________________________________ 1.2 Data extraction functions_________________________________________________________#

#Extraction of cereal data 
extract.data.cereals<-function(otex,crops,data){
  data<-as.data.frame(data)
  #Extraction
  pesticides <- data$CHRPH[data$otexe==otex]         #Pesticides
  fertilizer <- data$CHREN[data$otexe==otex]         #Fertilizer
  index1<-paste(c("PRQ3",crops),collapse = "")
  yield <- data[,index1][data$otexe==otex]           #Yield
  index2<-paste(c("PBQ3",crops),collapse = "") 
  GrossProductQuintal <- data[,index2][data$otexe==otex]
  index3<-paste(c("PBV3",crops),collapse = "") 
  GrossProductEuro <- data[,index3][data$otexe==otex] 
  region <- data$REGIO[data$otexe==otex]             #Region  
  altitude <- data$ZALTI[data$otexe==otex]           #Altitude     
  index4<-paste(c("SUT3",crops),collapse = "")
  my.area <- data[,index4][data$otexe==otex]         #Cultivated area  
  retail <- data$VDETAIL[data$otexe==otex]           #Retail sale
  assurance<- data$ASSRE[data$otexe==otex]           #Insurance premiums
  indem<- data$INDAS[data$otexe==otex]               #Insurance claims
  subventions.exp<-data$SUBEX[data$otexe==otex]      #Farming subsidies
  charge.ferti<-data$CHREN[data$otexe==otex]         #Real fertiliser charge: Opening inventory + Purchases - Closing inventory
  charge.sem<-data$CHRSE[data$otexe==otex]           #Real cost of seeds and seedlings
  travaux.sc<-data$TCULT[data$otexe==otex]           #Purchases of works and services for crops
  assurance.autres<-data$ASSAU[data$otexe==otex]     #Other insurance premiums
  taxes.exp<-data$TXPRO[data$otexe==otex]            #Taxes on farm products  
  charges.soc<-data$CHSOX[data$otexe==otex]          #Farmer's personal social charges
  charges.autres<- data$ACHEX[data$otexe==otex]      #Other exceptional expenses
  
  
  #Remove rows with GrossProductQuintal==0
  yield <- yield[GrossProductQuintal != 0.00]
  region <- region[GrossProductQuintal != 0.00]
  altitude <- altitude[GrossProductQuintal != 0.00]
  GrossProductEuro <- GrossProductEuro[GrossProductQuintal != 0.00]
  my.area <- my.area[GrossProductQuintal != 0.00]
  retail <- retail[GrossProductQuintal != 0.00]
  assurance<- assurance[GrossProductQuintal != 0.00]
  indem<- indem[GrossProductQuintal != 0.00]
  fertilizer <- fertilizer[GrossProductQuintal != 0.00]
  pesticides <- pesticides[GrossProductQuintal != 0.00]
  subventions.exp<-subventions.exp[GrossProductQuintal != 0.00]     
  charge.ferti<-charge.ferti[GrossProductQuintal != 0.00]        
  charge.sem<-charge.sem[GrossProductQuintal != 0.00]         
  travaux.sc<-travaux.sc[GrossProductQuintal != 0.00]         
  assurance.autres<-assurance.autres[GrossProductQuintal != 0.00]    
  taxes.exp<-taxes.exp[GrossProductQuintal != 0.00]           
  charges.soc<-charges.soc[GrossProductQuintal != 0.00]         
  charges.autres<- charges.autres[GrossProductQuintal != 0.00]      
  GrossProductQuintal <- GrossProductQuintal[GrossProductQuintal != 0.00]
  
  #Price crop
  price = GrossProductEuro / GrossProductQuintal #gross product of one quintal
  #Cultivated area 
  for (i in c(1:10)) {
    my.area[which(my.area == i)] <- (i - 1) * 5 + 2.5
  }
  for (i in c(11:25)) {
    my.area[which(my.area == i)] <- (2 * i - 11) * 5
  }
  for (i in c(26:30)) {
    my.area[which(my.area == i)] <- (2 * i - 43) * 5
  }
  #Yield/Cultivated area 
  yield_adj <- yield 
  yield_adj[which(my.area != 0)] <- yield[which(my.area != 0)] / my.area[which(my.area != 0)]           #Price per hectar
  #Insurance premiums/Cultivated area 
  assurance.adj<- assurance
  assurance.adj[which(my.area != 0)] <-assurance[which(my.area != 0)] / my.area[which(my.area != 0)]  #Pesticide per hectar 
  #Insurance indemnisation/Cultivated area 
  indem.adj<- indem
  indem.adj[which(my.area != 0)] <-indem[which(my.area != 0)] / my.area[which(my.area != 0)]  #Pesticide per hectar 
  #Pesticides/Cultivated area 
  pesticides_adj<-pesticides
  pesticides_adj[which(my.area != 0)] <-pesticides[which(my.area != 0)] / my.area[which(my.area != 0)]  #Pesticide per hectar 
  pesticides_adj_log <-log(pesticides_adj)
  #Fertilizer/Cultivated area
  fertilizer_adj<-fertilizer
  fertilizer_adj[which(my.area != 0)] <-fertilizer[which(my.area != 0)] / my.area[which(my.area != 0)]
  fertilizer_adj_log<-log(fertilizer_adj)
  #Subsidies/Cultivated area
  subventions.exp.adj<-subventions.exp
  subventions.exp.adj[which(my.area != 0)]<-subventions.exp[which(my.area != 0)]/my.area[which(my.area != 0)]
  #Fertilizer cost/Cultivated area
  charge.ferti.adj<- charge.ferti
  charge.ferti.adj[which(my.area != 0)]<-charge.ferti[which(my.area != 0)]/my.area[which(my.area != 0)]
  #Seeds and seedlings cost/Cultivated area
  charge.sem.adj<-  charge.sem
  charge.sem.adj[which(my.area != 0)]<-charge.sem[which(my.area != 0)]/my.area[which(my.area != 0)]
  #Purchases of works and services/Cultivated area
  travaux.sc.adj<-travaux.sc     
  travaux.sc.adj[which(my.area != 0)]<-travaux.sc[which(my.area != 0)]/my.area[which(my.area != 0)]
  #Other insurance premiums/Cultivated area
  assurance.autres.adj<-assurance.autres
  assurance.autres.adj[which(my.area != 0)]<-assurance.autres[which(my.area != 0)]/my.area[which(my.area != 0)]
  #Taxes/Cultivated area
  taxes.exp.adj<- taxes.exp
  taxes.exp.adj[which(my.area != 0)]<-taxes.exp[which(my.area != 0)]/my.area[which(my.area != 0)]
  #social charges/Cultivated area
  charges.soc.adj<-charges.soc 
  charges.soc.adj[which(my.area != 0)]<-charges.soc[which(my.area != 0)]/my.area[which(my.area != 0)]
  #Other expenses/Cultivated area
  charges.autres.adj<- charges.autres   
  charges.autres.adj[which(my.area != 0)]<-charges.autres[which(my.area != 0)]/my.area[which(my.area != 0)]
  
  return(data.frame(pesticides,pesticides_adj,pesticides_adj_log,fertilizer,fertilizer_adj,
                    fertilizer_adj_log,my.area,price,yield_adj,yield,retail,region,altitude,
                    assurance,assurance.adj,indem,indem.adj,subventions.exp,subventions.exp.adj,
                    charge.ferti, charge.ferti.adj, charge.sem, charge.sem.adj,
                    travaux.sc,travaux.sc.adj,  assurance.autres, assurance.autres.adj,
                    taxes.exp, taxes.exp.adj, charges.soc, charges.soc.adj,
                    charges.autres,charges.autres.adj))
}
#Extraction of wheat data and weather variable fusion
process.data.wheat.meteo<-function(data,temp,precip,sun){
  data.wheat<-extract.data.cereals(15,"BLET",data)
  #Merge with temperature 
  data.wheat<-merge(data.wheat,temp,by.x=c("region","altitude"),by.y=c("REGIO","ZALTI"),all.x=T)
  #Merge with precipitation
  data.wheat<-merge(data.wheat,precip[,c("REGIO","ZALTI","P","PD")],by.x=c("region","altitude"),by.y=c("REGIO","ZALTI"),all.x=T)
  #Merge with sunshine
  data.wheat<-merge(data.wheat,sun[,c("REGIO","ZALTI","S","SD")],by.x=c("region","altitude"),by.y=c("REGIO","ZALTI"),all.x=T)
  data.wheat<-as.data.table(data.wheat)
  ##Outliers processing : removal according to the agricultural expert
  data.wheat<-data.wheat[!(price>30)]
  data.wheat<-data.wheat[(price>=0)]
  data.wheat<-data.wheat[!(yield_adj>200)] 
  data.wheat<-data.wheat[!(yield_adj==0)] 
  data.wheat<-data.wheat[fertilizer_adj>0] 
  data.wheat<-data.wheat[pesticides_adj>0] 
  return(as.data.frame(data.wheat))
}

#_____________________________________ 1.3 Database _________________________________________________________#

#Extraction
List_data_wheat<-NULL  
year_i<-2014
year_f<-2016
for(i in year_i:year_f){
  List_data_wheat[[i]]<-process.data.wheat.meteo(List_data[[i]],List_temp[[i]],List_precip[[i]],List_sun[[i]])
  List_data_wheat[[i]]<-as.data.table(List_data_wheat[[i]])
  List_data_wheat[[i]][,year:=i]
}
#Selected year
i=2014
data<-List_data_wheat[[year_i]]
#Selected variables
var_keep <- c("pesticides_adj","fertilizer_adj","price","yield_adj","assurance.adj","indem.adj",
              "subventions.exp.adj", "charge.sem.adj", 
              "travaux.sc.adj", "assurance.autres.adj", "taxes.exp.adj", "charges.soc.adj","T") 
data.selec <- data[,var_keep,with=F]

#-------------------------------------------------------------------------------------------------------
#                                      2. Assumption testing
#           visual checks of whether the heavy-tailed assumption makes sense for wheat yield data 
#-------------------------------------------------------------------------------------------------------

#_____________________________________ 2.1 Descriptive statistics _________________________________________________________#

#### Inverse yield
Y<-data.selec$yield_adj
Y<-sort(1/Y) 
summary(Y)

#_____________________________________ 2.2 Inverse yield histogram _________________________________________________________#

plot(c(0.005,0.067),c(0,150),type = "n",xlab="",ylab="",font.lab=4, font.sub=4,font.axis=2)
#Add horizontal grid  
axis(2, at = c(0,25,50,75,100,125,150), tck = 1, lty = "dotted", lwd = par("lwd"), col = "grey", labels = NA)
#Add vertical grid
axis(1, at = c(0,0.01,0.02,0.03,0.04,0.05,0.06,0.07), tck = 1, lty = "dotted", lwd = par("lwd"), col = "grey", labels = NA)
#Add box around plot
box()
hist(Y, col="steelblue4", border="black", prob = TRUE, xlab = " ",ylab=" ",
     main = " ",breaks=seq(from=min(Y), to=max(Y),length.out=50),add = TRUE)

#_____________________________________ 2.2 Inverse yield Quantile-Quantile plot _________________________________________________________#

#### Quantile-Quantile plot function 
qqplt<-function(Y,k,n){
  ll<-cbind(log(k/1),log(Y[n])-log(Y[n-k]))
  plot(log(k/1),log(Y[n])-log(Y[n-k]),xlim=c(0,log(k)),pch=4,xlab="",ylim=c(0,1.5),ylab="", #ylim=c(0,2.5),
       font.lab=4, font.sub=4,font.axis=2)
  #grid(8, lty = "dotted", lwd = par("lwd"),col="grey")
  #Add horizontal grid  
  axis(2, at = c(0,0.25,0.5,0.75,1,1.25,1.5), tck = 1, lty = "dotted", lwd = par("lwd"), col = "grey", labels = NA)
  #Add vertical grid
  axis(1, at = seq(0,5,1), tck = 1, lty = "dotted", lwd = par("lwd"), col = "grey", labels = NA)
  #Add box around plot
  box()
  for (i in 2:k) {
    li<-cbind(log(k/i),log(Y[n-i+1])-log(Y[n-k]))
    points(log(k/i),log(Y[n-i])-log(Y[n-k]),pch=4)
    ll<-rbind(ll,li)
  }
  abline(coefficients(lm(ll[,2]~ll[,1])),col="red")
}

#### Quantile-Quantile plot
n<-length(Y)
k<-150
qqplt(Y,k,length(Y))

#_____________________________________ 2.3 Inverse yield Hill plot  _________________________________________________________#

#### Hill estimator
logY=rev(log(Y))
n=length(Y)
hi=1/(1:n)*cumsum(logY)-logY
hi.IC=1.96/sqrt(1:n)*hi

#### Hill plot

plot(1:n,hi,type="l",ylim=range(c(hi+hi.IC,hi-hi.IC)),xlab="",xlim=c(0,150),ylab="",
     font.sub=4,font.lab=4,font.axis=2)
#grid(8, lty = "dotted", lwd = par("lwd"),col="grey")
#Add horizontal grid  
axis(2, at = c(0,0.2,0.4,0.6,0.8,1), tck = 1, lty = "dotted", lwd = par("lwd"), col = "grey", labels = NA)
#Add vertical grid
axis(1, at = c(0,25,50,75,100,125,150), tck = 1, lty = "dotted", lwd = par("lwd"), col = "grey", labels = NA)
#Add box around plot
box()
lines(1:n,hi+hi.IC,col="blue",lwd=1.5,lty=3)
lines(1:n,hi-hi.IC,col="blue",lwd=1.5,lty=3)
lines(1:n,hi,lwd=1.5)

#-------------------------------------------------------------------------------------------------------
#                                      3. Application of Extreme-PLS model
#            
#-------------------------------------------------------------------------------------------------------

#_____________________________________ 3.1 cor(w X,X|Y>y)  _________________________________________________________#
#Conditional correlation between the projected covariate and each variable 

#### Input
#Covariates matrix X
X<-as.data.frame(data.selec) 
#Observations number
n<-dim(X)[1]
#Number of exceedences
nb_obs<-150
dep<-n-nb_obs
#Inverse yield variable
X$Y<-1/(X$yield_adj)
#Deletion of yield variable
X<-X[,-4]
#Sort data
X<-X[order(X$Y),]
#Inverse yield
Y<-X$Y  
#Deletion of Inverse yield variable from X
X<-X[,-13]

#Results storage table
res<- data.table(coeff=character(),sample_size=integer(),pesticides_adj_log=numeric(),
                 fertilizer_adj_log=numeric(),price=numeric(),
                 assurance=numeric(),indem=numeric(), subventions.exp=numeric() ,
                 charge.sem=numeric(),
                 travaux.sc=numeric(),
                 assurance.autres=numeric(),
                 taxes.exp=numeric(),charges.soc=numeric(),temp=numeric())

for(i in dep:n){ 
  #Threshold
  y<-Y[i] 
  #Number of exceedances
  k<-(1-(i/n))*n
  #Weight estimate using Extreme-PLS model
  w<-epls(as.matrix(X),Y,y)
  #Correlation between w X and X given Y>y
  corr<-vector("numeric")
  for(j in 1:dim(X)[2]){
    corr[j]<-cov(as.matrix(X[Y>y,])%*%w,as.matrix(X[Y>y,j]))/(sd(as.matrix(X[Y>y,])%*%w)*sd(as.matrix(X[Y>y,j])))
  }
  
  
  res<-rbind(res,data.table(coeff="w",sample_size=k,pesticides_adj_log=w[1],
                            fertilizer_adj_log= w[2],price=w[3],
                            assurance=w[4],indem=w[5],subventions.exp=w[6],
                            charge.sem=w[7],
                            travaux.sc=w[8],
                            assurance.autres=w[9],taxes.exp=w[10],charges.soc=w[11],
                            temp=w[12]))
  res<-rbind(res,data.table(coeff="corr",sample_size=k,pesticides_adj_log=corr[1],
                            fertilizer_adj_log= corr[2],price=corr[3],
                            assurance=corr[4],indem=corr[5],subventions.exp=corr[6],
                            charge.sem=corr[7],
                            travaux.sc=corr[8],
                            assurance.autres=corr[9],taxes.exp=corr[10],charges.soc=corr[11],
                            temp=corr[12]))
}

par(mar=c(5.1, 2.1, 4.1, 10.1), xpd=FALSE)
plot(res[coeff=="corr" & sample_size>= 50,sample_size],abs(res[coeff=="corr" & sample_size>= 50,pesticides_adj_log]),type='l',
     ylim=c(0,1),xlab="",ylab="",col="brown",font.lab=4,
     font.sub=4,font.lab=4, font.sub=4,lty=1,lwd=2,font.axis=2, xaxp = c(50, 150, 5))
#Add horizontal grid  
axis(2, at = c(0,0.2,0.4,0.6,0.8,1), tck = 1, lty = "dotted", lwd = par("lwd"), col = "grey", labels = NA)
#Add vertical grid
axis(1, at = seq(50,150,20), tck = 1, lty = "dotted", lwd = par("lwd"), col = "grey",labels = NA)
#Add box around plot
box()
lines(res[coeff=="corr" & sample_size>= 50,sample_size],abs(res[coeff=="corr" & sample_size>= 50,fertilizer_adj_log]),type='l',col="green",font.sub=4,lty=1,lwd=2)
lines(res[coeff=="corr"& sample_size>= 50,sample_size],abs(res[coeff=="corr"& sample_size>= 50,price]),type='l',col="purple",font.sub=4,lty=1,lwd=2)
lines(res[coeff=="corr"& sample_size>= 50,sample_size],abs(res[coeff=="corr"& sample_size>= 50,assurance]),type='l',col="yellow",font.sub=4,lty=1,lwd=2)
lines(res[coeff=="corr"& sample_size>= 50,sample_size],abs(res[coeff=="corr"& sample_size>= 50,indem]),type='l',col="orange",font.sub=4,lty=1,lwd=2)
lines(res[coeff=="corr"& sample_size>= 50,sample_size],abs(res[coeff=="corr"& sample_size>= 50,subventions.exp]),type='l',col="darkblue",font.sub=4,lty=1,lwd=2)
lines(res[coeff=="corr"& sample_size>= 50,sample_size],abs(res[coeff=="corr"& sample_size>= 50,charge.sem]),type='l',col="deeppink",font.sub=4,lty=1,lwd=2)
lines(res[coeff=="corr"& sample_size>= 50,sample_size],abs(res[coeff=="corr"& sample_size>= 50,travaux.sc]),type='l',col="deepskyblue",font.sub=4,lty=1,lwd=2)
lines(res[coeff=="corr"& sample_size>= 50,sample_size],abs(res[coeff=="corr"& sample_size>= 50,assurance.autres]),type='l',col="darkgreen",font.sub=4,lty=1,lwd=2)
lines(res[coeff=="corr"& sample_size>= 50,sample_size],abs(res[coeff=="corr"& sample_size>= 50,taxes.exp]),type='l',col="red",font.sub=4,lty=1,lwd=2)
lines(res[coeff=="corr"& sample_size>= 50,sample_size],abs(res[coeff=="corr"& sample_size>= 50,charges.soc]),type='l',col="deeppink4",font.sub=4,lty=1,lwd=2)
lines(res[coeff=="corr"& sample_size>= 50,sample_size],abs(res[coeff=="corr"& sample_size>= 50,temp]),type='l',col="steelblue4",font.sub=4,lty=1,lwd=2)
legend("topright", inset=c(-0.45,0),xpd=TRUE,legend = c("Pesticides","Fertilizers","Selling prices",
                                                        "Crop insurance purchase","Insurance claims","Farm subsidies","Seeds & seedlings costs",
                                                        "Works & services purchase","Other insurance premiums", "Income taxes",
                                                        "Personal social security costs","Temperature"),
       col = c("brown","green","purple","yellow","orange","darkblue","deeppink","deepskyblue",
               "darkgreen","red","deeppink4","steelblue4"),lty = 1,lwd=2,cex = 0.7,bty = "n")

#_____________________________________ 3.2 cor(w X,Y|Y>y)  _________________________________________________________#
#     Description : Conditional correlation between the projected covariate w X and Y on all directions

##### Conditional correlation function
correlation<-function(X,Y,y){
  corXY<-vector("numeric") 
  for(j in 1:ncol(X)){
    corXY[j]<-cor(Y[Y>y], X[Y>y,j],method = c("pearson"))
  }
  return(corXY)
}

##### Conditional correlation
#Covariates matrix X
X<-as.data.frame(data.selec) 
#Observations number
n<-dim(X)[1]
#Number of exceedences
nb_obs<-150
dep<-n-nb_obs
#Inverse yield variable
X$Y<-1/(X$yield_adj)
#Deletion of yield variable
X<-X[,-4]
#Sort data
X<-X[order(X$Y),]
#Inverse yield
Y<-X$Y  
#Deletion of Inverse yield variable from X
X<-X[,-13]
#Results storage table
res.residu<- data.table(residu=integer(),sample_size=integer(),pesticides_adj_log=numeric(),fertilizer_adj_log=numeric(),
                        price=numeric(),assurance=numeric(),indem=numeric(),subventions.exp=numeric(),
                        charge.sem=numeric(),travaux.sc=numeric(),assurance.autres=numeric(),
                        taxes.exp=numeric(),charges.soc=numeric(),temp=numeric(),corr=numeric(),corr2=numeric(),
                        cor.pest.Y=numeric(),cor.fert.Y=numeric(),cor.price.Y=numeric(),cor.ass.Y=numeric(),
                        cor.indem.Y=numeric(),cor.sub.Y=numeric(),cor.chargesem.Y=numeric(),cor.chargessc.Y=numeric(),
                        cor.assautres.Y=numeric(),cor.taxesexp.Y=numeric(),cor.chargessoc.Y=numeric(),cor.temp.Y=numeric()
)
#Conditional correlation for the first direction l=1
for(i in dep:n){ 
  #Threshold
  y<-Y[i] 
  #Number of exceedances
  k<-(1-(i/n))*n
  #Weight estimate using Extreme-PLS model  
  w<-epls(as.matrix(X),Y,y)
  #Conditional correlation between w X and Y given Y>y
  cor1<-cov(as.matrix(X[Y>y,])%*%w,Y[Y>y])/(sd(as.matrix(X[Y>y,])%*%w)*sd(Y[Y>y]))
  cor2<-cor(as.matrix(X[Y>y,])%*%w, Y[Y>y], method = c("pearson"))
  #Conditional correlation between X and Y given Y>y
  corXY<-correlation(X,Y,y)
  res.residu<-rbind(res.residu,data.table(residu=0,sample_size=k,pesticides_adj_log=w[1],
                                          fertilizer_adj_log= w[2],price=w[3],assurance=w[4],
                                          indem=w[5],subventions.exp=w[6] ,charge.sem=w[7],travaux.sc=w[8],
                                          assurance.autres=w[9],taxes.exp=w[10],charges.soc=w[11],temp=w[12],
                                          corr=cor1[1,1],corr2=cor2[1,1],cor.pest.Y=corXY[1],
                                          cor.fert.Y=corXY[2],cor.price.Y=corXY[3],cor.ass.Y=corXY[4],cor.indem.Y=corXY[5],
                                          cor.sub.Y=corXY[6],cor.chargesem.Y=corXY[7], 
                                          cor.chargessc.Y=corXY[8],cor.assautres.Y=corXY[9],
                                          cor.taxesexp.Y=corXY[10],cor.chargessoc.Y=corXY[11],
                                          cor.temp.Y=corXY[12]))
  
}
#Directions number -1
K<-11
#Initialization
R0<-X
#Conditional correlation for other directions l   {2,...,11}
for(j in 1:K){
  #Table of conditional correlation cor(w X,Y|Y>y) results 
  filt<-res.residu[residu==j-1,c("sample_size","corr")]
  #Classe 1 : k ∈  {0,...,50}  and Classe 2 :   k ∈  {50,...,150}
  filt[sample_size<50,class:=1]
  filt[sample_size>=50,class:=2]
  #Maximum conditional correlation per class
  filt[,max:=max(corr,na.rm = TRUE),class]
  #Maximum conditional correlation in class 2
  samp.max<-filt[class==2 & corr==max ,]$sample_size #0.5876041 corresponding to number of exceedances k = 97
  #Corresponding threshold 
  y<-unname(quantile(Y,1-(samp.max/length(Y))))
  #Weight 
  w.res1<-epls(as.matrix(R0),Y,y)
  #w.res1*R0
  Rbeta<-vector("numeric")
  for(i in 1:dim(R0)[1]){
    Rbeta[i]<- scalar(R0[i,],w.res1)}
  #Residu
  R<- R0-Rbeta%*%t(w.res1)
  #PCA of R to project it into perpendicular space at the first direction
  pca<-prcomp(R)
  Rpca<-pca$x 
  for(i in dep:n){   
    #Threshold
    y<-Y[i] 
    #Number of exceedances
    k<-(1-(i/n))*n
    #Weight estimate using Extreme-PLS model  
    w<-epls(as.matrix(Rpca),Y,y)
    #Conditional correlation between w R and Y given Y>y
    cor1<-cov(Rpca[Y>y,]%*%w,Y[Y>y])/(sd(Rpca[Y>y,]%*%w)*sd(Y[Y>y]))
    cor2<-cor(Rpca[Y>y,]%*%w, Y[Y>y], method = c("pearson"))
    #Conditional correlation between R and Y given Y>y
    corXY<-correlation(Rpca,Y,y)
    res.residu<-rbind(res.residu,data.table(residu=j,sample_size=k,pesticides_adj_log=w[1],
                                            fertilizer_adj_log= w[2],price=w[3],assurance=w[4],
                                            indem=w[5],subventions.exp=w[6] ,charge.sem=w[7],travaux.sc=w[8],
                                            assurance.autres=w[9],taxes.exp=w[10],charges.soc=w[11],temp=w[12],
                                            corr=cor1[1,1],corr2=cor2[1,1],cor.pest.Y=corXY[1],
                                            cor.fert.Y=corXY[2],cor.price.Y=corXY[3],cor.ass.Y=corXY[4],cor.indem.Y=corXY[5],
                                            cor.sub.Y=corXY[6],cor.chargesem.Y=corXY[7], 
                                            cor.chargessc.Y=corXY[8],cor.assautres.Y=corXY[9],
                                            cor.taxesexp.Y=corXY[10],cor.chargessoc.Y=corXY[11],
                                            cor.temp.Y=corXY[12]))
    
  }
  R0<-Rpca
}

#### Graphical representation of conditional correlation
par(mar=c(5.1, 2.1, 4.1, 10.1), xpd=FALSE)
plot(res.residu[residu==0 &sample_size>=50,sample_size],abs(res.residu[residu==0&sample_size>=50,corr]),type='l',
     ylim=c(0,1),xlab="",ylab="",col="blue",font.lab=4,
     font.sub=4,font.lab=4, font.sub=4,lty=1,lwd=2,font.axis=2,main="", xaxp = c(50, 150, 5))
#Add horizontal grid  
axis(2, at = c(0,0.2,0.4,0.6,0.8,1), tck = 1, lty = "dotted", lwd = par("lwd"), col = "grey", labels = NA)
#Add vertical grid
axis(1, at = seq(50,150,20), tck = 1, lty = "dotted", lwd = par("lwd"), col = "grey", labels = NA)
#Add box around plot
box()
lines(res.residu[residu==1&sample_size>=50,sample_size],abs(res.residu[residu==1&sample_size>=50,corr]),type='l',
      ylim=c(0,1),xlim=c(0,150),xlab="",ylab="",col="red",font.lab=4,
      font.sub=4,font.lab=4, font.sub=4,lty=1,lwd=2,font.axis=2,main="")
lines(res.residu[residu==2&sample_size>=50,sample_size],abs(res.residu[residu==2&sample_size>=50,corr]),type='l',
      ylim=c(0,1),xlim=c(0,150),xlab="",ylab="",col="darkgreen",font.lab=4,
      font.sub=4,font.lab=4, font.sub=4,lty=1,lwd=2,font.axis=2,main="")
lines(res.residu[residu==3&sample_size>=50,sample_size],abs(res.residu[residu==3&sample_size>=50,corr]),type='l',
      ylim=c(0,1),xlim=c(0,150),xlab="",ylab="",col="orange",font.lab=4,
      font.sub=4,font.lab=4, font.sub=4,lty=1,lwd=2,font.axis=2,main="")
lines(res.residu[residu==4&sample_size>=50,sample_size],abs(res.residu[residu==4&sample_size>=50,corr]),type='l',
      ylim=c(0,1),xlim=c(0,150),xlab="",ylab="",col="purple",font.lab=4,
      font.sub=4,font.lab=4, font.sub=4,lty=1,lwd=2,font.axis=2,main="")
lines(res.residu[residu==5&sample_size>=50,sample_size],abs(res.residu[residu==5&sample_size>=50,corr]),type='l',
      ylim=c(0,1),xlim=c(0,150),xlab="",ylab="",col="black",font.lab=4,
      font.sub=4,font.lab=4, font.sub=4,lty=1,lwd=2,font.axis=2,main="")
lines(res.residu[residu==6&sample_size>=50,sample_size],abs(res.residu[residu==6&sample_size>=50,corr]),type='l',
      ylim=c(0,1),xlim=c(0,150),xlab="",ylab="",col="deeppink",font.lab=4,
      font.sub=4,font.lab=4, font.sub=4,lty=1,lwd=2,font.axis=2,main="")
lines(res.residu[residu==7&sample_size>=50,sample_size],abs(res.residu[residu==7&sample_size>=50,corr]),type='l',
      ylim=c(0,1),xlim=c(0,150),xlab="",ylab="",col="yellow",font.lab=4,
      font.sub=4,font.lab=4, font.sub=4,lty=1,lwd=2,font.axis=2,main="")
lines(res.residu[residu==8&sample_size>=50,sample_size],abs(res.residu[residu==8&sample_size>=50,corr]),type='l',
      ylim=c(0,1),xlim=c(0,150),xlab="",ylab="",col="darkblue",font.lab=4,
      font.sub=4,font.lab=4, font.sub=4,lty=1,lwd=2,font.axis=2,main="")
lines(res.residu[residu==9&sample_size>=50,sample_size],abs(res.residu[residu==9&sample_size>=50,corr]),type='l',
      ylim=c(0,1),xlim=c(0,150),xlab="",ylab="",col="green",font.lab=4,
      font.sub=4,font.lab=4, font.sub=4,lty=1,lwd=2,font.axis=2,main="")
lines(res.residu[residu==10&sample_size>=50,sample_size],abs(res.residu[residu==10&sample_size>=50,corr]),type='l',
      ylim=c(0,1),xlim=c(0,150),xlab="",ylab="",col="deepskyblue",font.lab=4,
      font.sub=4,font.lab=4, font.sub=4,lty=1,lwd=2,font.axis=2,main="")
lines(res.residu[residu==11&sample_size>=50,sample_size],abs(res.residu[residu==11&sample_size>=50,corr]),type='l',
      ylim=c(0,1),xlim=c(0,150),xlab="",ylab="",col="steelblue4",font.lab=4,
      font.sub=4,font.lab=4, font.sub=4,lty=1,lwd=2,font.axis=2,main="")

legend("topright", inset=c(-0.17,0),xpd=TRUE, legend=c("\u2113=1","\u2113=2","\u2113=3","\u2113=4",
                                                       "\u2113=5","\u2113=6","\u2113=7","\u2113=8",
                                                       "\u2113=9","\u2113=10","\u2113=11","\u2113=12"),
       col = c("blue","red","darkgreen","orange","purple","black","deeppink","yellow","darkblue",
               "green","deepskyblue","steelblue4"),
       lty = 1,lwd=2,cex = 0.7,bty = "n")

###### Graphical representation of maximum conditionnal correlation for each direction
#Table of conditional correlation cor(w X,Y|Y>y) results 
filt<-res.residu[,c("residu","sample_size","corr")]
#Classe 1 : k ∈  {0,...,50}  and Classe 2 :   k ∈  {50,...,150}
filt[sample_size<=50,class:=1]
filt[sample_size>50,class:=2]
#Maximum conditional correlation per class
filt[,max:=max(corr,na.rm = TRUE),c("residu","class")]
#Maximum conditional correlation in class 2 for each direction
filt[class==2 & corr==max ,]

plot(1:12,filt[class==2 & corr==max,]$corr,type='b',
     xlab=" ",ylab=" ",font.lab=4, font.sub=4,font.axis=2,
     lty=1,lwd=2)  
#Add horizontal grid  
axis(2, at = seq(0,1,0.1), tck = 1, lty = "dotted", lwd = par("lwd"), col = "grey", labels = NA)
#Add vertical grid
axis(1, at = seq(0,11,1), tck = 1, lty = "dotted", lwd = par("lwd"), col = "grey", labels = NA)
#Add box around plot
box()
lines(1:12,filt[class==2 & corr==max,]$corr,type='b',
      xlab=" ",ylab=" ",font.lab=4, font.sub=4,font.axis=2,
      lty=1,lwd=2)  

#_____________________________________ 3.3 Conditional quantiles estimation _________________________________________________________#

#### Choice of threshold 
#Table of conditional correlation cor(w X, Y|Y>y) results 
filt<-res[,c("sample_size","corr")]
#Classe 1 : k ∈  {0,...,50}  and Classe 2 :   k ∈  {50,...,150}
filt[sample_size<=50,class:=1]
filt[sample_size>50,class:=2]
#Maximum conditional correlation per class
filt[,max:=max(corr,na.rm = TRUE),class]
#Maximum conditional correlation
filt[class==2 & corr==max ,] #0.5876041 corresponding to number of exceedances k = 97
#Selected y
y<-unname(quantile(Y,1-(97/length(Y))))

#### Weight 
beta<-epls(as.matrix(X),Y,y)
Xbeta<-vector("numeric")
for(i in 1:dim(X)[1]){
  Xbeta[i]<- scalar(X[i,],beta)
}

#### Useful functions
# Quadratic Kernel
quadratic.kernel<-function(u){
  if(u<=1 & u>=-1){return(15*(1-u^2)^2/16)}
  else{return(0)}
}
# F bar estimation
Fbar.est<-function(x,h,X,kernel.function,Y,y){
  f1 <-0
  f2<-0
  for(i in c(1:length(X))){   
    noyau2 <-  kernel.function((x-X[i])/h)
    noyau1 <-  noyau2*ind(Y[i],y)
    f1 <-f1 + noyau1
    f2 <-f2 + noyau2
  }
  f<-f1/f2
  return(f)
}

### Conditional quantile  
logY<-log(Y)
logXbeta<-log(Xbeta)
x.vec<-seq(min(logXbeta),max(logXbeta),length.out = 100)
y.vec<-seq(min(logY),max(logY),length.out = 100)
band<-0.15*(max(logXbeta)-min(logXbeta))
tab<- data.table(i=integer(),j=integer(),x=numeric(),y=numeric(),Fbar=numeric())
for(i in 1:length(x.vec)){
  for(j in 1:length(y.vec)){
    Fun<-Fbar.est(x.vec[i],band,logXbeta,quadratic.kernel,logY,y.vec[j])
    tab<-rbind(tab,data.table(i=i,j=j,x=x.vec[i],y=y.vec[j],Fbar=Fun))
  }
}

# alpha 5% 
res1<-tab
s1<-50/949 
res1<-res1[Fbar<s1,ind:=1]
res1<-res1[,maxFbar := max(Fbar),by=c("i","ind")]
res1<-res1[(maxFbar<1)& (Fbar==maxFbar),ind2 := 11]
res1<-res1[,miny := min(y),by=c("i","ind","ind2")]
res1<-res1[ind==1 & ind2==11,]
res1<-res1[,c("x","miny","maxFbar")] %>% unique
# alpha 15% 
res2<-tab
s2<-150/949 
res2<-res2[Fbar<s2,ind_2:=1]
res2<-res2[,maxFbar2 := max(Fbar),by=c("i","ind_2")]
res2<-res2[(maxFbar2<1)& (Fbar==maxFbar2),ind2_2 := 11]
res2<-res2[,miny2 := min(y),by=c("i","ind_2","ind2_2")]
res2<-res2[ind_2==1 & ind2_2==11,]
res2<-res2[,c("x","miny2","maxFbar2")] %>% unique


#### Scatter plot
plot(logXbeta, logY, xlab="  ",ylab=" ",col="black",font.lab=4,
     font.sub=4,lty=1,font.axis=2,pch='*',main="",xlim=c(6,10.2))
#grid(8, lty = "dotted", lwd = par("lwd"),col="grey")
#Add horizontal grid  
axis(2, at = seq(-5,3,0.5) , tck = 1, lty = "dotted", lwd = par("lwd"), col = "grey", labels = NA)
#Add vertical grid
axis(1, at = seq(6,10,1), tck = 1, lty = "dotted", lwd = par("lwd"), col = "grey", labels = NA)
#Add box around plot
box()
#Estimated conditional quantile 5%
lines(res1$x[-c(90:100)],res1$miny[-c(90:100)],col="red",xlab="  ",ylab=" ",lty=1,lwd=2)
#Estimated conditional quantile 15%
lines(res2$x[-c(90:100)],res2$miny2[-c(90:100)],col="blue",xlab="  ",ylab=" ",lty=1,lwd=2)

#_____________________________________ 3.4 Estimated link function using kernal estimator _________________________________________________________#

#### Link function estimation
G<-function(h,X,kernel.function,Y,y){
  f1 <-0
  f2<-0
  for(i in c(1:length(Y))){    
    noyau2 <- kernel.function((y-Y[i])/h)
    noyau1 <-  noyau2*X[i]
    f1 <-f1 + noyau1
    f2 <-f2 + noyau2
  }
  f<-f1/f2
  return(f)
}

#### inverse linear regression 
lm.model <- lm(logXbeta~logY)
#### G estimation 
band<-0.15*(max(logY)-min(logY))
y.vec<-seq(min(logY),max(logY),length.out = 200)
tab<- data.table(y=numeric(),G=numeric())
for(i in 1:length(y.vec)){
  Gest<-G(band,logXbeta,quadratic.kernel,logY,y.vec[i])
  tab<-rbind(tab,data.table(y=y.vec[i],G=Gest))
}
#### Scatter plot
plot(logY, logXbeta, xlab="  ",ylab=" ",col="black",font.lab=4,
     font.sub=4,lty=1,font.axis=2,pch='*',main="",ylim=c(6,10.2))
#grid(8, lty = "dotted", lwd = par("lwd"),col="grey")
#Add horizontal grid  
axis(2, at = seq(6,10,1), tck = 1, lty = "dotted", lwd = par("lwd"), col = "grey", labels = NA)
#Add vertical grid
axis(1, at = seq(-5,3,0.5), tck = 1, lty = "dotted", lwd = par("lwd"), col = "grey", labels = NA)
#Add box around plot
box()
#Inverse linear regression 
abline(lm.model,col="red" ,lty=1,lwd=2)
#kernel estimate of the link function
lines(tab$y,tab$G,col="blue",lty=1,lwd=2)
