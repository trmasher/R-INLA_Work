# Setup
remove(list=ls())

primary_dir <- paste("D:","`Grad_Research","Rudy_RA_Work","Ohio_data",sep="\\")
setwd(paste(primary_dir,"master_data",sep="\\"))

required_packages <- c("INLA","splancs","sp","fields","maptools","lattice","abind","spdep","fastmatch")

lapply(required_packages, require, character.only=TRUE)
inla.setOption(scale.model.default=FALSE)



## Making adjacency matrix
ohio_shape <- readShapePoly("tl_2010_39_county00.shp")
ohio_nb <- poly2nb(ohio_shape)
nb2INLA("Ohio.graph",ohio_nb)
ohio_graph <- inla.read.graph("Ohio.graph")
image(inla.graph2matrix(ohio_graph),xlab="Is Neighbor?",ylab="County")



## Bym Model
voronois_ohio <- read.csv("OhioRespMort.csv")

# We need to add an ID column for the model to run
unique_names <- unique(voronois_ohio$NAME)
ref_names <- seq(1:length(unique_names))
ref_assign <- vector("numeric",length=length(voronois_ohio$NAME))
count <- 1
for (item in voronois_ohio$NAME){
  ref_assign[count] <- which(unique_names == item)
  count <- count + 1
}
voronois_ohio$ID <- ref_assign
write.csv(voronois_ohio,file = "OhioRespMort_addID.csv")

voronois_ohio <- read.csv("OhioRespMort_addID.csv")

# Finishing up model
bym_formula <- y ~ 1 + f(ID, model="bym",graph=paste(primary_dir,"master_data","Ohio.graph",sep="\\"))
y <- voronois_ohio$y
E <- voronois_ohio$E

model_bym_ohio <- inla(bym_formula, family="poisson", data= voronois_ohio, E=E, control.compute=list(dic=TRUE))
summary(model_bym_ohio)



## BYM Space-Time Model
counties <- voronois_ohio$county
spacetimebym_formula <- y ~ 1 + f(county,model="bym",graph=ohio_graph, constr=TRUE) +
  f(counties,year,model="iid", constr=TRUE) + year

model_param <- inla(spacetimebym_formula,family="poisson",data=voronois_ohio,E=E,
                  control.predictor=list(compute=TRUE),
                  control.compute=list(dic=TRUE,cpo=TRUE))

round(model_param$summary.fixed[,1:5],3)


# Plotting nonspatially or temporally structured interactions
years <- seq(1,21)
plot(years,model_param$summary.fixed[2,1]*years, type="l", main="",xlab="t",ylab=expression(beta*t), ylim=c(-0.007,0.1))
  # Pro-tip: If the above plot isn't working, try par(mar=c(1,1,1,1))' into the console a couple of times!

lines(model_param$summary.fixed[2,3]*years,lty=2)
lines(model_param$summary.fixed[2,5]*years,lty=2)
