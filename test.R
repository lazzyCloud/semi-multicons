source('semi-multicons.R')

data_path <- "../Semi-MultiCons/dataset/iris.csv"
n_ml <- 30
n_cl <- 60
kList <- c(2,3,4,5,6)
MT <- 0.5

Dataset <- read.csv(file = data_path, header = F)
data <- Dataset[,c(1:(ncol(Dataset)-1))]
class <- Dataset[,ncol(Dataset)]

tmp <- generate_constraints(class, n_ml, n_cl)
ml <- tmp[[1]]
cl <- tmp[[2]]

BaseClusts <- matrix(0, nrow=nrow(Dataset), ncol=length(kList))
count = 1
cnames = c()
for (K in kList) {
  BaseClusts[,count] <- mpckmeans(data,K,ml, cl)
  cnames <- cbind(cnames, paste("mpckmeans",K,sep=""))
  count = count + 1
}

BaseClusts <- as.data.frame(BaseClusts)
names(BaseClusts) <- cnames

mc_res <- multicons(BaseClusts, MT)
smc_res <- semi_multicons(BaseClusts, MT, ml, cl)