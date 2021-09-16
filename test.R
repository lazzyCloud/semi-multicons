# import semi-multicons script, do not forget to set workspace as current repo
source('semi-multicons.R')

# set params to specify dataset path, number of must-link constraints, number of cannot-link constraints, list of k param and merge threshold MT
data_path <- "iris.csv" # dataset included in this repo
n_ml <- 30 # 30 must-link constraints
n_cl <- 60 # 60 cannot-link constraints
kList <- c(2,3,4,5,6) # list of k params
MT <- 0.5 # merge threshold, by default should be set to 0.5

# read dataset file and separate dataset with ground truth class
Dataset <- read.csv(file = data_path, header = F)
# normalize dataset
data <- as.data.frame(scale(Dataset[,c(1:(ncol(Dataset)-1))]))
class <- Dataset[,ncol(Dataset)] # last column is the ground truth class (y) for iris

# use util function in semi-multicons script to generate pairwise constraints
# otherwise you can use your own function to do that
# if you use your own function to generate constraints, be attention that ml and cl should be *matrix*, not dataframe
# where each row represent a pairwise constraint
# an example of cannot-link constraints cl
#      [,1] [,2]
# [1,]   18  103
# [2,]  136   58
# [3,]  137   89
# [4,]   20  136
# [5,]    5  127
# the first row represent that data instance no.18 has a cannot-link with no.103, where index of data instance starts from *1* (not 0! R index starts with 1)
tmp <- generate_constraints(class, n_ml, n_cl)
# if you use our util function to generate constraints, param class represents ground truth class of dataset (y)
# it should be *values* in R
# an example of class
# [1] Iris-setosa     Iris-setosa     Iris-setosa     Iris-setosa     Iris-setosa     Iris-setosa    
# [7] Iris-setosa     Iris-setosa     Iris-setosa     Iris-setosa     Iris-setosa     Iris-setosa
# n_ml and n_cl represent the number of must-link constraints and the number of cannot-link constraints
ml <- tmp[[1]]
cl <- tmp[[2]]

# use util function in semi-multicons to generate base clusterings
# you can use your own function to do that
# here we use our mpckm implementation to generate base clusterings results of mpckmeans, with k from 2 to 6
BaseClusts <- matrix(0, nrow=nrow(Dataset), ncol=length(kList))
count = 1
cnames = c()
for (K in kList) {
  BaseClusts[,count] <- mpckmeans(data,K,ml, cl)
  cnames <- cbind(cnames, paste("mpckmeans",K,sep=""))
  count = count + 1
}
# the final BaseClusts result should be a *dataframe*
# each column represents a clustering result, and must be assigned a unique name as column name
BaseClusts <- as.data.frame(BaseClusts)
names(BaseClusts) <- cnames

# call multicons function and semi-multicons function to calculate consensus clustering results
# in result, each column is a clustering result in the generated hierarchy
# X1 is the bottom level
mc_res <- multicons(BaseClusts, MT)
smc_res <- semi_multicons(BaseClusts, MT, ml, cl)