# semi-multicons
R implementation of Semi-MultiCons approach

The util functions are contained in semi-multicons.R, an example of calculating Semi-MultiCons result on Iris dataset with base clusterings MPC-Kmeans (K=[2,6]) is shown in test.R

Please ensure that you have library **hash** and **arules** installed on your local

To use the functions in semi-multicons.R, in your R script, you can source("semi-multicons.R"), as shown in test.R. Please do not forget to move semi-multicons.R in your workspace before source it

In semi-multicons.R, only four functions are frequently used

## generate_constraints(data, n_ml, n_cl)
Randomly generate pairwise constraints based on ground-truth class. Param data is the ground truth class of dataset (aka y), it should be **values** in R

An example of class
```
[1] Iris-setosa     Iris-setosa     Iris-setosa     Iris-setosa     Iris-setosa     Iris-setosa    
[7] Iris-setosa     Iris-setosa     Iris-setosa     Iris-setosa     Iris-setosa     Iris-setosa
```
n_ml and n_cl respectively represent the number of must-link constraints and cannot-link constraints

## mpckmeans(data, K, ml, cl, maxIter = 10)

Util function to implement MPC-Kmeans

Param data is the dataset, without ground-truth class (aka X), K is the number of clusters we want

ml and cl represents the must-link constraints and cannot-link constraints. ml and cl must be **matrix** type, where each row represent a pairwise constraint

An example of cannot-link constraints cl
```
       [,1] [,2]
  [1,]   18  103
  [2,]  136   58
  [3,]  137   89
  [4,]   20  136
  [5,]    5  127
```
The first row represent that data instance no.18 has a cannot-link with no.103, where index of data instance starts from **1** (not 0! R index starts with 1)

## semi_multicons(BaseClusts, MT, Must_Link, Cannot_Link)
Function to implement semi-multicons, BaseClusts represents the base clusterings calculated from our util function, or directly input from user

The BaseClusts should be a **dataframe**

An example of BaseClusts
```
     mpckmeans2 mpckmeans3 mpckmeans4 mpckmeans5 mpckmeans6
1            1          2          1          4          5
2            1          2          2          2          2
3            1          2          2          2          2
4            1          2          2          2          2
5            1          2          2          4          4
6            1          2          1          4          5
7            1          2          2          2          4
8            1          2          2          2          4
9            1          2          2          2          2
10           1          2          2          2          2
```
Each column represents a clustering result, and must be assigned a unique name as column name. Each row represents a data instance

In the example, the first row means the data instance no.1 belongs to cluster 1 for mpckmeans2, cluster 2 for mpckmeans3 etc..
The first column represent the clustering results of mpckmeans2 for data instance no.1 to no.10

MT by default set to 0.5

Must_Link/Cannot_Link represent the must-link/cannot-link constraints generated from our util function or directly input from user

The format of Must_Link and Cannot_Link should be the same as ml and cl for mpckemans function, that means

Must_Link and Cannot_Link must be **matrix** type, where each row represent a pairwise constraint

An example of cannot-link constraints Cannot_Link
```
       [,1] [,2]
  [1,]   18  103
  [2,]  136   58
  [3,]  137   89
  [4,]   20  136
  [5,]    5  127
```
The first row represent that data instance no.18 has a cannot-link with no.103, where index of data instance starts from *1* (not 0! R index starts with 1)

## multicons(BaseClusts, MT)
Function to impelemnt multicons, BaseClusts represents the base clusterings calculated from our util function, or directly input from user

The BaseClusts should be a **dataframe**

An example of BaseClusts
```
     mpckmeans2 mpckmeans3 mpckmeans4 mpckmeans5 mpckmeans6
1            1          2          1          4          5
2            1          2          2          2          2
3            1          2          2          2          2
4            1          2          2          2          2
5            1          2          2          4          4
6            1          2          1          4          5
7            1          2          2          2          4
8            1          2          2          2          4
9            1          2          2          2          2
10           1          2          2          2          2
```
Each column represents a clustering result, and must be assigned a unique name as column name. Each row represents a data instance

In the example, the first row means the data instance no.1 belongs to cluster 1 for mpckmeans2, cluster 2 for mpckmeans3 etc..

The first column represent the clustering results of mpckmeans2 for data instance no.1 to no.10

MT by default set to 0.5
## An example of calculating semi-multicons results on MPCKMeans (k=[2,6]) on Iris dataset
```
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
data <- Dataset[,c(1:(ncol(Dataset)-1))]
class <- Dataset[,ncol(Dataset)] # last column is the ground truth class (y) for iris

# use util function in semi-multicons script to generate pairwise constraints

tmp <- generate_constraints(class, n_ml, n_cl)
ml <- tmp[[1]]
cl <- tmp[[2]]

# use util function in semi-multicons to generate base clusterings
# here we use our mpckm implementation to generate base clusterings results of mpckmeans, with k from 2 to 6

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

# call multicons function and semi-multicons function to calculate consensus clustering results
# in result, each column is a clustering result in the generated hierarchy
# X1 is the bottom level

mc_res <- multicons(BaseClusts, MT)
smc_res <- semi_multicons(BaseClusts, MT, ml, cl)
```
