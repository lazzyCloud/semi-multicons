library(hash)
library("arules")

# randomly generate pairwise constraints based on ground-truth class
# data is the ground truth class of dataset (aka y)
# it should be *values* in R
# an example of class
# [1] Iris-setosa     Iris-setosa     Iris-setosa     Iris-setosa     Iris-setosa     Iris-setosa    
# [7] Iris-setosa     Iris-setosa     Iris-setosa     Iris-setosa     Iris-setosa     Iris-setosa
# n_ml and n_cl respectively represent the number of must-link constraints and cannot-link constraints
generate_constraints <- 
  function(data, n_ml, n_cl) {
    Data_Size <- length(data)
    
    must_link = matrix(0L, nrow=n_ml, ncol=2)
    cannot_link = matrix(0L, nrow=n_cl, ncol=2)
    size_m = 0
    size_c = 0
    
    while ((size_m < n_ml) || (size_c < n_cl)) {
      i <- sample(1:Data_Size, size = 1)
      j <- sample(1:Data_Size, size = 1)
      if (i != j) {
        if (data[i] == data[j]) {
          if (size_m < n_ml) {
            size_m = size_m + 1
            must_link[size_m,1] = i
            must_link[size_m,2] = j
          }
          
        } else {
          if (size_c < n_cl) {
            size_c = size_c + 1
            cannot_link[size_c,1] = i
            cannot_link[size_c,2] = j
          }
        }
      }
    }
    return(list(must_link,cannot_link))
  }

# util function to implement MPC-Kmeans
# data is the dataset, without ground-truth class (aka X)
# K is the number of clusters we want
# ml and cl represents the must-link constraints and cannot-link constraints
# ml and cl must be *matrix* type
# where each row represent a pairwise constraint
# an example of cannot-link constraints cl
#      [,1] [,2]
# [1,]   18  103
# [2,]  136   58
# [3,]  137   89
# [4,]   20  136
# [5,]    5  127
# the first row represent that data instance no.18 has a cannot-link with no.103, where index of data instance starts from *1* (not 0! R index starts with 1)
mpckmeans <-
  function(data, K, ml, cl, maxIter = 10) {
    dist0 <- function(x,y,A) {
      tmp <- x - y
      tmp %*% A %*% tmp
    }
    dist1 <- function(x,y) {
      (x - y)^2
    }
    dist3 <- function(X,y,A) {
      tmp <- sweep(X,2,y,FUN='-')
      tmpD <- tmp %*% A
      tmpD <- rowSums(tmpD * tmp)
      tmpD
    }
    
    data <- as.matrix(data)
    N <- nrow(data)
    D <- ncol(data)
    cons_n <- unique(c(ml[,1],ml[,2],cl[,1],cl[,2]))  
    
    n_cl <- 0
    if (!is.null(cl))
      n_cl <- nrow(cl)
    n_ml <- 0
    if (!is.null(ml))
      n_ml <- nrow(ml)
    
    
    # initialization of constriant graph
    ml_graph <- list()
    cl_graph <- list()
    
    for (i in c(1:N)) {
      ml_graph[i] <- list(NULL)
      cl_graph[i] <- list(NULL)
    }
    # using adjacent graph to represent constraints
    if (n_ml > 0) {                                                                     
      for (i in c(1:n_ml)) {                                                              
        a <- ml[i,1]                                       
        b <- ml[i,2]                                         
        ml_graph[[a]] <- c(ml_graph[[a]],b)
        ml_graph[[b]] <- c(ml_graph[[b]],a)
      }                                                                                               
    }   
    if (n_cl > 0) {                                                                      
      for (i in c(1:n_cl)) {                                                              
        a <- cl[i,1]                                        
        b <- cl[i,2]                                        
        cl_graph[[a]] <- c(cl_graph[[a]],b)
        cl_graph[[b]] <- c(cl_graph[[b]],a)
      }                                                                                              
    }
    # remove duplicated points
    for (i in c(1:N)) {
      if (length(ml_graph[[i]]) > 1)
        ml_graph[[i]] <- unique(ml_graph[[i]])
      if (length(cl_graph[[i]]) > 1)
        cl_graph[[i]] <- unique(cl_graph[[i]])
    }
    # using depth first search to find transitive component in ml graph
    visited = rep(FALSE, N)
    neighborhoods = list()
    for (i in c(1:N)) {
      if (!visited[i] && length(ml_graph[[i]]) > 0) {
        component = c()
        
        stack <- c()
        # push the current node.  
        stack <- c(stack, i)  
        while (length(stack) > 0) {
          # pop a vertex from stack
          s <- stack[length(stack)]
          component <- c(component, s)
          stack <- head(stack, length(stack) - 1)
          # depth first search other nodes
          if (!visited[s]) {
            visited[s] <- TRUE
            if (length(ml_graph[[s]]) > 0)
              for (j in ml_graph[[s]])
                if (!visited[j])
                  stack <- c(stack, j)
          }
        }
        # add current component to the neighborhood results
        neighborhoods[[length(neighborhoods) + 1]] <- component
        for (a in component) 
          for (b in component) 
            if (a != b) 
              ml_graph[[a]] <- c(ml_graph[[a]], b)
      }
    }
    # remove duplicated points
    for (i in c(1:N)) {
      if (length(ml_graph[[i]]) > 1)
        ml_graph[[i]] <- unique(ml_graph[[i]])
    }
    # update cl 
    if (n_cl > 0) {                                                                      
      for (i in c(1:n_cl)) {                                                              
        a <- cl[i,1]                                        
        b <- cl[i,2]                                        
        for (x in ml_graph[[a]]) {
          cl_graph[[x]] <- c(cl_graph[[x]],b)
          cl_graph[[b]] <- c(cl_graph[[b]],x)
        }
        for (y in ml_graph[[b]]) {
          cl_graph[[y]] <- c(cl_graph[[y]],a)
          cl_graph[[a]] <- c(cl_graph[[a]],y)
        }
        for (x in ml_graph[[a]]) 
          for (y in ml_graph[[b]]) {
            cl_graph[[y]] <- c(cl_graph[[y]],x)
            cl_graph[[x]] <- c(cl_graph[[x]],y)
          }
      }                                                                                              
    }
    # remove duplicated points
    for (i in c(1:N)) {
      if (length(cl_graph[[i]]) > 1)
        cl_graph[[i]] <- unique(cl_graph[[i]])
    }
    ####################################centroid initialization###############################################
    lambda <- length(neighborhoods)
    # if lambda greater than k, choose k farthest-first traveral Np as centroids
    # otherwise, if an instance has cannot-link with all other Np, choose it as centroid
    # otherwise random select centroid
    fft <- matrix(nrow=lambda, ncol=3)
    fft[,2] <- rep(0,lambda)
    fft[,3] <- rep(.Machine$integer.max, lambda)
    centroids <- matrix(0, nrow=lambda, ncol=D)
    for (i in 1:lambda) {
      selectedData <- data[neighborhoods[[i]],]
      centroids[i,] <- colSums(selectedData) / nrow(selectedData)
      fft[i,1] <- length(neighborhoods[[i]])
    }
    tt <- sample(which(fft[,1]==max(fft[,1])),1)
    fft[tt,2] <- 1
    fft[tt,3] <- 0
    i <- 1
    while (i < min(lambda,K)) {
      i <- i+1
      selected <- which(fft[,2] == 1)
      candidates <- which(fft[,2] == 0)
      for (j in candidates) 
        for (tmpC in selected) 
          fft[j,3] <- min(fft[j,3], sum((centroids[j,]-centroids[tmpC])^2) * fft[j,1]) 
      
      tt <- which.max(fft[,3])
      fft[tt,2] <- 1
      fft[tt,3] <- 0
    }
    centroid <- centroids[which(fft[,2] == 1),]
    if (lambda < K) 
      centroid <- rbind(centroid, data[sample(1:N, K-lambda),])
    
    #####################################mpckm####################################################
    # initialize marix A
    
    A <- diag(D)
    label <- rep(0,N)
    Dist <- NULL
    for (i in 1:K) {
      tmpD <- dist3(data, centroid[i,], A)
      Dist <- cbind(Dist, tmpD)
    }
    Dist <- Dist - log(det(A)) / log(2)
    label <- apply(Dist, 1, which.min)
    
    for (iter in 1:maxIter) {
      #cat(iter,'\n')
      fC1 = matrix(0,nrow = K, ncol = 3)  
      # find maximally separated pair of points according to lith metric
      for (i in 1:N) 
        for (j in i:N) 
          if (label[i] == label[j]) {
            tmpL <- label[i]
            tmpD <- dist0(data[i,], data[j,], A)
            if (tmpD > fC1[tmpL,3]) 
              fC1[tmpL,] <- c(i,j,tmpD)
          }
      
      #cat('farthest pairs found\n')
      Dist <- NULL
      for (i in 1:K) {
        tmpD <- dist3(data, centroid[i,], A)
        Dist <- cbind(Dist, tmpD)
      }
      Dist <- Dist - log(det(A)) / log(2)
      
      label <- apply(Dist, 1, which.min)
      label[cons_n] <- 0
      #cat('without constraints points done\n')
      ord = sample(1:length(cons_n), length(cons_n))
      # start to assign data points related to constraints
      for (i in 1:length(cons_n)) {
        ss <- cons_n[ord[i]]
        
        if (length(cl_graph[[ss]]) > 0) 
          for (j in cl_graph[[ss]]) 
            if (label[j] != 0) {
              tmpD <- dist0(data[j,], data[ss,], A)
              
              Dist[ss, label[j]] <- Dist[ss, label[j]] + max(fC1[label[j],3] - tmpD,0)
            }
        if (length(ml_graph[[ss]]) > 0)
          for (j in ml_graph[[ss]])
            if (label[j] != 0) {
              tmpD <- dist0(data[j,], data[ss,], A)
              Dist[ss,] <- Dist[ss,] + rep(tmpD,K)
              Dist[ss,label[j]] <- Dist[ss,label[j]] - tmpD
            }
        
        # assign to smallest k 
        label[ss] = which.min(Dist[ss,])
      }
      if (iter == maxIter) 
        return(label)
      #cat('assignment done\n')
      # recalculate centroids
      for (i in c(1:K)) {
        temp <- data[which(label == i),]
        if (length(temp) > 0) {
          if (!is.null(nrow(temp)) )
            centroid[i,] <- colSums(temp) / nrow(temp)
          else
            centroid[i,] <- temp
        } else {
          centroid[i,] <- NA
        }
      }
      
      #cat('centroid caculated \n')
      tmp <- rep(0, D)
      # update metrics
      for (i in c(1:N)) {
        tmp <- tmp + dist1(data[i,],centroid[label[i],])
        if (length(cl_graph[[i]]) > 0) 
          for (j in cl_graph[[i]]) 
            if ((i < j) && (label[j] == label[i])) 
              if (fC1[label[j], 1] != 0) {
                temp <- dist1(data[fC1[label[j],1],],data[fC1[label[j],2],])
                temp[which(temp <0)] <- 0
                tmp <- tmp + temp
                
              }
        
        if (length(ml_graph[[i]]) > 0)
          for (j in ml_graph[[i]])
            if ((i < j) && (label[j] != label[i])) {
              tmp <- tmp + dist1(data[i,],data[j,])* 0.5
            }
      }
      
      A <- matrix(0, nrow = ncol(data), ncol = ncol(data))
      
      for (i in 1:ncol(data))
        A[i,i] <- 1* N / max(tmp[i],1e-9)
      #cat('metric updated\n')
    }
    return(label)
  }

# function to remove bad clusterings from base clusterings
# bad means clustering only contains one cluster
# or there is a super large cluster that containing more than 90% data instances
refine_base_clusts <- 
  function(BaseClusts) {
  rmv_Indx <- NULL
  n_instance <- nrow(BaseClusts)
  
  for (col in colnames(BaseClusts)) {
    Sizes <- table(BaseClusts[[col]])
    if (any(Sizes / n_instance >= 0.9))
      rmv_Indx <- c(rmv_Indx, col)
  }
  if (length(rmv_Indx) > 0)
  {
    cat("The following base clustering(s) is(are) bad and removed: ")
    for (col in rmv_Indx) {
      cat(col, table(BaseClusts[[col]]), '\n')
    }
    BaseClusts <- BaseClusts[,-which(names(BaseClusts) %in% rmv_Indx)]
  }
  return(BaseClusts)
}

# function to generate binary membership matrix based on base clusterings
binary_membership_matrix <- 
  function(BaseClusts) {
  rnum <- nrow(BaseClusts)
  rnames <- c(1:rnum)
  Clstrings <- ncol(BaseClusts)
  All_Matrix <- as.data.frame(matrix(nrow = rnum,ncol = 0,dimnames = list(rnames)))
  for(col in colnames(BaseClusts))
  {
    Cluster <- BaseClusts[[col]]
    CL <- unique(Cluster)
    Clust_Labls <- sort(CL)
    N_Clustrs <- length(CL)
    cnames <- paste(col,"C",Clust_Labls,sep="")
    BMatrx <- as.data.frame(matrix(data = F, nrow = rnum, ncol = N_Clustrs, dimnames = list(rnames,cnames)))
    for(i in 1:N_Clustrs)
    {BMatrx[Cluster==Clust_Labls[i],i]<-T}
    All_Matrix <- cbind(All_Matrix,BMatrx)
  }
  rm(i,j,rnum,rnames,Cluster,cnames,BMatrx,CL,N_Clustrs,Clust_Labls)
  
  return(All_Matrix)
}

# function to calculate closed pattern based on binary membership matrix
closed_frequent_pattern <- 
  function(All_Matrix, Clstrings) {
  
  options(stringsAsFactors = FALSE)
  n_instance <- nrow(All_Matrix)
  T1 <- Sys.time()
  
  Transactional <- as(All_Matrix, "transactions")
  
  A_Rules <- apriori(data = Transactional, parameter = list(support= 1/n_instance , target="closed frequent itemsets", confidence=1, maxlen=Clstrings))
  Combined <- supportingTransactions(x = A_Rules, transactions = Transactional)
  
  FCI <- data.frame(inspect(Combined),arules::size(Combined))
  colnames(FCI)<-c("ClosedSet","Object.List","Support")
  row.names(FCI) <- NULL
  FCI$ClosedSet <- substr(FCI$ClosedSet,2,nchar(FCI$ClosedSet)-1)   # remove { } from the string
  FCI$Object.List <- substr(FCI$Object.List,2,nchar(FCI$Object.List)-1)   # remove { } from the string
  FCI$ClosedSet <- strsplit(FCI$ClosedSet,split = ",")   # split the string into items
  FCI$Object.List <- strsplit(FCI$Object.List,split = ",")   # split the string into objects
  FCI$Object.List <- sapply(FCI$Object.List,as.integer)   # convert the objects into integers
  FCI$CS_size <- sapply(FCI$ClosedSet,length)
  Clstrings <- max(FCI$CS_size)
  FCI <- FCI[order(FCI$Support),]
  T2 <- Sys.time()
  rm(A_Rules,Combined,Transactional)
  cat("# of patterns=",nrow(FCI),'\n')
  Tme_FCI <- difftime(T2,T1,units = "sec")
  cat("FCP generation time",Tme_FCI,"sec.",'\n')
  rm(T1,T2)
  
  return(FCI)
}

# consensus function of multicons
consensus <- 
  function(FCI, Clstrings, Data_Size, MT=0.5) {
  cat("Consensus process starts, may need long time...\n")
  T1 <- Sys.time()
  # build first consensus
  Cons_Vctrs <- list()
  
  Bi_Clust <- FCI[FCI$CS_size == max(FCI$CS_size),"Object.List"]
  N_Row <- length(Bi_Clust)
  ClustV <- rep(NA, times= Data_Size)
  for(i in 1:N_Row)
  {
    ClustV[Bi_Clust[[i]]] <- i
  }
  Cons_Vctrs[[Clstrings]] <- ClustV
  cat("First consensus is built!\n")
  # build multiple consensuses:
  DT_range <- c((Clstrings-1):1)
  for (DT in DT_range)
  {
    
    
    Bi_Clust <- c(Bi_Clust,FCI[FCI$CS_size == DT,"Object.List"])
    N_Row <- length(Bi_Clust)
    
    for(i in 1:(N_Row-1))
    {
      Xi <- Bi_Clust[[i]]
      XiL <- length(Xi)
      if (XiL==0){next}
      for (j in (i+1):N_Row)
      {
        if (j==i){next}
        Xj <- Bi_Clust[[j]]
        XjL <- length(Xj)
        Intrs_Size <-length(intersect(Xi,Xj))
        if ( Intrs_Size==0 || XjL==0 ){next}
        Xkk <- intersect(Xi,Xj)
        Xkj <- setdiff(Xj,Xkk)
        Xki <- setdiff(Xi,Xkk)
        if (Intrs_Size==XiL)
        {
          Bi_Clust[[i]]<-integer(0)   # set Xi to empty to remove it later
          break
        }
        else if (Intrs_Size==XjL)
        {
          Bi_Clust[[j]]<-integer(0)   # set Xj to empty to remove it later
          next
        }
        else if ((Intrs_Size >= XiL*MT)||(Intrs_Size >= XjL*MT))
        { # merge bi_clusters Xi and Xj
          Bi_Clust[[j]] <- sort(union(Xi,Xj))
          Bi_Clust[[i]] <- integer(0)
          break
        }  else  
        { # split bi_clusters Xi and Xj
          if (XiL <= XjL)
          { 
            Bi_Clust[[j]] <- Xkj
          }
          else
          {
            Xi <- Bi_Clust[[i]] <- Xki
            XiL <- length(Xi)
          }
        }
      } # end for j
    } # end for i
    # Build final clusters
    Bi_Clust <- unique(Bi_Clust)
    Bi_Clust <- Bi_Clust[sapply(Bi_Clust,length)>0]
    N_Row <- length(Bi_Clust)
    ClustV <- NA
    for(i in 1:N_Row)
    {
      Indx <- Bi_Clust[[i]]
      if (all(is.na(ClustV[Indx]))) {ClustV[Indx] <- i} else {cat("error! cluster vector overwritten at DT =",DT,", final Bi-clusters are not unique.",'\n')}
    }
    Cons_Vctrs[[DT]] <- ClustV
    cat("One consensus is built! Current: ", Clstrings - DT, " Total: ", Clstrings, "\n")
  }
  
  rm(Bi_Clust,ClustV,DT,Intrs_Size,Xi,Xj,XiL,XjL,N_Row,Rpt_Chk,Indx)
  rm(Xk, Xki, Xkj, Xkk)
  
  T2 <- Sys.time()
  Tme_consensus <- difftime(T2,T1,units = "sec")
  cat("Consensus process time",Tme_consensus,"sec.",'\n')
  Cons_Vctrs <- Cons_Vctrs[lapply(Cons_Vctrs,length)>0]
  res <- data.frame(matrix(unlist(Cons_Vctrs), ncol=length(Cons_Vctrs), byrow=FALSE))
  for (i in 1:ncol(res)) {
    colnames(res)[i] <- paste("X",6-i,sep="")
  }
  return(res)
}

# consensus function of semi-multicons
constrained_consensus <- 
  function(FCI, Clstrings, Data_Size, MT=0.5, Must_Link, Cannot_Link) {
  cons_n <- c()
  if (!is.null(Must_Link) && (nrow(Must_Link) > 0)) {
    cons_n <- c(Must_Link[,1],Must_Link[,2])
  }
  if (!is.null(Cannot_Link) && (nrow(Cannot_Link) > 0)) {
    cons_n <- c(cons_n,Cannot_Link[,1],Cannot_Link[,2] )
  }
  cons_n <- unique(cons_n)   
  
  idx_list = hash()                                                                                 #
  .set(idx_list, keys=as.character(c(1:Data_Size)),values=length(cons_n)+1)                         #
  .set(idx_list, keys=as.character(cons_n), values=1:length(cons_n))                                #
  #
  constraints = matrix(0L, nrow=(length(cons_n) + 1), ncol=(length(cons_n) + 1))                    #
  #
  ml_list = hash()                                                                                  #
  n_ml <- 1                                                                                         #
  #
  if (!is.null(Must_Link)) {                                                                        #
    for (i in c(1:nrow(Must_Link))) {                                                               #
      #
      a <- toString(idx_list[[toString(Must_Link[i,1])]])                                         #
      b <- toString(idx_list[[toString(Must_Link[i,2])]])                                         #
      #
      if (is.null(ml_list[[a]])) {                                                                  #
        if (is.null(ml_list[[b]])) {                                                                #
          .set(ml_list, keys=b,values=n_ml)                                                         #
          .set(ml_list, keys=a,values=n_ml)                                                         #
          n_ml <- n_ml + 1                                                                          #
        } else {                                                                                    #
          .set(ml_list, keys=a,values=ml_list[[b]])                                                 #
        }                                                                                           #
      } else {                                                                                      #
        if (is.null(ml_list[[b]])) {                                                                #
          .set(ml_list, keys=b,values=ml_list[[a]])                                                 #
        } else {                                                                                    #
          ll <- c()                                                                                 #
          for (k in keys(ml_list)) {                                                                #
            if (ml_list[[k]] == ml_list[[a]]) {                                                     #
              ll <- rbind(ll,k)                                                                     #
            }                                                                                       #
          }                                                                                         #
          .set(ml_list, keys=ll, values=ml_list[[b]])                                               #
        }                                                                                           #
      }                                                                                             #
    }                                                                                               #
  }                                                                                                 #
  ml_list <- invert(ml_list)                                                                        #
  for (k in keys(ml_list)) {                                                                        #
    closed_set <- as.numeric(ml_list[[k]])                                                          #
    for (i in closed_set) {                                                                         #
      for (j in closed_set) {                                                                       #
        constraints[i,j] <- 1                                                                       #
      }                                                                                             #
    }                                                                                               #
  }                                                                                                 #
  for (i in c(0:length(cons_n) + 1))                                                                #
    constraints[i,i] <- 1                                                                           #
  if (!is.null(Cannot_Link))                                                                        #
    for (i in c(1:nrow(Cannot_Link))) {                                                             #
      a <- idx_list[[toString(Cannot_Link[i,1])]]                                                 #
      b <- idx_list[[toString(Cannot_Link[i,2])]]                                                 #
      #
      tmp1 <- constraints[b,]==1                                                                    #
      tmp2 <- constraints[a,]==1                                                                    #
      constraints[tmp2, tmp1] <- -1                                                                 #
      constraints[tmp1, tmp2] <- -1                                                                 #
    }                                                                                               #
  for (i in c(0:length(cons_n) + 1))                                                                #
    constraints[i,i] <- 0 
  
  
  cat("Consensus process starts, may need long time...\n")
  T1 <- Sys.time()
  # build first consensus
  Cons_Vctrs <- list()
  
  Bi_Clust <- FCI[FCI$CS_size == max(FCI$CS_size),"Object.List"]
  N_Row <- length(Bi_Clust)
  ClustV <- rep(NA, times= Data_Size)
  for(i in 1:N_Row)
  {
    ClustV[Bi_Clust[[i]]] <- i
  }
  Cons_Vctrs[[Clstrings]] <- ClustV
  cat("First consensus is built!\n")
  # build multiple consensuses:
  DT_range <- c((Clstrings-1):1)
  for (DT in DT_range)
  {
    Bi_Clust <- c(Bi_Clust,FCI[FCI$CS_size == DT,"Object.List"])
    N_Row <- length(Bi_Clust)
    #repeat
    #{
    # Rpt_Chk <- F
    for(i in 1:N_Row)
    {
      Xi <- Bi_Clust[[i]]
      XiL <- length(Xi)
      if (XiL==0){next}
      for (j in 1:N_Row)
      {
        if (j==i){next}
        Xj <- Bi_Clust[[j]]
        XjL <- length(Xj)
        Intrs_Size <-length(intersect(Xi,Xj))
        if ( Intrs_Size==0 || XjL==0 ){next}
        
        scorek <- 0
        scorekj <- 0
        scoreki <- 0
        
        Xkk <- intersect(Xi,Xj)
        Xkj <- setdiff(Xj,Xkk)
        Xki <- setdiff(Xi,Xkk)
        
        temp_kkk <- intersect(Xkk, cons_n)
        temp_jjj <- intersect(Xkj, cons_n)
        temp_iii <- intersect(Xki, cons_n)
        
        scorek <- 0
        scoreki <- 0
        scorekj <- 0
        
        if ((length(temp_iii) > 0) && (length(temp_kkk) > 0) && (length(Xki) > 0) && (length(Xkk) > 0)) {
          scoreki = -sum(constraints[unname(values(idx_list, keys=as.character(Xkk))),
                                     unname(values(idx_list, keys=as.character(Xki)))])
        }
        if ((length(temp_jjj) > 0) && (length(temp_kkk) > 0) && (length(Xkj) > 0) && (length(Xkk) > 0)) {
          scorekj = -sum(constraints[unname(values(idx_list, keys=as.character(Xkk))),
                                     unname(values(idx_list, keys=as.character(Xkj)))]) 
        }
        if ((length(temp_iii) > 0) && (length(temp_jjj) > 0) && (length(Xki) > 0) && (length(Xkj) > 0)) {
          scorek = sum(constraints[unname(values(idx_list, keys=as.character(Xki))),
                                   unname(values(idx_list, keys=as.character(Xkj)))])
        }
        if (scorek == 0 && scoreki == 0 && scorekj == 0) {
          if (Intrs_Size==XiL)
          {
            Bi_Clust[[i]]<-integer(0)   # set Xi to empty to remove it later
            break
          }
          else if (Intrs_Size==XjL)
          {
            Bi_Clust[[j]]<-integer(0)   # set Xj to empty to remove it later
            next
          }
          else if ((Intrs_Size >= XiL*MT)||(Intrs_Size >= XjL*MT))
          { # merge bi_clusters Xi and Xj
            #Rpt_Chk <- T
            Bi_Clust[[j]] <- sort(union(Xi,Xj))
            Bi_Clust[[i]] <- integer(0)
            break
          }  else  
          { # split bi_clusters Xi and Xj
            #Rpt_Chk <- T
            if (XiL <= XjL)
            { 
              Bi_Clust[[j]] <- Xkj
            }
            else
            {
              Xi <- Bi_Clust[[i]] <- Xki
              XiL <- length(Xi)
            }
          }
        } else {
          scorek <- scorek / (length(Xki) + length(Xkj))
          scoreki <- scoreki / length(Xi)# / length(Xkk) * 0.5
          scorekj <- scorekj / length(Xj) #/ length(Xkk) * 0.5
          if ((scorek  >= scoreki)&& (scorek >= scorekj))
          { # merge bi_clusters Xi and Xj
            Bi_Clust[[j]] <- sort(union(Xi,Xj))
            Bi_Clust[[i]] <- integer(0)
            break
          }  
          else  
          { # split bi_clusters Xi and Xj
            #print("split")
            if (scoreki <= scorekj)
            { 
              Bi_Clust[[j]] <- Xkj
            }
            else
            {
              Xi <- Bi_Clust[[i]] <- Xki
              XiL <- length(Xi)
            }
          }
        }
        
      } # end for j
    } # end for i
    # Build final clusters
    Bi_Clust <- unique(Bi_Clust)
    Bi_Clust <- Bi_Clust[sapply(Bi_Clust,length)>0]
    N_Row <- length(Bi_Clust)
    ClustV <- NA
    for(i in 1:N_Row)
    {
      Indx <- Bi_Clust[[i]]
      if (all(is.na(ClustV[Indx]))) {ClustV[Indx] <- i} else {cat("error! cluster vector overwritten at DT =",DT,", final Bi-clusters are not unique.",'\n')}
    }
    Cons_Vctrs[[DT]] <- ClustV
    cat("One consensus is built! Current: ", Clstrings - DT, " Total: ", Clstrings, "\n")
    
  }
  
  rm(Bi_Clust,ClustV,DT,Intrs_Size,Xi,Xj,XiL,XjL,N_Row,Rpt_Chk,Indx)
  rm(Xk, Xki, Xkj, Xkk)
  
  T2 <- Sys.time()
  Tme_consensus <- difftime(T2,T1,units = "sec")
  cat("Consensus process time",Tme_consensus,"sec.",'\n')
  Cons_Vctrs <- Cons_Vctrs[lapply(Cons_Vctrs,length)>0]
  res <- data.frame(matrix(unlist(Cons_Vctrs), ncol=length(Cons_Vctrs), byrow=FALSE))
  for (i in 1:ncol(res)) {
    colnames(res)[i] <- paste("X",6-i,sep="")
  }
  return(res)
}

# function to implement semi-multicons, BaseClusts represents the base clusterings calculated
# from our util function, or directly input from user
# the BaseClusts should be a *dataframe*
# an example of BaseClusts
#    mpckmeans2 mpckmeans3 mpckmeans4 mpckmeans5 mpckmeans6
#1            1          2          1          4          5
#2            1          2          2          2          2
#3            1          2          2          2          2
#4            1          2          2          2          2
#5            1          2          2          4          4
#6            1          2          1          4          5
#7            1          2          2          2          4
#8            1          2          2          2          4
#9            1          2          2          2          2
#10           1          2          2          2          2
# each column represents a clustering result, and must be assigned a unique name as column name
# each row represents a data instance
# in the example, the first row means the data instance no.1 belongs to cluster 1 for mpckmeans2, cluster 2 for mpckmeans3 etc..
# the first column represent the clustering results of mpckmeans2 for data instance no.1 to no.10
# MT by default set to 0.5
# Must_Link/Cannot_Link represent the must-link/cannot-link constraints generated from our util function
# or directly input from user
# the format of Must_Link and Cannot_Link should be the same as ml and cl for mpckemans function
# that means
# Must_Link and Cannot_Link must be *matrix* type
# where each row represent a pairwise constraint
# an example of cannot-link constraints Cannot_Link
#      [,1] [,2]
# [1,]   18  103
# [2,]  136   58
# [3,]  137   89
# [4,]   20  136
# [5,]    5  127
# the first row represent that data instance no.18 has a cannot-link with no.103, where index of data instance starts from *1* (not 0! R index starts with 1)
semi_multicons <- function(BaseClusts, MT, Must_Link, Cannot_Link) {

  BaseClusts <- refine_base_clusts(BaseClusts)
  
  All_Matrix <- binary_membership_matrix(BaseClusts)
  
  FCI <- closed_frequent_pattern(All_Matrix, ncol(BaseClusts))
  
  res <- constrained_consensus(FCI, ncol(BaseClusts), nrow(BaseClusts), MT, Must_Link, Cannot_Link)
  
  return(res)
}

# function to impelemnt multicons, BaseClusts represents the base clusterings calculated
# from our util function, or directly input from user
# the BaseClusts should be a *dataframe*
# an example of BaseClusts
#    mpckmeans2 mpckmeans3 mpckmeans4 mpckmeans5 mpckmeans6
#1            1          2          1          4          5
#2            1          2          2          2          2
#3            1          2          2          2          2
#4            1          2          2          2          2
#5            1          2          2          4          4
#6            1          2          1          4          5
#7            1          2          2          2          4
#8            1          2          2          2          4
#9            1          2          2          2          2
#10           1          2          2          2          2
# each column represents a clustering result, and must be assigned a unique name as column name
# each row represents a data instance
# in the example, the first row means the data instance no.1 belongs to cluster 1 for mpckmeans2, cluster 2 for mpckmeans3 etc..
# the first column represent the clustering results of mpckmeans2 for data instance no.1 to no.10
# MT by default set to 0.5
multicons <- function(BaseClusts, MT) {
  BaseClusts <- refine_base_clusts(BaseClusts)
  
  All_Matrix <- binary_membership_matrix(BaseClusts)
  
  FCI <- closed_frequent_pattern(All_Matrix, ncol(BaseClusts))
  
  res <- consensus(FCI, ncol(BaseClusts), nrow(BaseClusts), MT)
  
  return(res)
}