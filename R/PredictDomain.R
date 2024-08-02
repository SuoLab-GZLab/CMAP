#' Data preprocess
#'
#' @param sc_norm
#' @param st_norm
#' @param spatial_genes
#' @param batch
#' @param npc
#' @param nclust
#' @param pca_method
#' @param max.iter.harmony
#' @param max.iter.cluster
#' @param cosine_norm
#'
#' @return
#' @export
#'
#' @examples
data_to_transform = function(sc_norm,
                             st_norm,
                             spatial_genes,
                             batch=TRUE,
                             npc= 100,nclust= 100,pca_method = c("prcomp", "prcomp_irlba"),
                             max.iter.harmony = 50,
                             max.iter.cluster = 6, cosine_norm = FALSE

){
  pca_method = match.arg(pca_method, choices = c('prcomp','prcomp_irlba'))
  genes <- intersect(spatial_genes,rownames(sc_norm))
  st_qm = preprocessCore::normalize.quantiles(as.matrix(st_norm[genes,,drop=FALSE]),keep.names = TRUE)
  sc_qm = preprocessCore::normalize.quantiles(as.matrix(sc_norm[genes,,drop=FALSE]),keep.names = TRUE)
  if(batch){
    matrix <- harmony_embedding(st_norm = st_qm,
                                sc_norm = sc_qm,
                                genes = genes,
                                batch = batch,
                                npc = min(npc,dim(st_qm)[2],dim(sc_qm)[2]),
                                nclust=min(nclust,ceiling(dim(sc_qm)[2]/10),ceiling(dim(st_qm)[2]/10)),
                                pca_method = pca_method,
                                max.iter.harmony = max.iter.harmony,
                                max.iter.cluster = max.iter.cluster,
                                cosine_norm = cosine_norm)
  }else{
    matrix <- cbind(sc_qm,st_qm)
  }
  return(matrix)
}





#' @title Scale
#' @description Scale data to [0,1] as SVM input
#' @param train_set spot*gene
#' @param test_set cell*gene
#' @examples
#' @export
#' @return
# a list including train and test matrics
scale_data_by_column <- function(train_set,test_set){
  tmp.train = train_set
  tmp.test = test_set

  for(i in 1:dim(test_set)[2]){
    min = min(train_set[,i,drop=FALSE])
    max = max(train_set[,i,drop=FALSE])
    tmp.test[,i] = (1/(max-min))*(test_set[,i]-min)
    tmp.train[,i] = (1/(max-min))*(train_set[,i]-min)
  }
  data <- list(tmp.train,tmp.test)
  return(data)
}


#'  Predict the spatial domain
#'
#' @param train_set spot*(gene+label)
#' @param test_set cell*gene
#' @param scale whether scale data to [0,1]
#' @param class.weight whether considering data imbalance
#' @param cost
#' @param gamma
#' @param verbose
#' @param kernel
#' @param st_svm
#'
#' @return
#' @export
#'
#' @examples
PredictDomain <- function(train_set, test_set, scale=TRUE, class.weight=TRUE,
                          cost=1, gamma=1/ncol(test_set), kernel = "radial",
                          st_svm=FALSE, verbose=FALSE){
  # scale to [0,1]
  if(scale){
    scale_data = scale_data_by_column(train_set,test_set)
    tmp.train = scale_data[[1]]
    tmp.test = scale_data[[2]]
  }

  # consider the unbalanced data
  if(class.weight){
    num_cluster <- length(unique(tmp.train$label))
    wts = dim(tmp.train)[1]/table(tmp.train$label)/length(unique(tmp.train$label))
  }else{
    wts = rep(1,length(unique(tmp.train$label)))
    names(wts) <- names(table(tmp.train$label))
  }

  # set parameters
  set.seed(321)
  svmfit = e1071::svm(label~., data=tmp.train,scale = FALSE,gamma = gamma,cost = cost,
                      probability = TRUE,kernel=kernel,class.weight=wts)
  # compute decision values and probabilites
  pred_st_svm <- predict(svmfit,tmp.train,scale = FALSE,probability = T)
  # compute decision values and probabilites
  pred_sc_svm <- predict(svmfit,tmp.test,scale = FALSE,probability = T)
  if(verbose){
    print("The prediction table of training data")
    print(table(pred_st_svm))
    print(paste0("The accuracy of spatial data:",sum(diag(table(pred_st_svm,train_set$label)))/dim(train_set)[1]))
    print("The prediction table of test data")
    print(table(pred_sc_svm))
  }
  if(st_svm){
    return(pred_st_svm)
  }
  return(pred_sc_svm)
}


#' Tune the parameters
#'
#' @param train_set spot*(gene+label)
#' @param test_set cell*gene
#' @param scale whether scale data to [0,1]
#' @param class.weight whether considering data imbalance
#' @param verbose
#' @param kernel
#' @param cross_para
#'
#' @return
#' @export
#'
#' @examples
tune_parameter <- function(train_set,test_set,scale=TRUE,class.weight=TRUE,
                           kernel="radial",verbose=FALSE, cross_para=c(4,6,8,10)){
  parameter <- list()

  # scale to [0,1]
  if(scale){
    scale_data = scale_data_by_column(train_set,test_set)
    tmp.train = scale_data[[1]]
    tmp.test = scale_data[[2]]
  }else{
    print("Please set Scale TRUE")
    tmp.train = train_set
    tmp.test = test_set
  }

  # consider the unbalanced data
  if(class.weight){
    num_cluster <- length(unique(tmp.train$label))
    wts = dim(tmp.train)[1]/table(tmp.train$label)/length(unique(tmp.train$label))
  }else{
    wts = rep(1,length(unique(tmp.train$label)))
    names(wts) <- names(table(tmp.train$label))
  }

  # tune the parameters
  for(cross_i in cross_para){
    set.seed(321)
    # C is the penalty coefficient, i.e., the tolerance for error. The higher the c, the less tolerant the error is, and the more prone to overfitting. The smaller the C, the easier it is to underfit. If C is too large or too small, the generalization ability becomes poor
    tune.out.corse=e1071::tune(svm, label~.,data=tmp.train ,kernel=kernel,
                               ranges =list(cost=2^seq(-5,15,2),
                                            gamma=2^seq(-15,3,2)), probability = TRUE,class.weight=wts,
                               scale = FALSE, tunecontrol=tune.control(cross=cross_i))
    bestmodel <- tune.out.corse$best.model
    parameter[[paste0('cross_',cross_i)]][['cost']] <- bestmodel$cost
    parameter[[paste0('cross_',cross_i)]][['gamma']] <- bestmodel$gamma
    if(verbose){
      print(paste0("SVM: cross=",cross_i))
      print("best parameters:")
      print(paste0("cost:",bestmodel$cost))
      print(paste0("gamma:",bestmodel$gamma))
    }
  }
  return(parameter)
}
