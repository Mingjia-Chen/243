#' Class mfa
#'
#' Class \code{mfa} defines a mfa model
#'
#' @name mfa-class
#' @rdname mfa-class
#' @exportClass mfa
setClass(
  Class="mfa",
  slots=list(
    sets="list",
    weights="numeric",
    eigenvalues="numeric",
    common_factor_score="matrix",
    partial_factor_score="list",
    loadings="matrix")
)

############ constructor function: building the model #################
#' Constructor method of mfa Class.
#'
#' @name mfa
#' @rdname mfa-class
#' @title MFA
#' @description Creates an object of class \code{mfa}.
#' @param data could be a numeric matrix or data.frame, should be in the same order of sets
#' @param sets list of vector contains vector of indices or variable names of each group
#' @param ncomps positive integer indicating how many number of components are to be extracted
#' @param center either a logical value or a numeric vector of length equal to the number of active variables in the analysis
#' @param scale either a logical value or a numeric vector of length equal to the number of active variables in the analysis
#' @return an object of class \code{mfa} with the following elements:
#' @return \item{sets}{list of sets}
#' @return \item{weights}{vector of weights}
#' @return \item{eigenvalues}{vector of eigenvalues}
#' @return \item{common_factor_score}{matrix of common factor scores}
#' @return \item{partial_factor_score}{list of matrix of partial factor scores}
#' @return \item{loadings}{matrix of loadings}
#' @export
#' @examples
#' # default
#' test<-mfa(wine_data,sets=list(1:6,7:12,13:18,19:23,24:29,30:34,35:38,39:44,45:49,50:53))
#'
#' # use your own scale method
#' #scale with vector
#' test<-mfa(wine_data,sets=list(1:6,7:12,13:18,19:23,24:29,30:34,35:38,39:44,45:49,50:53),
#'           center = apply(wine_data,2,mean),scale = apply(wine_data,2,sd))
#' # scale in advance
#' ndatas<-apply(wine_data,2,function(x){ (x-mean(x))/norm(x-mean(x),type="2")})
#' test<-mfa(ndatas,sets=list(1:6,7:12,13:18,19:23,24:29,30:34,35:38,39:44,45:49,50:53),center=FALSE,scale=FALSE)
#'
#' # only print the first two components
#' ndatas<-apply(wine_data,2,function(x){ (x-mean(x))/norm(x-mean(x),type="2")})
#' test<-mfa(ndatas,sets=list(1:6,7:12,13:18,19:23,24:29,30:34,35:38,39:44,45:49,50:53),ncomp=2,center=FALSE,scale=FALSE)
#'
#' # character sets
#' ndatas<-apply(wine_data,2,function(x){ (x-mean(x))/norm(x-mean(x),type="2")})
#' # Use the first and last variable names of each group:
#' test<-mfa(ndatas,sets=list(c("V1","V6"),c("V1.1","V8"),c("V1.2","V10"),c("V1.3","V8.1"),
#'                            c("V1.4","V12"),c("V1.5","V13"),c("V1.6","V4.6"),
#'                            c("V1.7","V5.1"),c("V1.8","V15"),c("V1.9","V4.9")
#'                            ),scale = F, center = F)
#' # or use the full list of variable names of each group:
#' test<-mfa(ndatas,sets=list(c("V1","V2","V3","V4","V5","V6"),
#'                            c("V1.1","V2.1","V3.1","V4.1","V7","V8"),
#'                            c("V1.2","V2.2","V3.2","V4.2","V9","V10"),
#'                            c("V1.3","V2.3","V3.3","V4.3","V8.1"),
#'                            c("V1.4","V2.4","V3.4","V4.4","V11","V12"),
#'                            c("V1.5","V2.5","V3.5","V4.5","V13"),
#'                            c("V1.6","V2.6","V3.6","V4.6"),
#'                            c("V1.7","V2.7","V3.7","V4.7","V14","V5.1"),
#'                            c("V1.8","V2.8","V3.8","V4.8","V15"),
#'                            c("V1.9","V2.9","V3.9","V4.9")
#'                            ),scale = F, center = F)



# constructor function: to construct 'mfa' and run the model to get attributes
# parameter: data: could be a matrix or a data.frame, should be in the same order of sets
# parameter: sets: list of vector contains vector of indices of each group
# eg. sets=c(1:3,4:5), means the 1:3 columns of data is Group 1 and the next 4:5 columns is Group2
# center and scale: the same parameters as in the function scale(), logical values or a numeric vector
mfa<-function(data,sets,ncomps=NULL,center=TRUE,scale=TRUE){
  #check numeric data
  data<-as.matrix(data)
  if(!is.numeric(data)) {stop("data should be numeric matrix or data.frame")}
  
  #check center and scale
  if(!is.logical(center)&& !(is.numeric(center) && (length(center) == ncol(data)))) {
    stop("center should be either a logical value or a numeric vector of length equal to the number of columns of 'data'")
  }
  if(!is.logical(scale)&& !(is.numeric(scale) && (length(scale) == ncol(data)))) {
    stop("scale should be either a logical value or a numeric vector of length equal to the number of columns of 'data'")
  }
  if(!is.logical(scale)&&any(scale==0)) {stop("scale vector can't contain zero values.")}
  
  #check sets
  if(is.numeric(sets[[1]])){
    check_sets<-NULL
    for (i in 1:length(sets)) {
      if(!is.numeric(sets[[i]])) {stop("sets should be a list of numeric vectors or character vectors.")}
      check_sets<-c(check_sets,sets[[i]])
    }
    if(length(check_sets)!=ncol(data)) {stop("The sum of sets lengths does not equal to the number of columns.")}
    if(any(!check_sets%in%1:ncol(data))) {stop("sets out of bounds")}
    if(!identical(1:ncol(data),check_sets)) {warning("sets contain some overlapped and skipped columns.")}
  }else{
    if(is.character(sets[[1]])){
      check_sets<-NULL
      check_names<-NULL
      for (i in 1:length(sets)){
        if(!is.character(sets[[i]])) {stop("sets should be a list of character vectors or numeric vectors.")}
        if(any(!sets[[i]]%in%colnames(data))) {stop("sets contain wrong variable names.")}
        check_names<-c(check_names,sets[[i]])
        check_sets<-c(check_sets,c(which(colnames(data)==sets[[i]][1]):which(colnames(data)==sets[[i]][length(sets[[i]])])))
      }
      if(length(check_sets)!=ncol(data)) {stop("The sum of sets lengths does not equal to the number of columns, or the variable names are in the wrong order.")}
      if(!identical(colnames(data),check_names)) {warning("sets contain some overlapped and skipped columns.")}
    }else{
      stop("sets should be a list of numeric vectors or character vectors.")
    }
  }
  
   
  datarownames<-row.names(data)
  
  # scale and center
  data<-scale(data,center,scale)
  
  # check singularity and ncomps
  rank<-Matrix::rankMatrix(data)[1]
  if (!is.null(ncomps)){
    if(ncomps<=0||ncomps%%1!=0) {
      stop("ncomps should be positive integer.")
    }
    if ((ncomps)>rank){
      warning(paste0("ncomps larger than rank: outputing ",rank," dimensions."))
    }else{
      rank<-ncomps
    }
  }else{
    if (!rank%in%c(dim(data))){
      message(paste0("Matrix is singular: outputing ",rank," dimensions."))
    }
  }
  
  # if sets is character: turn sets into indicies acccording to rownames of data
  osets<-sets
  if (!is.numeric(sets[[1]])){
    newlist<-list()
    for (i in 1:length(sets)){
      newlist[[i]]<-c(which(colnames(data)==sets[[i]][1]):which(colnames(data)==sets[[i]][length(sets[[i]])]))
    }
    sets<-newlist
  }
  
  # divide data into several group according to values in sets
  # store the ith group of data to variable "Groupi"
  for (i in 1:length(sets)) {
    assign(paste0("Group",i),data.matrix(data[,min(sets[[i]]):max(sets[[i]])]))
  }
  
  
  # for each data groups conduct svd
  # store the first singular values in singularvalues
  singularvalues<-c(rep(1,dim(data)[1]))
  for (i in 1:length(sets)) {
    singularvalues[i]<-max(svd(eval(parse(text=paste0("Group",i))))$d)
  }
  
  
  # construct A to compute Q: QAQt=I
  # A's diagonal elements are the inverse of the first square singular values
  # each first square singular value is expanded to the same dimension of each data group
  expanded<-c()
  for (i in 1:length(sets)){
    expanded<-c(expanded,rep(singularvalues[i],max(sets[[i]])-min(sets[[i]])+1))
  }
  weights<-1/expanded^2
  A<-diag(x = 1/expanded^2,length(expanded),length(expanded))
  A_half<-diag(x = 1/expanded,length(expanded),length(expanded))
  A_half_inv<-diag( x = expanded,length(expanded),length(expanded))
  
  
  
  # construct M to compute P: PMPt=I
  M<-diag(x=1/(dim(data)[1]),dim(data)[1],dim(data)[1])
  M_half<-diag(x=1/sqrt((dim(data)[1])),dim(data)[1],dim(data)[1])
  M_half_inv<-diag(x=sqrt((dim(data)[1])),dim(data)[1],dim(data)[1])
  
  
  #  X: the whole data
  X<-data.matrix(data)
  # Construct S=XAXt
  S<-as(X %*% A %*% t(X),"matrix")
  
  # do spectral decomposition on S: S=P*LAMBDA*Pt, PtMP=I
  # construct inverse delta: delta^2=LAMBDA
  eigens<-eigen(S)
  d<-matrix(0,dim(X)[1],dim(X)[1])
  for (i in 1:length(eigens$values)){
    d[i,i]<-eigens$values[i]
  }
  u<-eigens$vectors
  lambda<-as(M_half %*% d %*% M_half,"matrix")
  
  delta_value<-diag(as(sqrt(lambda),"matrix"))
  
  delta_inv<-1/sqrt(lambda)
  delta_inv[is.infinite(delta_inv)]<-0
  
  # P is PMPt=I FOR S=P*LAMBDA*Pt
  P <- as(u %*% M_half_inv,"matrix")
  # Q FOR Q=Xt*M*P*DELTA_inverse
  Q <- as(t(X) %*% M %*% P %*% delta_inv, "matrix")[,1:rank]
  dimnames(Q) <- list(rownames(Q),colnames(Q, do.NULL = FALSE, prefix = "Dim"))
  
  
  # build a list: 'partial_factor_score' to store partial factor score
  # build a matrix: 'common_factor_score' to store common factor score
  # the partial score of group i is named "Partial Score: Group i"
  # partial factor score i = no. of group * A_i* data group i* Q_i
  # common factor score = sum of partial factor score i
  
  
  partial_factor_score<-list()
  common_factor_score<-0
  for (i in 1:length(sets)){
    datai<- data.matrix(eval(parse(text=paste0("Group",i))))
    score<-length(sets) * (1/singularvalues[i]^2) * datai %*% t(datai)  %*% M %*% P %*% delta_inv[,1:rank]
    dimnames(score) <- list(datarownames,colnames(score, do.NULL = FALSE, prefix = "Dim"))
    partial_factor_score[[paste0("Partial Score: Group ",i)]]=as(score,"matrix")
    common_factor_score<-score+common_factor_score
  }
  common_factor_score<-common_factor_score/length(sets)
  
  
  # loading uses Q
  new (Class = "mfa",
       sets=osets,
       eigenvalues = c(delta_value^2)[1:rank],
       weights= weights,
       common_factor_score = as(common_factor_score,"matrix"),
       partial_factor_score = partial_factor_score,
       loadings = as(Q,"matrix")
  )
}

