net.fun=function(interest, interest.label, res1,res2,res3){
  feature_id=c(interest,res1$feature_id, res2$feature_id, res3$feature_id)
  label=c(interest.label, res1$feature_desc, res2$feature_desc, res3$feature_desc)
  nodes <- data.frame(id = 1:length(feature_id), feature_id=feature_id)
  nodes$id=as.numeric(nodes$id)
  nodes$feature_id=as.character(nodes$feature_id)
  nodes$font.size = 25
  nodes[which(nodes$feature_id %in% interest),]$font.size = 40
  nodes$font.bold = TRUE
  nodes$shape <- "box"
  nodes$shadow <- TRUE
  nodes$size <- 50
  nodes$color.background = NA
  nodes[which(substr(nodes$feature_id,1,3) == "Phe"),]$color.background = "pink"
  nodes[which(substr(nodes$feature_id,1,3) == "RXN"),]$color.background = "yellow"
  nodes[which(substr(nodes$feature_id,1,3) %in% c("LOI","Sho", "Oth")),]$color.background = "orange"
  nodes[which(substr(nodes$feature_id,1,3) == "CCS"),]$color.background = "skyblue"
  nodes[which(nodes$feature_id %in% interest),]$color.background = "red"
  nodes$title=label
  nodes$label=label
  nodes$borderWidth <- 0.3
  nodes$color.border <- "black"
  
  #nodes[which(nodes$feature_id %in% res2$feature_id),]$color.background = "white"
  #nodes[which(nodes$feature_id %in% res3$feature_id),]$color.background = "white"
  
  nodes$group="Both"
  nodes[nodes$feature_id%in%c(res2$feature_id)==1,"group"]="RPDR_only"
  nodes[nodes$feature_id%in%c(res3$feature_id)==1,"group"]="VA_only"

  
  nodes$site=NA
  nodes$site[nodes$group=="Both"]=paste("RPDR", "VA", "Both", sep=",")
  nodes$site[nodes$group=="VA_only"]=paste("VA", "Both", sep=",")
  nodes$site[nodes$group=="RPDR_only"]=paste("RPDR", "Both", sep=",")
 
  code.keep=rbind(res1,res2,res3)[,c(1,3)]
  code.keep=cbind(interest, code.keep)
  edges=cbind(1,2:length(feature_id),rbind(res1,res2,res3)[,3])
  edges=data.frame(edges)
  colnames(edges)=c("from", "to", "weight")
  edges=edges[edges$weight!=0,]
  edges$width <- 1 
  edges$color <- "black"
  edges$shadow <- FALSE
  
  edges[edges$to%in%nodes[nodes$color.background=="white", "id"]|edges$from%in%nodes[nodes$color.background=="white", "id"], "color"]="lightblue"
  nodes=nodes[nodes$id%in%c(edges$from, edges$to),]
  
  color.list = cbind(color = c("red", "skyblue", "yellow", "orange", "pink"), type = c("interest", "PROC", "MED","LAB", "DX"))
  lnodes <- data.frame(color = color.list[,1],
                       label = color.list[,2],
                       shape = "box")
  #, selected = "Both")
  # Network
  net <- visNetwork(nodes, edges, width = "100%", height = 670) %>% visLayout(18) %>% 
    visLegend(addNodes = lnodes, useGroups = F, width = 0.07) %>%
    visGroups(groupname = "Site", color="white")%>%
    visOptions(highlightNearest = list(enabled = T, labelOnly = F),
               selectedBy = list(variable = c("site"), multiple=T)) %>%
    visPhysics(barnesHut = list(avoidOverlap = 0)) %>%
    visInteraction(navigationButtons = T,dragNodes = FALSE) 
  net
}

norm.function=function(x1,x2){
  x1.new=x1/sqrt(sum(x1^2))
  x2.new=x2/sqrt(sum(x2^2))
  sqrt(sum((x1.new-x2.new)^2))
}

cos.function=function(x1,x2){
  x1.new=x1/sqrt(sum(x1^2))
  x2.new=x2/sqrt(sum(x2^2))
  crossprod(x1.new, x2.new)
}

level.fun=function(x){
  if(grepl("[.]",x)==0){level=0}else{level=nchar(strsplit(x,"[.]")[[1]][2])}
  level
}

integer.fun=function(x){
  if(grepl("[.]",x)==0){integer=x}else{integer=strsplit(x,"[.]")[[1]][1]}
  integer
}
integer.fun2=function(x){
  x=gsub("PheCode.","", x)
  if(grepl("[.]",x)==0){integer=x}else{integer=strsplit(x,"[.]")[[1]][1]}
  paste0("PheCode.",integer)
}


cos.in.FUN=function(phecode.keep){
  rxnorm.in=phecode.rxnorm[phecode.rxnorm$phecode==phecode.keep,"rxnorm"]
  cos.in= unlist(lapply(1:length(rxnorm.in), function(ll) crossprod(P.est.A[,phecode.keep],P.est.A[,rxnorm.in[ll]])/
                          sqrt(crossprod(P.est.A[,phecode.keep])*crossprod(P.est.A[,rxnorm.in[ll]]))))
  cos.in
}

cos.out.FUN=function(phecode.keep){
  rxnorm.in=phecode.rxnorm[phecode.rxnorm$phecode==phecode.keep,"rxnorm"]
  rxnorm.out=unique(setdiff(phecode.rxnorm[,"rxnorm"], rxnorm.in))
  cos.out= unlist(lapply(1:length(rxnorm.out), function(ll) crossprod(P.est.A[,phecode.keep],P.est.A[,rxnorm.out[ll]])/
                           sqrt(crossprod(P.est.A[,phecode.keep])*crossprod(P.est.A[,rxnorm.out[ll]]))))
  cos.out
}

cos.in.dis.FUN=function(phecode.keep, phecode.keep.int, phecode.keep.mtx){
  phecode.in=setdiff(phecode.keep.mtx[phecode.keep.mtx[,"int"]==phecode.keep.int,"code"],phecode.keep)
  cos.in= unlist(lapply(1:length(phecode.in), function(jj) crossprod(P.est.A[,phecode.keep],P.est.A[,phecode.in[jj]])/
                          sqrt(crossprod(P.est.A[,phecode.keep])*crossprod(P.est.A[,phecode.in[jj]]))))
  cos.in
}

cos.out.dis.FUN=function(phecode.keep, phecode.keep.int, phecode.keep.mtx){
  phecode.in=setdiff(phecode.keep.mtx[phecode.keep.mtx[,"int"]==phecode.keep.int,"code"],phecode.keep)
  phecode.out=unique(setdiff(phecode.keep.mtx[,"code"], phecode.in))
  cos.out= unlist(lapply(1:length(phecode.out), function(ll) crossprod(P.est.A[,phecode.keep],P.est.A[,phecode.out[ll]])/
                           sqrt(crossprod(P.est.A[,phecode.keep])*crossprod(P.est.A[,phecode.out[ll]]))))
  cos.out
}


digit.fun=function(x){
  if(grepl("[.]",x)==0){digit=NA}else{integer=strsplit(x,"[.]")[[1]][1]; tmp=strsplit(x,"[.]")[[1]][2]; tmp2=substr(tmp,1,1); digit=paste(integer, tmp2, sep=".")}
  digit
}

digit.fun2=function(x){
  x=gsub("PheCode.","", x)
  if(grepl("[.]",x)==0){digit=NA}else{integer=strsplit(x,"[.]")[[1]][1]; tmp=strsplit(x,"[.]")[[1]][2]; tmp2=substr(tmp,1,1); digit=paste(integer, tmp2, sep=".")}
  paste0("PheCode.",digit)
}

cos.fun=function(interest, feature.list1, P.est.A, cos.thrsh){
  interest.idx=which(feature.list1==interest)
  all.cos <- top.cos <- list()
  cos.sims <- unlist(lapply(c(1:ncol(P.est.A)),function(j) {
    crossprod(P.est.A[,interest.idx],P.est.A[,j])/
      sqrt(crossprod(P.est.A[,interest.idx])*crossprod(P.est.A[,j]))}))
  
  all.cos=data.frame(feature.list1,cos.sims)
  all.cos[,1]=as.character(all.cos[,1])
  top.cos=cos.sims>cos.thrsh
  return(list(all.cos=all.cos, top.cos=top.cos))
}


filter.fun=function(interest, feature.list1, P.est.A, all.new, explain.match, cos.thrsh, freq.thrsh,coocur.thrsh){
  interest.idx=which(feature.list1==interest)
  junk=cos.fun(interest, feature.list1, P.est.A, cos.thrsh)
  top.cos=junk$top.cos
  all.cos=junk$all.cos
  all_phen <- all.new[all.new$code1==interest,]
  
  exclude.indx <- unique(c(interest.idx,which(!top.cos),which(is.na(abs(all.cos$cos.sims))),
                           which(feature.list1 %in% all_phen[all_phen$count<coocur.thrsh,'code2']),
                           which(feature.list1 %in% explain.match[explain.match$freq_single<freq.thrsh,'feature_id'])
                           ))
  exclude.indx
}

filter.fun.new=function(interest, feature.list1, P.est.A, cos.thrsh){
  interest.idx=which(feature.list1==interest)
  junk=cos.fun(interest, feature.list1, P.est.A, cos.thrsh)
  top.cos=junk$top.cos
  all.cos=junk$all.cos

  exclude.indx <- unique(c(interest.idx,which(!top.cos),which(is.na(abs(all.cos$cos.sims)))))
  exclude.indx
}

B.lasso.fun=function(interest, P.est.A, P.est.sam, P.est.sam2, all.cos, feature.list1, exclude.indx, nlambda, n_words, alpha){
  interest.idx=which(feature.list1==interest)
  ## step 1: grid.lambda
  set.seed(1234)
  obj1 <- cv.glmnet(x = P.est.A[,-exclude.indx], y = P.est.A[,interest.idx],
                   alpha = alpha,intercept=F,standardize=F,
                   penalty.factor = 1 / abs(all.cos$cos.sim)[-exclude.indx]
  )
  max.lambda=obj1$lambda.min*1e2; min.lambda=obj1$lambda.min*1e-6
  grid.lambda = seq(from = max.lambda, to = min.lambda, length.out = nlambda)
  obj.A <- glmnet(x = P.est.A[,-exclude.indx], y = P.est.A[,interest.idx],
                alpha = alpha, lambda = grid.lambda,intercept=F,standardize=F,
                penalty.factor = 1 / abs(all.cos$cos.sim[-exclude.indx]))
  obj.sam <- glmnet(x = P.est.sam[,-exclude.indx], y = P.est.sam[,interest.idx],
                alpha = alpha, lambda = grid.lambda,intercept=F,standardize=F,
                penalty.factor = 1 / abs(all.cos$cos.sim[-exclude.indx]))
  
  B.sam = Matrix(0, n_words, nlambda)
  B.sam[-exclude.indx,] = as.matrix(coef(obj.sam))[-1,]
  
  B.A = Matrix(0, n_words, nlambda)
  B.A[-exclude.indx,] = as.matrix(coef(obj.A))[-1,]
 
  MSE = unlist(lapply(1:nlambda, function(ll) mean((P.est.sam2[,interest.idx] - P.est.sam2 %*% matrix(B.sam[,ll],ncol=1))^2)))
  B.lasso = as.matrix(B.A[,which.min(MSE)])
  list(B.lasso=B.lasso, MSE=MSE, grid.lambda=grid.lambda)
}

drop.dat.fun=function(x, y, drop.rate, up.rate=10){
  tmpind = sample(1:nrow(x),up.rate*nrow(x),replace=T)
  ynew = y[tmpind]; xnew = x[tmpind,]; 
  colnames(xnew) = colnames(x)
  col.drop = 1: ncol(xnew)
  xnew[,col.drop] = sapply(col.drop,function(kk){tmpx=xnew[,kk]; tmpx[sample(1:length(tmpx),round(drop.rate*length(tmpx)))] = mean(tmpx); tmpx})
  list(ynew=ynew, xnew=xnew)
}

B.dropout.fun=function(interest, P.est.A, P.est.sam, P.est.sam2, all.cos, feature.list1, exclude.indx, nlambda, n_words, alpha, up.rate=10, drop.rate=0.5){
  interest.idx=which(feature.list1==interest)
  ## step 1: grid.lambda
  set.seed(1234)
  junk1=drop.dat.fun(P.est.A[,-exclude.indx], P.est.A[,interest.idx], drop.rate, up.rate)
  xnew.A=junk1$xnew
  ynew.A=junk1$ynew
  cv.fit1 = cv.glmnet(xnew.A,ynew.A,family="gaussian",intercept=F,alpha=alpha,standardize=F,
                      penalty.factor = 1 / abs(all.cos$cos.sim[-exclude.indx]))
  max.lambda=cv.fit1$lambda.min*1e2; min.lambda=cv.fit1$lambda.min*1e-6
  grid.lambda = seq(from = max.lambda, to = min.lambda, length.out = nlambda)
  
  ## step 2: training
  set.seed(1234)
  junk2.sam=drop.dat.fun(P.est.sam[,-exclude.indx], P.est.sam[,interest.idx], drop.rate, up.rate)
  xnew.sam=junk2.sam$xnew
  ynew.sam=junk2.sam$ynew
  fit.sam = glmnet(xnew.sam,ynew.sam,family="gaussian",intercept=F,alpha=alpha,standardize=F,
                     penalty.factor = 1 / abs(all.cos$cos.sim[-exclude.indx]), lambda=grid.lambda)

  junk2.A=drop.dat.fun(P.est.A[,-exclude.indx], P.est.A[,interest.idx], drop.rate, up.rate)
  xnew.A=junk2.A$xnew
  ynew.A=junk2.A$ynew
  fit.A = glmnet(xnew.A,ynew.A,family="gaussian",intercept=F,alpha=alpha,standardize=F,
                   penalty.factor = 1 / abs(all.cos$cos.sim[-exclude.indx]), lambda=grid.lambda)

  B.sam = Matrix(0, n_words, nlambda)
  B.sam[-exclude.indx,] = as.matrix(coef(fit.sam))[-1,]
  
  B.A = Matrix(0, n_words, nlambda)
  B.A[-exclude.indx,] = as.matrix(coef(fit.A))[-1,]
  
  ## step 3: validation
  MSE = unlist(lapply(1:nlambda, function(ll) mean((P.est.sam2[,interest.idx] - P.est.sam2 %*% matrix(B.sam[,ll],ncol=1))^2)))
  B.dropout = as.matrix(B.A[,which.min(MSE)])
  list(B.dropout=B.dropout, MSE=MSE, grid.lambda=grid.lambda)
}

B.lasso.fun.backup2=function(interest, P.est.A, P.est.sam, P.est.sam2, all.cor, feature.list1, exclude.indx, nlambda, n_words){
  interest.idx=which(feature.list1==interest)
  #max.lambda.ori=1e6
  grid.lambda1 = seq(from = 0.01, to = 0, length.out = 100)
  obj1 <- cv.glmnet(x = P.est.A[,-exclude.indx], y = P.est.A[,interest.idx],
                    alpha = 0.8,intercept=F,standardize=F,
                    penalty.factor = 1 / abs(all.cor)[-exclude.indx]
  )
  #max.lambda = min(obj1$lambda.min*100, max.lambda.ori)
  max.lambda=obj1$lambda.min*10; min.lambda=obj1$lambda.min*0.01
  #max.lambda=0.00022; min.lambda=0.0002
  
  ###
  # Lasso estimate (1st lasso)
  #grid.lambda = seq(from = max.lambda, to = 0.0001*max.lambda, length.out = nlambda)
  grid.lambda = seq(from = max.lambda, to = min.lambda, length.out = nlambda)
  obj.A <- glmnet(x = P.est.A[,-exclude.indx], y = P.est.A[,interest.idx],
                  alpha = 0.8, lambda = grid.lambda,intercept=T,standardize=F,
                  penalty.factor = 1 / abs(all.cor[-exclude.indx]))
  obj.sam <- glmnet(x = P.est.sam[,-exclude.indx], y = P.est.sam[,interest.idx],
                    alpha = 0.8, lambda = grid.lambda,intercept=T,standardize=F,
                    penalty.factor = 1 / abs(all.cor[-exclude.indx]))
  
  beta.sam = Matrix(0, n_words, nlambda)
  beta.sam[-exclude.indx,] = as.matrix(coef(obj.sam))[-1,]
  
  beta.A = Matrix(0, n_words, nlambda)
  beta.A[-exclude.indx,] = as.matrix(coef(obj.A))[-1,]
  ###
  B.A = NULL
  beta.list = unlist(beta.A)
  for(order in 1:nlambda){
    B.temp = beta.list[,order]
    names(B.temp) = feature.list1
    B.A = cbind(B.A, B.temp)
  }
  
  B.sam = NULL
  beta.list = unlist(beta.sam)
  for(order in 1:nlambda){
    B.temp = beta.list[,order]
    names(B.temp) = feature.list1
    B.sam = cbind(B.sam, B.temp)
  }
  # calculate AIC (not accurate here, the penalty is L0 norm) -----------------
  AIC.lasso = matrix(NA, nlambda, 1)
  min.AIC.lasso = 1e10
  order.min.lambda = 0
  MSE_reg.mat <- MSE.mat <- AIC.mat <- rep(NA, nlambda)
  for(i in 1:nlambda){
    B.final.sam = B.sam[,i]
    n.select=sum(B.final.sam!=0)
    obj = picasso(P.est.sam[,B.final.sam!=0], P.est.sam[,interest.idx],
                  lambda = grid.lambda[i], type.gaussian="naive",
                  intercept=F)
    
    MSE_reg.mat[i] <- mean((P.est.sam2[,interest.idx]-P.est.sam2[,B.final.sam!=0] %*% obj$beta)^2)
    res = P.est.sam2[,interest.idx] - P.est.sam2 %*% matrix(B.final.sam,ncol=1)
    MSE.mat[i] <- mean(res^2)
    SSR = sum(res %*% t(res))
    AIC.lasso[i,1] = d*log(SSR) #+ 0.5*length(Select.name)
    AIC.mat[i] <- d*log(MSE.mat[i]) + 0.5*n.select# checa bien la formula de AIC
  }
  B.lasso = as.matrix(B.A[,which.min(MSE.mat)])
  list(B.lasso=B.lasso, MSE=MSE.mat, MSE_reg.mat=MSE_reg.mat, AIC.mat=AIC.mat, AIC.lasso=AIC.lasso, grid.lambda=grid.lambda)
}

B.lasso.fun.backup=function(interest, P.est.A, P.est.sam, all.cor, feature.list1, exclude.indx, nlambda, n_words){
  interest.idx=which(feature.list1==interest)
  #max.lambda.ori=1e6
  grid.lambda1 = seq(from = 0.01, to = 0, length.out = 100)
  obj1 <- cv.glmnet(x = P.est.A[,-exclude.indx], y = P.est.A[,interest.idx],
                    alpha = 1,intercept=F,standardize=F,
                    penalty.factor = 1 / abs(all.cor)[-exclude.indx]
  )
  #max.lambda = min(obj1$lambda.min*100, max.lambda.ori)
  max.lambda=obj1$lambda.min*100; min.lambda=obj1$lambda.min
  #max.lambda=0.00022; min.lambda=0.0002
  
  ###
  # Lasso estimate (1st lasso)
  #grid.lambda = seq(from = max.lambda, to = 0.0001*max.lambda, length.out = nlambda)
  grid.lambda = seq(from = max.lambda, to = min.lambda, length.out = nlambda)
  obj <- glmnet(x = P.est.A[,-exclude.indx], y = P.est.A[,interest.idx],
                alpha = 1, lambda = grid.lambda,intercept=T,standardize=F,
                penalty.factor = 1 / abs(all.cor[-exclude.indx]))
  
  beta = Matrix(0, n_words, nlambda)
  ###
  beta[-exclude.indx,] = as.matrix(coef(obj))[-1,]
  ###
  B = NULL
  beta.list = unlist(beta)
  for(order in 1:nlambda){
    B.temp = beta.list[,order]
    names(B.temp) = feature.list1
    B = cbind(B, B.temp)
  }
  # calculate AIC (not accurate here, the penalty is L0 norm) -----------------
  AIC.lasso = matrix(NA, nlambda, 1)
  min.AIC.lasso = 1e10
  order.min.lambda = 0
  MSE_reg.mat <- MSE.mat <- AIC.mat <- rep(NA, nlambda)
  for(i in 1:nlambda){
    B.final = B[,i]
    n.select=sum(B.final!=0)
    obj = picasso(P.est.sam[,B.final!=0], P.est.sam[,interest.idx],
                  lambda = grid.lambda[i], type.gaussian="naive",
                  intercept=F)
    
    MSE_reg.mat[i] <- mean((P.est.sam[,interest.idx]-P.est.sam[,B.final!=0] %*% obj$beta)^2)
    res = P.est.sam[,interest.idx] - P.est.sam %*% matrix(B.final,ncol=1)
    MSE.mat[i] <- mean(res^2)
    SSR = sum(res %*% t(res))
    AIC.lasso[i,1] = d*log(SSR) #+ 0.5*length(Select.name)
    AIC.mat[i] <- d*log(MSE.mat[i]) + 0.5*n.select# checa bien la formula de AIC
  }
  B.lasso = as.matrix(B[,which.min(MSE.mat)])
  list(B.lasso=B.lasso, MSE=MSE.mat, grid.lambda=grid.lambda)
}

