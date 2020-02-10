## hyperbinning
## by Jan Stuchly

HBmodel<-function(flowframe,depth=4,channels,plus=1000,log=TRUE){
    if (log){
        data<-exprs(flowframe)[,channels]+plus
        data[data<1]<-1
        data<-log10(data)
    } else {
        data<-exprs(flowframe)[,channels]
    }

    pres<-new.env()
    pres$val<-matrix(NA,nrow=length(channels)+1,ncol=2^depth-1)
    HB<-function(data,node,pres){
        N<-nrow(data)

        if (node>=2^depth) return(0)
        covdata<-cov(data)
        direction<-eigen(covdata)$vectors[,1]
        proj<-data%*%direction
        op<-order(proj)

        if (N%%2==0){
            x<-0.5*(data[op[N%/%2],]+data[op[N%/%2+1],])

            b<--t(direction)%*%x
            pres$val[,node]<-c(direction,b)

            HB(data[op[1:(N%/%2)],],2*node,pres)
            HB(data[op[(N%/%2+1):N],],2*node+1,pres)

        } else {
            x<-data[op[N%/%2+1],]
            b<--t(direction)%*%x
            pres$val[,node]<-c(direction,b)
            HB(data[op[1:(N%/%2)],],2*node,pres)
            HB(data[op[(N%/%2+1):N],],2*node+1,pres)
        }
    }
    HB(data,1,pres)
    return(list(hyperplanes=pres$val,channels=channels,nbbins=2^depth,plus=plus,log=log))
}

HBclassify<-function(model,flowframe){
    channels<-model$channels
    pres<-new.env()
    if (model$log){
        pres$data<-exprs(flowframe)[,channels]+model$plus
        pres$data[pres$data<1]<-1
        pres$data<-log10(pres$data)
    } else {
        pres$data<-exprs(flowframe)[,channels]
    }

    nc<-ncol(pres$data)
    pres$viz<-rep(0,nrow(pres$data))
    pres$val<-rep(0,model$nbbins)
    ## SS1<<-SS<<-0
    HBC<-function(model,pres,node,sel){
         ## classifier<-(model$hyperplanes[1:nc,node])%*%t(pres$data[sel,]) + model$hyperplanes[nc+1,node]
         ##    sel1<-which(classifier<0)

         ##    print(length(sel))
         ##    print(length(sel1))
         ##    print("--------------")

        if ((2*node)>=(model$nbbins) | (2*node+1)>=(model$nbbins)) {

           if (length(sel)==1) classifier<-(model$hyperplanes[1:nc,node])%*%(pres$data[sel,]) + model$hyperplanes[nc+1,node] else classifier<-(model$hyperplanes[1:nc,node])%*%t(pres$data[sel,]) + model$hyperplanes[nc+1,node]
            sel1<-which(classifier<0)
            #SS1<<-SS1+length(sel)
            if (length(sel)!=(length(sel1)+ifelse(length(sel1)>0,length(sel[-sel1]),length(sel)))){
                ## print(sel)
                ## print(sel1)
                ## uu<<-(pres$data[sel,])
                ## print(classifier)
            }
            #SS<<-SS+length(sel1)+ifelse(length(sel1)>0,length(sel[-sel1]),length(sel))

            pres$viz[sel[sel1]]<-2*node-model$nbbins+1
            if (length(sel1)>0) pres$viz[sel[-sel1]]<-2*node-model$nbbins+2 else pres$viz[sel]<-2*node-model$nbbins+2

            pres$val[2*node-model$nbbins+1]<-length(sel1)
            pres$val[2*node-model$nbbins+2]<-ifelse(length(sel1)>0,length(sel[-sel1]),length(sel))
            return(0)
        }

        if (length(sel)==1) classifier<-(model$hyperplanes[1:nc,node])%*%(pres$data[sel,]) + model$hyperplanes[nc+1,node] else classifier<-(model$hyperplanes[1:nc,node])%*%t(pres$data[sel,]) + model$hyperplanes[nc+1,node]
        classifier<-(model$hyperplanes[1:nc,node])%*%t(pres$data[sel,]) + model$hyperplanes[nc+1,node]
        sel1<-which(classifier<0)



        if (length(sel1)>0) HBC(model,pres,2*node,sel[sel1])
        if (length(sel1)==0) HBC(model,pres,2*node+1,sel) else HBC(model,pres,2*node+1,sel[-sel1])


    }
    sel<-1:nrow(pres$data)
    HBC(model,pres,1,sel)

    return(list(visualize=pres$viz,finger=pres$val))

}

power_k<-function(k){
    res<-matrix(0,2^k,k)
    for (i in 1:(2^k)){
        iloc<-i
        for (j in (k-1):0){
            if (iloc>2^j){
                res[i,j+1]<-1
                iloc<-iloc-2^j
            }


        }
    }
    ## res<-lapply(1:2^k,FUN=function(i) which(res[i,]==1))
    res[res==0]<--1
    return(res)

}

xor_power<-function(power){
    N<-nrow(power)
    done=FALSE
    while (!done){

        powloc<-t(t(power[2:N,])+power[1,])

        rem<-which(apply((powloc==1),MARGIN=1,FUN=all))
        print(powloc)
        print(rem)
        if (length(rem)==0) done=TRUE else power<-power[-1,]
        N<-nrow(power)
    }
    power<-lapply(1:N,FUN=function(i) which(power[i,]==1))
    return(power)
}

pops_vs_pac<-function(power,pac){
    Npo<-length(power)
    Npa<-pac
    res<-list()[1:Npa*Npo]
    for (i in 1:Npa){
        for (j in 1:Npo){
            res[[j+(i-1)*Npo]]<-c(i,j)
        }
    }
    return(res)

}

preprocess_data<-function(flowset,stains,stratification,maxevents,plus=100){
    if (is.character(flowset)){
        Nsamples<-length(flowset)
        clusters<-rep(NA,Nsamples)
        res<-NULL
        i<-1
        for (file in flowset){
            fcs<-read.FCS(file)
            if (is.character(stratification)){
                j<-which(names(fcs@description)==stratification)
                clusters[i]<-fcs@description[[j]]
            }
            fcs<-exprs(fcs)
            fcs<-fcs[1:min(nrow(fcs),maxevents),stains]
            fcs<-fcs+plus
            fcs[which(fcs<1)]<-1
            fcs<-log10(fcs)
            fcs<-data.frame(fcs,sample=i)
            i<-i+1
            res<-rbind(res,fcs)
        }
    }
    return(list(events=res,clusters=clusters))
}

f_objective<-function(params,data,vars){
    N<-ncol(data)-1
    pops<-vars$pops
    pac<-unique(data[,N+1])

    dirs<-matrix(params,nrow=N+1)
    k<-ncol(dirs)
    b<-dirs[N+1,]
    dirs<-dirs[-(N+1),]

    val<-t(t(data[,-(N+1)]%*%dirs)+b)
    val[val<0]<--1
    val[val>0]<-1
    val<-val%*%t(pops)
    val[val!=k]<-0
    val<-val/k

    res<-lapply(pac,FUN=function(p){
        sel<-which(data[,N+1]==p)

        r<-apply(val[sel,],MARGIN=2,FUN=function(x) sum(x)/length(sel))


    })
    res<-do.call("rbind",res)

    res<-data.frame(clusters=vars$clusters,res)

    ## res<-J48(clusters~.,data=res)
    ## correct<-length(which(res$prediction==vars$clusters))

    ## res<-correct/length(vars$clusters)
    res<-dunn_mod(Data=res[,-1],clusters=as.numeric(res[,1]),method="maximum")
    print(res)
    return(-res)

}

f_objective_fd<-function(params,data,vars){
    N<-ncol(data)-1
    pops<-vars$pops
    pac<-unique(data[,N+1])

    dirs<-matrix(vars$params,nrow=N+1)
    k<-ncol(dirs)
    b<-params
    dirs<-dirs[-(N+1),]

    val<-t(t(data[,-(N+1)]%*%dirs)+b)
    val[val<0]<--1
    val[val>0]<-1
    val<-val%*%t(pops)
    val[val!=k]<-0
    val<-val/k

    res<-lapply(pac,FUN=function(p){
        sel<-which(data[,N+1]==p)

        r<-apply(val[sel,],MARGIN=2,FUN=function(x) sum(x)/length(sel))


    })
    res<-do.call("rbind",res)

    res<-data.frame(clusters=vars$clusters,res)

    ## res<-J48(clusters~.,data=res)
    ## correct<-length(which(res$prediction==vars$clusters))

    ## res<-correct/length(vars$clusters)
    res<-dunn_mod(Data=res[,-1],clusters=as.numeric(res[,1]),method="maximum")
    print(res)
    return(-res)

}
f_objective_fb<-function(params,data,vars){
    N<-ncol(data)-1
    pops<-vars$pops
    pac<-unique(data[,N+1])

    dirs<-matrix(params,nrow=N)
    k<-ncol(dirs)
    b<-vars$params[N+1,]


    val<-t(t(data[,-(N+1)]%*%dirs)+b)
    val[val<0]<--1
    val[val>0]<-1
    val<-val%*%t(pops)
    val[val!=k]<-0
    val<-val/k

    res<-lapply(pac,FUN=function(p){
        sel<-which(data[,N+1]==p)

        r<-apply(val[sel,],MARGIN=2,FUN=function(x) sum(x)/length(sel))


    })
    res<-do.call("rbind",res)

    res<-data.frame(clusters=vars$clusters,res)

    ## res<-J48(clusters~.,data=res)
    ## correct<-length(which(res$prediction==vars$clusters))

    ## res<-correct/length(vars$clusters)
    res<-dunn_mod(Data=res[,-1],clusters=as.numeric(res[,1]),method="maximum")
    print(res)
    return(-res)

}

f_objective_eval<-function(params,data,vars){
    N<-ncol(data)-1
    pops<-vars$pops
    pac<-unique(data[,N+1])

    dirs<-matrix(params,nrow=N+1)
    k<-ncol(dirs)
    b<-dirs[N+1,]
    dirs<-dirs[-(N+1),]

    val<-t(t(data[,-(N+1)]%*%dirs)+b)
    val[val<0]<--1
    val[val>0]<-1
    val<-val%*%t(pops)
    val[val!=k]<-0
    val<-val/k

    res<-lapply(pac,FUN=function(p){
        sel<-which(data[,N+1]==p)

        r<-apply(val[sel,],MARGIN=2,FUN=function(x) sum(x)/length(sel))


    })
    res<-do.call("rbind",res)

    res<-data.frame(clusters=vars$clusters,res)


    return(res)

}
prepare_vars<-function(data,k){
    clusters<-data$clusters
    data<-data$events
    Nchan<-ncol(data)-1
    Npac<-length(unique(data[,Nchan+1]))
    pops<-power_k(k)


    pvp<-pops_vs_pac(pops,Npac)
    return(list(pops=pops,pvp=pvp,clusters=clusters))

}

dunn_mod<-function (distance = NULL, clusters, Data = NULL, method = "euclidean",epsilon=0.1)
{
    if (is.null(distance) & is.null(Data))
        stop("One of 'distance' or 'Data' is required")
    if (is.null(distance))
        distance <- as.matrix(dist(Data, method = method))
    if (class(distance) == "dist")
        distance <- as.matrix(distance)
    nc <- max(clusters)
    interClust <- matrix(NA, nc, nc)
    intraClust <- rep(NA, nc)
    for (i in 1:nc) {
        c1 <- which(clusters == i)
        for (j in i:nc) {
            if (j == i)
                intraClust[i] <- max(distance[c1, c1])
            if (j > i) {
                c2 <- which(clusters == j)
                interClust[i, j] <- min(distance[c1, c2])
            }
        }
    }
    dunn <- min(interClust, na.rm = TRUE)/(max(intraClust)+epsilon)
    return(dunn)
}

max_dist<-function (distance = NULL, clusters, Data = NULL, method = "euclidean")
{
    if (is.null(distance) & is.null(Data))
        stop("One of 'distance' or 'Data' is required")
    if (is.null(distance))
        distance <- as.matrix(dist(Data, method = method))
    if (class(distance) == "dist")
        distance <- as.matrix(distance)
    nc <- max(clusters)
    interClust <- matrix(NA, nc, nc)
    intraClust <- rep(NA, nc)
    for (i in 1:nc) {
        c1 <- which(clusters == i)
        for (j in i:nc) {
            if (j == i)
                intraClust[i] <- max(distance[c1, c1])
            if (j > i) {
                c2 <- which(clusters == j)
                interClust[i, j] <- min(distance[c1, c2])
            }
        }
    }
    dunn <- min(interClust, na.rm = TRUE)
    return(dunn)
}

init_perp_hyperplanes<-function(data,k){
    N<-ncol(data)-1
    data<-data[,-(N+1)]
    med<-rep(NA,N)
    for (i in 1:N) med[i]<-median(data[,i])
    covdata<-cov(data)
    direction<-eigen(covdata)$vectors
    k<-min(k,ncol(direction))
    direction<-direction[,1:k]
    med<- matrix(med,nrow=1)
    print(dim(med))
    print(dim(direction))
    b<--med%*%direction
    return(rbind(direction,b))


}

write.data.sne<-function(data,theta=0.5,perplexity=30){
    d<-ncol(data)-1
    N<-nrow(data)
    X<-as.double(t(as.matrix(data[,-(d+1)])))
    file<-file("data.dat","wb")
    writeBin(as.integer(N),file)
    writeBin(as.integer(d),file)
    writeBin(as.double(theta),file)
    writeBin(as.double(perplexity),file)
    writeBin(X,file)
    close(file)
}

read.data.sne<-function(file){
    print(file)
    file<-file(file,"rb")
    print(file)
    N<-readBin(file,"integer",1)
    d<-readBin(file,"integer",1)
    X<-readBin(file,"double",N*d)
    landmarks<-readBin(file,"integer",N)
    costs<-readBin(file,"double",N)
    X<-t(matrix(X,d,N))
    close(file)
    return(list(X=X,landmarks=landmarks+1,costs=costs))

}

run_tsne_bh<-function(events,theta=0.5,perplexity=30,fcs=NULL){
    write.data.sne(events$events,theta=theta,perplexity=perplexity)
    system("./bh_tsne")
    X<-read.data.sne("result.dat")$X
    X<-X+abs(min(X))+50
    colnames(X)<-c("1","2")
    res<-cbind(X,events$orig,sample=events$events[,ncol(events$events)]+5)
    if (!is.null(fcs)){
        fcsf<-flowFrame(res)
        write.FCS(fcsf,file=fcs)
    }
    return(res)
}

preprocess_data_sne<-function(flowset,stains,stratification,maxevents,original=NULL){
    if (is.character(flowset)){
        Nsamples<-length(flowset)
        clusters<-rep(NA,Nsamples)
        res<-NULL
        i<-1
        orig<-NULL
        jo<-1
        for (file in flowset){
            fcs<-read.FCS(file)
            ## fcs<-compensate(fcs,fcs@description$SPILL)
            if (is.null(original)) orig<-rbind(orig,exprs(fcs)[1:min(nrow(fcs),maxevents),stains]) else {
                fcso<-read.FCS(original[jo])
                orig<-rbind(orig,exprs(fcso)[1:min(nrow(fcso),maxevents),stains])
                jo<-jo+1
            }
            if (is.null(original)){
                ll<-estimateLogicle(x=fcs,channels=colnames(exprs(fcs))[stains])

                fcs<-transform(fcs,ll)
                if (is.character(stratification)){
                    j<-which(names(fcs@description)==stratification)
                    clusters[i]<-fcs@description[[j]]
                }
            } else {
                if (is.character(stratification)){
                    j<-which(names(fcso@description)==stratification)
                    clusters[i]<-fcso@description[[j]]
                }
            }
            fcs<-exprs(fcs)
            fcs<-fcs[1:min(nrow(fcs),maxevents),stains]

            fcs<-data.frame(fcs,sample=i)
            i<-i+1
            res<-rbind(res,fcs)
        }
    }
    return(list(events=res,clusters=clusters,orig=orig))
}

HBDmodel<-function(flowframe,depth=4,channels,targetN=100,plus=1000,log=TRUE){
    require(feature)
     if (log){
        data<-exprs(flowframe)[,channels]+plus
        data[data<1]<-1
        data<-log10(data)
    } else {
        data<-exprs(flowframe)[,channels]
    }

    pres<-new.env()
    pres$val<-matrix(NA,nrow=length(channels)+1,ncol=2^depth-1)



    HB<-function(data,node,pres,targetN){
        N<-nrow(data)
        if(is.vector(data) | (node>=2^depth) ) return(0)
        if (N<=targetN) return(0)
        covdata<-cov(data)
        direction<-eigen(covdata)$vectors[,1]
        proj<-data%*%direction


        split<-get_split_point(proj)
        left<-which(proj<=split)
        right<-which(proj>split)

        l<-which.max(proj[proj<=split])[1]
        p<-which.min(proj[proj>split])[1]
        x<-0.5*(data[left[l],]+data[right[p],])
        b<--t(direction)%*%x
        pres$val[,node]<-c(direction,b)
        HB(data[left,],2*node,pres,targetN)
        HB(data[right,],2*node+1,pres,targetN)

    }
    HB(data,1,pres,targetN)
    return(list(hyperplanes=pres$val,channels=channels,nbbins=2^depth,plus=plus,log=log))
}

HBDclassify<-function(model,flowframe){
    channels<-model$channels
    pres<-new.env()

    if (model$log){
        pres$data<-exprs(flowframe)[,channels]+model$plus
        pres$data[pres$data<1]<-1
        pres$data<-log10(pres$data)
    } else {
        pres$data<-exprs(flowframe)[,channels]
    }


    nc<-ncol(pres$data)
    pres$viz<-rep(0,nrow(pres$data))
    pres$val<-rep(0,model$nbbins)

    HBC<-function(model,pres,node,sel){

        if ((2*node)>=(model$nbbins) | (2*node+1)>=(model$nbbins) | ((2*node)<=ncol(model$hyperplanes)) && is.na(model$hyperplanes[1,2*node])) {

            classifier<-(model$hyperplanes[1:nc,node])%*%t(pres$data[sel,]) + model$hyperplanes[nc+1,node]
            sel1<-which(classifier<0)

            pres$viz[sel[sel1]]<-2*node-model$nbbins+1
            pres$viz[sel[-sel1]]<-2*node-model$nbbins+2

            pres$val[2*node-model$nbbins+1]<-length(sel1)
            pres$val[2*node-model$nbbins+2]<-length(sel[-sel1])
            return(0)
        }

        classifier<-(model$hyperplanes[1:nc,node])%*%t(pres$data[sel,]) + model$hyperplanes[nc+1,node]
        sel1<-which(classifier<0)
        HBC(model,pres,2*node,sel[sel1])
        HBC(model,pres,2*node+1,sel[-sel1])


    }
    sel<-1:nrow(pres$data)
    HBC(model,pres,1,sel)

    return(list(visualize=pres$viz,finger=pres$val))

}



get_split_point<-function(x,bw,tau=5,d=1,gridsize=401,signifLevel=0.05){

    n <- length(x)

    if (missing(bw))           ## b/w not specified -> interactive
        {
            ##bw.range <- dfltBWrange(x,gridsize,tau, scale.fac=1.5)
            bw.range <- feature:::dfltBWrange(as.matrix(x),tau)

            bw <- matrix(unlist(bw.range), nrow=2, byrow=FALSE)
            dfltCounts.out <- dfltCounts(x,gridsize, apply(bw, 2, max))

            h.low <- bw[1,]
            h.upp <- bw[2,]
            hmix.prop <- 1/4
            h.init <- h.low^(hmix.prop)*h.upp^(1-hmix.prop) ##sqrt(h.low*h.upp)
            dfltCounts.out <- dfltCounts(x,gridsize, h.init)
            h <- h.init
        }
    else
        {
            dfltCounts.out <- dfltCounts(x,gridsize, bw)
            h <- bw
        }
    gcounts <- dfltCounts.out$counts
    range.x <- dfltCounts.out$range.x

    dest <- drvkde(gcounts, drv=0, bandwidth=h, binned=TRUE, range.x=range.x, se=FALSE, gridsize=gridsize)
    dest$est[dest$est<0] <- 0

    ESS <- n*dest$est*prod(h)*(sqrt(2*pi)^d)
    SigESS <- ESS >= 5

    Sig.scalar <- array(NA, dim=gridsize)
    Sig2.scalar <- array(NA, dim=gridsize)

    dest$est[dest$est<0] <- 0
    ## constant for variance of gradient estimate
    Sig.scalar <- 1/2*(2*sqrt(pi))^(-d)*n^(-1)*prod(h)^(-1)*dest$est
    Sig2.scalar <- (8*sqrt(pi)*n*prod(h))^(-1)*dest$est

                                        #gradient dignificance
    obj1 <- drvkde(gcounts, drv=1, bandwidth=h, binned=TRUE, range.x=range.x, se=FALSE)
    fhat1 <- obj1$est

    Sig.inv12 <- 1/sqrt(Sig.scalar * h^(-2))
    WaldGrad <- (Sig.inv12 * fhat1)^2

                                        #curvature significance
    obj2 <- drvkde(gcounts,drv=2,bandwidth=h,binned=TRUE,range.x=range.x, se=FALSE)
    fhat2 <- obj2$est

    Sig2.inv12 <- 1/sqrt(Sig2.scalar * 3*h^(-4))
    lambda1 <- Sig2.inv12 * fhat2
    WaldCurv <- lambda1^2
    local.mode <- (lambda1 < 0)

    ## multiple hypothesis testing - based on Hochberg's method
    ## - modified Bonferroni method using ordered p-values

    ## test statistic for gradient
    if (TRUE)
        {
            pval.Grad <- 1 - pchisq(WaldGrad, d)
            pval.Grad.ord <- pval.Grad[order(pval.Grad)]
            num.test <- sum(!is.na(pval.Grad.ord))

            if (num.test>=1)
                num.test.seq <- c(1:num.test, rep(NA, prod(gridsize) - num.test))
            else
                num.test.seq <- rep(NA, prod(gridsize))

            reject.nonzero <- ((pval.Grad.ord <= signifLevel/(num.test + 1 - num.test.seq)) &
                               (pval.Grad.ord > 0))
            reject.nonzero.ind <- which(reject.nonzero)

            ## p-value == 0 => reject null hypotheses automatically
            SignifGrad <- array(FALSE, dim=gridsize)
            SignifGrad[which(pval.Grad==0, arr.ind=TRUE)] <- TRUE

            ## p-value > 0 then reject null hypotheses indicated in reject.nonzero.ind
            for (i in reject.nonzero.ind)
                SignifGrad[which(pval.Grad==pval.Grad.ord[i], arr.ind=TRUE)] <- TRUE
        }

    ## test statistic for curvature
    if (TRUE)
        {
            pval.Curv <- 1 - pchisq(WaldCurv, d*(d+1)/2)
            pval.Curv.ord <- pval.Curv[order(pval.Curv)]
            num.test <- sum(!is.na(pval.Curv.ord))

            if (num.test>=1)
                num.test.seq <- c(1:num.test, rep(NA, prod(gridsize) - num.test))
            else
                num.test.seq <- rep(NA, prod(gridsize))
            reject.nonzero <- ((pval.Curv.ord <= signifLevel/(num.test + 1 - num.test.seq)) &(pval.Curv.ord > 0))
            reject.nonzero.ind <- which(reject.nonzero)

            SignifCurv <- array(FALSE, dim=gridsize)

            ## p-value == 0 => reject null hypotheses automatically
            SignifCurv[which(pval.Curv==0, arr.ind=TRUE)] <- TRUE

            ## p-value > 0 then reject null hypotheses indicated in reject.nonzero.ind
            for (i in reject.nonzero.ind)
                SignifCurv[which(pval.Curv==pval.Curv.ord[i], arr.ind=TRUE)] <- TRUE

            neg <- SignifCurv & !local.mode
            pos<- SignifCurv & local.mode
        }
    min.neg.gradpoint<-min(dest$x.grid[[1]][SignifGrad & (obj1$est<0)])
    max.pos.gradpoint<-max(dest$x.grid[[1]][SignifGrad & (obj1$est>0)])
    q10<-quantile(x,probs=c(0.1,0.9))
    feasible.neg.curv<-which(dest$x.grid[[1]]>=max(q10[1],min.neg.gradpoint) & dest$x.grid[[1]]<=min(q10[2],max.pos.gradpoint) & neg)
    split<-dest$x.grid[[1]][feasible.neg.curv][which.min(dest$est[feasible.neg.curv])]
    if (length(feasible.neg.curv)==0){


        op<-order(x)
        if (n%%2==0) split<-0.5*(x[op[n%/%2]]+x[op[n%/%2+1]]) else split<-x[op[n%/%2+1]]

     }
    ## print(n)
    ## print(length(feasible.neg.curv)==0)
    ## plot(dest$est~dest$x.grid[[1]],type="l")
    ## abline(v=split)
    ## readline()

    return(split)

    return(list(grad=SignifGrad, curv.neg=neg,curv.pos=pos,dest=dest,range=range.x,grad=obj1,curv=obj2,h=h,split=split))


}

vizualize.model<-function(fcs,model,pdf="model.pdf",bins=NULL){
    N<-length(model$channels)
    clas<-HBclassify(model,fcs)

    data<-exprs(fcs)+model$plus
    pch<-rep(".",nrow(data))
    if (!is.null(bins)) pch[which(clas$visualize %in% bins)]<-"\2"
    data[data<1]<-1
    data<-log10(data)
    pdf(pdf)
    for (i in 1:(N-1)){
        for(j in (i+1):N){
            plot(data[,c(i,j)],col=clas$visualize+1,pch=pch)
        }
    }
    dev.off()

}


batch.classify<-function(files, model){
    res1<-res2<-NULL
    for (i in files){
        print(i)
        fcs<-rename.extract(i,model$channels)
        fcs<-flowFrame(fcs)

        res<-HBclassify(flowframe=fcs,model=model)
        res1<-cbind(res1,res$visualize)
        res2<-cbind(res2,res$finger)
    }
    colnames(res1)<-colnames(res2)<-gsub(" ","_",gsub("^.*/","",files))
    return(list(finger=res2, viz=res1))

}

rename.extract<-function(path,vars){
    require(cytofCore)
    if (is.character(path)) fcs<-read.FCS(path) else fcs<-path
    fcs.exprs<-exprs(fcs)
    colnames(fcs.exprs)<-as.character(fcs@parameters@data$desc)
    fcs.exprs[,vars]
}

create.reference<-function(files,vars,n=20000){
    res<-NULL
    for (i in files){
        fcs<-read.FCS(i)
        fcs<-fcs[1:min(n,dim(fcs)[1])]
        print(i)
        res<-rbind(res,rename.extract(fcs,vars))
    }
    return(res)
}
