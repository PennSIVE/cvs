dtboundary<-function(mask) {
  get.d.boundary.exact.balloon<-function(v,mask,d.max=30){
    if (mask[v[1],v[2],v[3]]==0) print('ERROR! - voxel outside of mask...')
    inf<-1000
    balloon.empty<-TRUE
    r<-1
    #expand balloon
    while(balloon.empty){
      balloon<-1-mask[(v[1]-r):(v[1]+r),(v[2]-r):(v[2]+r),(v[3]-r):(v[3]+r)]
      #If balloon had reached edge
      if (sum(balloon>0)){
        which.outside<-which(balloon>0,arr.ind=TRUE)
        d.out<-min(sqrt((which.outside[,1]-(r+1))^2+(which.outside[,2]-(r+1))^2+(which.outside[,3]-(r+1))^2))
        balloon.empty<-FALSE
      }
      else {
        if (r<=d.max) {
          r<-r+1
        }
        else {
          balloon.empty<-FALSE
          d.out<-inf
        }
      }
    }
    return(d.out)
  }
  which.mask.arrind<-which(mask>0,arr.ind=TRUE)
  #For each voxel in the mask
  min.d<-rep(0,dim(which.mask.arrind)[1])
  for (i in 1:(dim(which.mask.arrind)[1])) {
    #Get minimum distance to boundary
    min.d[i]<-get.d.boundary.exact.balloon(which.mask.arrind[i,],mask)
  }
  mask[mask>0]<-min.d
  return(mask)
}
frangifilter=function(image,mask,radius=1,color="dark",parallel=FALSE,cores=2,
                c3d=F,min.scale=0.5,max.scale=0.5){
  
  if(c3d==F){
    eigvals=hessian(image,mask,radius,parallel,cores)
    
    print("Calculating vesselness measure")
    l1=eigvals$eigval1
    l2=eigvals$eigval2
    l3=eigvals$eigval3
    l1=as.vector(l1[mask==1])
    l2=as.vector(l2[mask==1])
    l3=as.vector(l3[mask==1])
    rm(eigvals)
    
    al1=abs(l1)
    al2=abs(l2)
    al3=abs(l3)
    
    Ra=al2/al3
    Ra[!is.finite(Ra)]<-0
    Rb=al1/sqrt(al2*al3)
    Rb[!is.finite(Rb)]<-0
    
    S=sqrt(al1^2 + al2^2 + al3^2)
    A=2*(.5^2)
    B=2*(.5^2)
    C=2*(.5*max(S))^2
    
    rm(al1,al2,al3)
    
    eA=1-exp(-(Ra^2)/A)
    eB=exp(-(Rb^2)/B)
    eC=1-exp(-(S^2)/C)
    
    rm(Ra,Rb,S,A,B,C)
    
    vness=eA*eB*eC
    
    rm(eA,eB,eC)
    
    if(color=="dark"){
      vness[l2<0 | l3<0] = 0
      vness[!is.finite(vness)] = 0
    }else if(color=="bright"){
      vness[l2>0 | l3>0] = 0
      vness[!is.finite(vness)] = 0
    }
    
    image[mask==1]<-vness
    return(image)
  }else{
    tempinv=tempfile(pattern="file", tmpdir=tempdir(), fileext=".nii.gz")
    if(color=="dark"){
      writenii(-1*image,tempinv)
    }else{
      writenii(image,tempinv)
    }
    tempvein=tempfile(pattern="file", tmpdir=tempdir(), fileext=".nii.gz")
    system(paste0("c3d ",tempinv," -hessobj 1 ",min.scale," ",max.scale," -oo ",tempvein))
    
    veinmask=readnii(tempvein)
    return(veinmask)
  }
}
gradient=function(image,mask=NULL,which="all",radius=1){
  if(radius>=min(dim(image))){stop("Radius larger than smallest image dimension")}
  if(is.nifti(image)){
    if(which=="all"){
      dx=image
      dx@.Data[1:radius,,]=0
      dx@.Data[dim(image)[1]-radius+1,,]=0
      dx@.Data[((1+radius):(dim(image)[1]-radius)),,]=
        (image@.Data[((1+2*radius):dim(image)[1]),,]-
           image@.Data[(1:(dim(image)[1]-2*radius)),,])/(2*radius)
      dx@.Data[mask@.Data==0]<-0
      
      dy=image
      dy@.Data[,1:radius,]=0
      dy@.Data[,dim(image)[2]-radius+1,]=0
      dy@.Data[,((1+radius):(dim(image)[2]-radius)),]=
        (image@.Data[,((1+2*radius):dim(image)[2]),]-
           image@.Data[,(1:(dim(image)[2]-2*radius)),])/(2*radius)
      dy@.Data[mask@.Data==0]<-0
      
      dz=image
      dz@.Data[,,1:radius]=0
      dz@.Data[,,dim(image)[3]-radius+1]=0
      dz@.Data[,,((1+radius):(dim(image)[3]-radius))]=
        (image@.Data[,,((1+2*radius):dim(image)[3])]-
           image@.Data[,,(1:(dim(image)[3]-2*radius))])/(2*radius)
      dz@.Data[mask@.Data==0]<-0
      
      return(list(Dx=dx,Dy=dy,Dz=dz))
    }else if(which=="x"){
      dx=image
      dx@.Data[1:radius,,]=0
      dx@.Data[dim(image)[1]-radius+1,,]=0
      dx@.Data[((1+radius):(dim(image)[1]-radius)),,]=
        (image@.Data[((1+2*radius):dim(image)[1]),,]-
           image@.Data[(1:(dim(image)[1]-2*radius)),,])/(2*radius)
      dx@.Data[mask@.Data==0]<-0
      return(dx)
    }else if(which=="y"){
      dy=image
      dy@.Data[,1:radius,]=0
      dy@.Data[,dim(image)[2]-radius+1,]=0
      dy@.Data[,((1+radius):(dim(image)[2]-radius)),]=
        (image@.Data[,((1+2*radius):dim(image)[2]),]-
           image@.Data[,(1:(dim(image)[2]-2*radius)),])/(2*radius)
      dy@.Data[mask@.Data==0]<-0
      return(dy)
    }else if(which=="z"){
      dz=image
      dz@.Data[,,1:radius]=0
      dz@.Data[,,dim(image)[3]-radius+1]=0
      dz@.Data[,,((1+radius):(dim(image)[3]-radius))]=
        (image@.Data[,,((1+2*radius):dim(image)[3])]-
           image@.Data[,,(1:(dim(image)[3]-2*radius))])/(2*radius)
      dz@.Data[mask@.Data==0]<-0
      return(dz)
    }
  }else if(is.array(image)){
    if(which=="all"){
      dx=image
      dx[1:radius,,]=0
      dx[dim(image)[1]-radius+1,,]=0
      dx[((1+radius):(dim(image)[1]-radius)),,]=
        (image[((1+2*radius):dim(image)[1]),,]-
           image[(1:(dim(image)[1]-2*radius)),,])/(2*radius)
      dx[mask==0]<-0
      
      dy=image
      dy[,1:radius,]=0
      dy[,dim(image)[2]-radius+1,]=0
      dy[,((1+radius):(dim(image)[2]-radius)),]=
        (image[,((1+2*radius):dim(image)[2]),]-
           image[,(1:(dim(image)[2]-2*radius)),])/(2*radius)
      dy[mask==0]<-0
      
      dz=image
      dz[,,1:radius]=0
      dz[,,dim(image)[3]-radius+1]=0
      dz[,,((1+radius):(dim(image)[3]-radius))]=
        (image[,,((1+2*radius):dim(image)[3])]-
           image[,,(1:(dim(image)[3]-2*radius))])/(2*radius)
      dz[mask==0]<-0
      
      return(list(Dx=dx,Dy=dy,Dz=dz))
    }else if(which=="x"){
      dx=image
      dx[1:radius,,]=0
      dx[dim(image)[1]-radius+1,,]=0
      dx[((1+radius):(dim(image)[1]-radius)),,]=
        (image[((1+2*radius):dim(image)[1]),,]-
           image[(1:(dim(image)[1]-2*radius)),,])/(2*radius)
      dx[mask==0]<-0
      return(dx)
    }else if(which=="y"){
      dy=image
      dy[,1:radius,]=0
      dy[,dim(image)[2]-radius+1,]=0
      dy[,((1+radius):(dim(image)[2]-radius)),]=
        (image[,((1+2*radius):dim(image)[2]),]-
           image[,(1:(dim(image)[2]-2*radius)),])/(2*radius)
      dy[mask==0]<-0
      return(dy)
    }else if(which=="z"){
      dz=image
      dz[,,1:radius]=0
      dz[,,dim(image)[3]-radius+1]=0
      dz[,,((1+radius):(dim(image)[3]-radius))]=
        (image[,,((1+2*radius):dim(image)[3])]-
           image[,,(1:(dim(image)[3]-2*radius))])/(2*radius)
      dz[mask==0]<-0
      return(dz)
    }
  }else{
    print("Image must be array or NifTI")
  }
}
hessian=function(image,mask,radius=1,parallel=FALSE,cores=2){
  
  print("Getting derivatives")
  grads=gradient(image,which="all",radius=radius)
  gx=grads$Dx
  gy=grads$Dy
  gz=grads$Dz
  rm(grads)
  
  gradsx=gradient(gx,which="all",radius=radius)
  gxx=gradsx$Dx
  gxy=gradsx$Dy
  gxz=gradsx$Dz
  rm(gx,gradsx)
  
  gradsy=gradient(gy,which="all",radius=radius)
  gyx=gradsy$Dx
  gyy=gradsy$Dy
  gyz=gradsy$Dz
  rm(gy,gradsy)
  
  gradsz=gradient(gz,which="all",radius=radius)
  gzx=gradsz$Dx
  gzy=gradsz$Dy
  gzz=gradsz$Dz
  rm(gz,gradsz)
  
  print("Creating hessian matrices")
  bigmat=cbind(as.vector(gxx[mask==1]),as.vector(gxy[mask==1]),as.vector(gxz[mask==1]),
               as.vector(gyx[mask==1]),as.vector(gyy[mask==1]),as.vector(gyz[mask==1]),
               as.vector(gzx[mask==1]),as.vector(gzy[mask==1]),as.vector(gzz[mask==1]))
  
  rm(gxx,gxy,gxz,gyx,gyy,gyz,gzx,gzy,gzz)
  
  biglist=split(bigmat,row(bigmat))
  biglist=lapply(biglist,matrix,nrow=3,byrow=T)
  
  rm(bigmat)
  
  getevals=function(matrix){
    thiseig=eigen(matrix)$values
    sort=order(abs(thiseig))
    return(thiseig[sort])
  }
  
  print("Calculating eigenvalues")
  if(parallel==TRUE){
    result=matrix(unlist(pbmclapply(biglist,getevals,mc.cores=cores)),
                  ncol=3,byrow=T)
  }else if(parallel==FALSE){
    result=matrix(unlist(lapply(biglist,getevals)),ncol=3,byrow=T)
  }
  e1=mask
  e1[mask==1]<-result[,1]
  e2=mask
  e2[mask==1]<-result[,2]
  e3=mask
  e3[mask==1]<-result[,3]
  
  return(list(eigval1=e1,eigval2=e2,eigval3=e3))
}
labelreg=function(fullimage,labelimage,fixedimage,typeofTransform="Rigid",
                  interpolator="lanczosWindowedSinc"){
  imtofix=registration(filename=fullimage,template.file=fixedimage,
                       typeofTransform=typeofTransform,remove.warp=FALSE,
                       outprefix="fun")
  labtofix<-antsApplyTransforms(fixed=oro2ants(fixedimage),moving=oro2ants(labelimage),
                                transformlist=imtofix$fwdtransforms,
                                interpolator=interpolator)
  return(list(image_reg=ants2oro(imtofix$outfile),label_reg=ants2oro(labtofix)))
}
lesioncenters=function(probmap,binmap,minCenterSize=10,radius=1,
                       parallel=F,cores=2){
  scale=ceiling((1/mean(probmap@pixdim[2:4]))^3)
  phes=hessian(probmap,mask=binmap,radius,parallel,cores)
  clusmap=ants2oro(labelClusters(oro2ants(binmap),minClusterSize=20*scale))
  
  les=clusmap
  les[les!=0]<-1
  les[phes$eigval1>0 | phes$eigval2>0 | phes$eigval3>0]<-0
  les=ants2oro(labelClusters(oro2ants(les),minClusterSize=minCenterSize*scale))
  
  return(list(lesioncenters=les,lesioncount=max(les)))
}
getnulldist=function(x,centsub,coords,frangsub){
  samp=sample(1:length(centsub))
  centsamp=centsub[samp]
  coordsamp=coords[samp,]
  sampprod=frangsub*centsamp
  return(sum(sampprod))
}
getcands=function(x,lables,nottiss){
  if(sum(nottiss[lables==x])>0){
    return(0)
  }else{
    return(1)
  }
}
