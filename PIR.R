
# CDC Color
 cdc.blue <- "#005DAB"
 ncezid.orange <- "#C33A32"
 
# Michigan Colors
 mi.blue <- "#00274C"
 mi.maze <- "#FFCB05"
     
"polygonBand" <-  function(x,lower,upper,f1=NULL,f2=NULL,fill=TRUE,col="lightgray",density,border=FALSE,...){
     if(!is.null(f1)) lower <- f1(x)
     if(!is.null(f2)) upper <- f2(x)
     if(fill) density <- -1
     if(!fill & missing(density)){
         density <- 10
         border <- TRUE
     }
     polygon(c(x,rev(x)), c(lower, rev(upper)), density=density, border=border,col=col,...)
}

"ellipse" <- function(a=1,b=1,h=0,k=0,from=0,to=2*pi,npts=1000){
  theta <- seq(from,to,length=npts)
  x <- cos(theta)
  y <- sin(theta)
  x <- x/a-h
  y <- y/b-k
  list(x=x,y=y)
}

"lemniscate" <- function(a=1,b=1,from=pi/2,to=3*pi/2,orientation=c("vertical","horizontal"),npts=1000){
   orientation <- match.arg(orientation)
   theta <- seq(from,to,length=npts)
   switch(orientation,
        vertical =  ans <- list(y=b*a*cos(theta)/(1+sin(theta)^2),x=a*sin(theta)*cos(theta)/(1+sin(theta)^2)),
        horizontal =  ans <- list(x=a*cos(theta)/(1+sin(theta)^2),y=b*a*sin(theta)*cos(theta)/(1+sin(theta)^2))
   )
   ans
}

#"lemniscate.f" <- function(x, a=1, b=1,orientation=c("vertical","horizontal")){
#    orientation <- match.arg(orientation)
#    x <- sort(x)
#    n <- length(x)
#    from <- x[1]
#    to <- x[n]
#    ans <- vector(length=n)
#    for(i in 1:n)
#        ans[i] <-  lemniscate(a,b,x[i],x[i],orientation,npts=1)$y
#    ans
#}

# figure built with these xlim=c(-5,5), ylim=c(-20,20)
"drawMosquito" <- function(scale=1,bg=c("transparent","white"),lwd=3,col="black"){
    scale <- 1/scale # so the interpretation at the call is bigger scale = bigger image
    xlim <- scale*c(-5,5)
    ylim <- scale*c(-20,20)
    #BG <- match.arg(bg)
    par(bg=bg[1], pty="s")
    
    plot(lemniscate(1,10,pi/2,3*pi/2,orient="v"),lwd=lwd,xlim=xlim,ylim=ylim,axes=FALSE,xlab="",ylab="",type="l",col=col)
        lmn <- lemniscate(1,10,pi/2,3*pi/2,orient="v")
        polygonBand(lmn$x,lmn$x,lmn$y,col=col,border=FALSE)
    
      lines(lemniscate(1,4,from=0,to=pi/2,orient="v"),lwd=lwd,col=col)
        lmn <- lemniscate(1,4,from=0,to=pi/2,orient="v")
        polygonBand(lmn$x,lmn$x,lmn$y,col=col,border=FALSE)
      
      lines(lemniscate(1,4,from=3*pi/2,to=2*pi,orient="v"),lwd=lwd,col=col)
        lmn <- lemniscate(1,4,from=3*pi/2,to=2*pi,orient="v")
        polygonBand(lmn$x,lmn$x,lmn$y,col=col,border=FALSE)
      
      lines(ellipse(a=5,b=1.5,h=0,k=-4.7),lwd=lwd,col=col)
        ell <- ellipse(a=5,b=1.5,h=0,k=-4.7)
        polygonBand(ell$x,ell$x,ell$y,col=col,border=FALSE)
        
      segments(0,0,1,-5,lwd=lwd,col=col)
        segments(1,-5,0.5,-11,lwd=lwd,col=col)
        segments(0.5,-11,1.2,-18,lwd=lwd,col=col)
      segments(0,0,-1,-5,lwd=lwd,col=col)
        segments(-1,-5,-0.5,-11,lwd=lwd,col=col)
        segments(-0.5,-11,-1.2,-18,lwd=lwd,col=col)
      segments(0,0,1.4,-4,lwd=lwd,col=col)
        segments(1.4,-4,1,-10,lwd=lwd,col=col)
        segments(1,-10,1.3,-15,lwd=lwd,col=col)
      segments(0,0,-1.4,-4,lwd=lwd,col=col)
        segments(-1.4,-4,-1,-10,lwd=lwd,col=col)
        segments(-1,-10,-1.3,-15,lwd=lwd,col=col)
      segments(0.35,2.8,1,4,lwd=lwd,col=col)
        segments(1,4,1.15,8,lwd=lwd,col=col)
        segments(1.15,8,2,11,lwd=lwd,col=col)
      segments(-0.35,2.8,-1,4,lwd=lwd,col=col)
        segments(-1,4,-1.15,8,lwd=lwd,col=col)
        segments(-1.15,8,-2,11,lwd=lwd,col=col)
      # proboscis & antennas
      segments(0,5.4,0,9.5,lwd=lwd,col=col)  
        segments(0,5.4,0.45,8,lwd=lwd,col=col)
        segments(0,5.4,-0.45,8,lwd=lwd,col=col)
      # wings
      segments(0.2,0.9,3,2.4,lwd=lwd,col=col)
        segments(0,0,2,-0.5,lwd=lwd,col=col)
        segments(2,-0.5,3.1,1.65,lwd=lwd,col=col)
        ell.right <- ellipse(a=5,b=2.5,h=-3,k=-2,from=-pi/3,to=pi/2) # wingtip
        lines(ell.right,lwd=lwd,col=col)
     segments(-0.2,0.9,-3,2.4,lwd=lwd,col=col)
        segments(0,0,-2,-0.5,lwd=lwd,col=col)
        segments(-2,-0.5,-3.1,1.65,lwd=lwd,col=col)
        ell.left <- ell.right
        ell.left$x <- -ell.right$x
        lines(ell.left,lwd=lwd,col=col) # wingtip
    invisible(scale)
}   

#tiff("mosquitoGlyph.tif",res=600,height=11,width=8.5,units="in")
png("mosquitoGlyph.png",res=600,height=11,width=8.5,units="in")
#pdf("mosquitoGlyph.pdf",height=11,width=8.5)

    drawMosquito(scale=1,col=cdc.blue)
    #drawMosquito(scale=1,col=ncezid.orange)

dev.off()   
 
 library(magick)
 img <- image_read("mosquitoGlyph.png")
 img <- image_background(img, color = "none")
 img <- image_trim(img)
# img <- image_resize(img, geometry_size_percent(width=100))

 library(hexSticker)
 # plot(sticker(img,package="PooledInfRate",
 #              s_width=0.9,s_height=1.25,
 #         h_color="dodgerblue4",h_fill="gold1",
 #         p_size=25,p_y=1.5, p_color="dodgerblue4",
 #         s_x=1,s_y=0.8,
 #         filename="PIR.png",
 #         asp = TRUE,
 #         dpi=600) )
 # 
 
 plot(sticker(img,package="PooledInfRate",
              s_width=0.9,s_height=1.25,
              h_color=mi.blue,h_fill=mi.maze,
              p_size=25,p_y=1.5, p_color=mi.blue,
              s_x=1,s_y=0.8,
              filename="PIR.png",
              asp = TRUE,
              dpi=600) )
 
 # CDC Colors
 plot(sticker(img,package="PooledInfRate",
              s_width=0.9,s_height=1.25,
              h_color=cdc.blue,h_fill="white",
              p_size=25,p_y=1.5, p_color=cdc.blue,
              s_x=1,s_y=0.8,
              h_size=1.2, # 1.2 is default
              filename="PIR.png",
              asp = TRUE,
              dpi=600,
              white_around_sticker = FALSE) )
 
 # CDC blue with NCEZID orange
 plot(sticker(img,package="PooledInfRate",
              s_width=0.9,s_height=1.25,
              h_color=ncezid.orange,h_fill=cdc.blue,
              p_size=25,p_y=1.5, p_color=ncezid.orange,
              s_x=1,s_y=0.8,
              filename="PIRBlueOrange.png",
              asp = TRUE,
              dpi=600) )
#  
# "hex" <- function(h=1.73,x.scale=1,y.scale=4,center=c(0,0),...){
#     scale <- c(x.scale, y.scale)
#     v1 <- center + scale * c(h/2, h/2) 
#     v2 <- center + scale * c(0, h/2 + h/2 ) 
#     v3 <- center + scale * c(-h/2,h/2) 
#     v4 <- center + scale * c(-h/2,-h/2) 
#     v5 <- center + scale * c(0,-(h/2 + h/2)) 
#     v6 <- center + scale * c(h/2,-h/2) 
#     segments(v1[1],v1[2], v2[1], v2[2],...)
#     segments(v2[1],v2[2], v3[1], v3[2],...)
#     segments(v3[1],v3[2], v4[1], v4[2],...)
#     segments(v4[1],v4[2], v5[1], v5[2],...)
#     segments(v5[1],v5[2], v6[1], v6[2],...)
#     segments(v6[1],v6[2], v1[1], v1[2],...)
#     invisible(h)
# }  

 
 
 
# scale <- 1
# lwd <- 3
# xlim <- scale*c(-5,5)
# ylim <- scale*c(-20,20)
# 
# plot(lemniscate(1,10,pi/2,3*pi/2,orient="v"),lwd=lwd,xlim=xlim,ylim=ylim,axes=FALSE,xlab="",ylab="",type="l")
# lmn <- lemniscate(1,10,pi/2,3*pi/2,orient="v")
# polygonBand(lmn$x,lower=lmn$x,upper=lmn$y,col='black')
 



"addTick" <- function(a=2,b=5,x.shift=0,y.shift=0,scale=1,bg=c("transparent","white"),lwd=3,col="black",add=TRUE){
    scale <- 1/scale # so the interpretation at the call is bigger scale = bigger image
    xlim <- scale*c(-5,5)
    ylim <- scale*c(-20,20)
    #BG <- match.arg(bg)
    par(bg=bg[1])
    
    lmn <- lemniscate(a,b,from=pi/2,to=3*pi/2,orientation="horizontal")
    lmn$x <- lmn$x + x.shift 
    lmn$y <- lmn$y + y.shift 
    if(add){
        lines(lmn,col=col)
    } else {
        plot(lmn,lwd=lwd,xlim=xlim,ylim=ylim,axes=FALSE,xlab="",ylab="",type="l",col=col)
    }
    polygonBand(lmn$x,lmn$x,lmn$y,col=col,border=FALSE)

    points(0+x.shift,0+y.shift,pch=16,col=cdc.blue,cex=1.5)
    
    segments(-0.25+x.shift,-1+y.shift,-0.2+x.shift,-2.5+y.shift,lwd=3,col=cdc.blue)
    segments(-0.2+x.shift,-2.5+y.shift,0.5+x.shift,-3.5+y.shift,lwd=3,col=cdc.blue)
    segments(-0.25+x.shift,1+y.shift,-0.2+x.shift,2.5+y.shift,lwd=3,col=cdc.blue)
    segments(-0.2+x.shift,2.5+y.shift,0.5+x.shift,3.5+y.shift,lwd=3,col=cdc.blue)
    
    segments(-0.5+x.shift,-2+y.shift,-0.3+x.shift,-5.5+y.shift,lwd=3,col=cdc.blue)
    segments(-0.5+x.shift,2+y.shift,-0.3+x.shift,5.5+y.shift,lwd=3,col=cdc.blue)
    
    segments(-0.8+x.shift,-3.1+y.shift,-1.5+x.shift,-6.5+y.shift,lwd=3,col=cdc.blue)
    segments(-0.8+x.shift,3.1+y.shift,-1.5+x.shift,6.5+y.shift,lwd=3,col=cdc.blue)
    
    segments(-1.3+x.shift,-3.6+y.shift,-2.5+x.shift,-6.25+y.shift,lwd=3,col=cdc.blue)
    segments(-1.3+x.shift,3.6+y.shift,-2.5+x.shift,6.25+y.shift,lwd=3,col=cdc.blue)
    
    
    
    
    invisible(scale)
}   

"addFlea" <- function(a=2,b=5,x.shift=0,y.shift=0,scale=1,bg=c("transparent","white"),lwd=3,col="black",add=TRUE){
    scale <- 1/scale # so the interpretation at the call is bigger scale = bigger image
    xlim <- scale*c(-5,5)
    ylim <- scale*c(-20,20)
    #BG <- match.arg(bg)
    par(bg=bg[1])
    
    # front side
    lmn <- lemniscate(a,b,from=-pi/2,to=pi/2,orientation="horizontal")
    lmn$x <- lmn$x + x.shift 
    lmn$y <- lmn$y + y.shift 
    if(add){
        lines(lmn,col=col)
    } else {
        plot(lmn,lwd=lwd,xlim=xlim,ylim=ylim,axes=FALSE,xlab="",ylab="",type="l",col=col)
    }
    polygonBand(lmn$x,lmn$x,lmn$y,col=col,border=FALSE)
    
    # back side
    lmn <- lemniscate(a,b,from=pi/2,to=3*pi/2,orientation="horizontal")
    lmn$x <- lmn$x + x.shift + 1.2*a
    lmn$y <- lmn$y + y.shift 
    lines(lmn,col=col)
    polygonBand(lmn$x,lmn$x,lmn$y,col=col,border=FALSE)
   
    # head
    #points(0.1+x.shift,-0.1+y.shift,pch=16,col=col,cex=1.5) 
    points(0.2+x.shift,0+y.shift,pch=16,col=col,cex=1.25) 
    points(0.1+x.shift,-0.1+y.shift,pch=16,col=col,cex=1.5) 
    points(0.075+x.shift,-0.3+y.shift,pch=16,col=col,cex=1.5) 
    points(-0.05+x.shift,-0.6+y.shift,pch=16,col=col,cex=1.4) 
    points(-0.03+x.shift,-0.5+y.shift,pch=16,col=col,cex=1.5) 
    points(-0.02+x.shift,-0.45+y.shift,pch=16,col=col,cex=1.5) 
    points(0.25+x.shift,0.1+y.shift,pch=16,col=col,cex=1.6) 
    #points(0.2+x.shift,0.9+y.shift,pch=16,col=col,cex=0.75)
    
    # legs
    # front
    segments(0.15+x.shift,-0.7+y.shift,0.25+x.shift,-2.5+y.shift,col=col,lwd=3)
    segments(0.25+x.shift,-2.5+y.shift,0.1+x.shift,-4.5+y.shift,col=col,lwd=3)
    
    segments(0.15+x.shift,-0.7+y.shift,0.4+x.shift,-3+y.shift,col=col,lwd=3)
    segments(0.4+x.shift,-3+y.shift,0.3+x.shift,-4.5+y.shift,col=col,lwd=3)
    
    
    # middle 
    segments(0.45+x.shift,-1.6+y.shift,0.6+x.shift,-3.25+y.shift,col=col,lwd=3)
    segments(0.6+x.shift,-3.25+y.shift,0.7+x.shift,-3.4+y.shift,col=col,lwd=3)
    segments(0.7+x.shift,-3.4+y.shift,0.6+x.shift,-5.5+y.shift,col=col,lwd=3)
    
    segments(0.55+x.shift,-1.8+y.shift,0.8+x.shift,-3.25+y.shift,col=col,lwd=3)
    segments(0.8+x.shift,-3.25+y.shift,0.9+x.shift,-3.4+y.shift,col=col,lwd=3)
    segments(0.9+x.shift,-3.4+y.shift,0.8+x.shift,-5.5+y.shift,col=col,lwd=3)
    
    
    # back 
    segments(0.95+x.shift,-2.7+y.shift,1.2+x.shift,-4+y.shift,col=col,lwd=3)
    segments(1.2+x.shift,-4+y.shift,1.4+x.shift,-4+y.shift,col=col,lwd=3)
    segments(1.4+x.shift,-4+y.shift,2.4+x.shift,-5.5+y.shift,col=col,lwd=3)
    
    segments(0.95+x.shift,-2.7+y.shift,1.1+x.shift,-4.3+y.shift,col=col,lwd=3)
    segments(1.1+x.shift,-4.3+y.shift,1.2+x.shift,-4.4+y.shift,col=col,lwd=3)
    segments(1.2+x.shift,-4.4+y.shift,2+x.shift,-5.8+y.shift,col=col,lwd=3)
    
    
    
    
    invisible(scale)
}   



#tiff("mosquitoGlyph.tif",res=600,height=11,width=8.5,units="in")
png("PIRGlyph.png",res=600,height=11,width=8.5,units="in")
    drawMosquito(scale=0.75,col=cdc.blue,bg="white")
    addTick(a=2,b=5,x.shift=-3,y.shift=17,col=cdc.blue)
    addFlea(a=2,b=4,x.shift=3,y.shift=17,col=cdc.blue,bg="white")
dev.off()   

library(magick)
img <- image_read("PIRGlyph.png")
img <- image_background(img, color = "none")
img <- image_trim(img)
# img <- image_resize(img, geometry_size_percent(width=100))


# CDC Colors
plot(sticker(img,package="PooledInfRate",
             s_width=0.9,s_height=1.25,
             h_color=cdc.blue,h_fill="white",
             p_size=25,p_y=1.5, p_color=cdc.blue,
             s_x=1,s_y=0.8,
             h_size=1.2, # 1.2 is default
             filename="PIRFinal.png",
             asp = TRUE,
             dpi=600,
             white_around_sticker = FALSE) )


