
# CDC Color
 cdc.blue <- "#005DAB"
 ncezid.orange <- "#C33A32"

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

png("mosquitoGlyph.png",res=600,height=11,width=8.5,units="in")

    drawMosquito(scale=1,col=cdc.blue)
    #drawMosquito(scale=1,col=ncezid.orange)

dev.off()

 library(magick)
 img <- image_read("mosquitoGlyph.png")
 img <- image_background(img, color = "none")
 img <- image_trim(img)

 library(hexSticker)

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

