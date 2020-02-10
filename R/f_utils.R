
fill_panel_plot <- function(model_res, raw_data, predictions, plot_pars) {
  
  par(plot_pars)
  
  for (funct in model_res$funct_trait) {
    
    tab   <- raw_data %>% filter(Funct == funct) %>% filter(value != -1)    
    mod   <- model_res[model_res$funct_trait == funct, ]
    pred <- predictions[[funct]]
    
    # plot data
    
    ylim  <- setrge(tab$value, 25)
    xlim  <- c(.5, num_q + 0.5)
    title <- paste(strsplit(funct, "_")[[1]], collapse = "\n")
    plot(tab$value ~ tab$Rank, type = "n", xaxt = 'n', main = "", ylim = ylim, 
         xlim = xlim, xlab = "", ylab = "")
    text(5, max(ylim) - 0.15 * diff(ylim), labels = title, cex = 1.3,
         col = colors_category[mod$funct_category], font = 2)
    pt_cex <- 0.5
    segments(xlim[1],1,xlim[2], 1, lty = 2)
    points(tab$value ~ tab$Rank, pch = 21, cex = pt_cex, col = "#75757530", bg  = "#75757530")
    
    # plot model
    lapply(names(pred), function(m) {
      if (m == mod$best_model) {
        lines(1:10, pred[[m]], lwd = 3, col = 'black')
      } else {
        lines(1:10, pred[[m]], lwd = 1, col = 'black', lty = 3)
      }
    })
 
    # legend('top', bty = 'n', cex = 1, text.font = 4,
    #        legend = paste('R2=', round_text(mod[, "r2m"], 2)))
    axis(1, 1:10, labels = 1:10, font = 1, cex = 0.1)
  }
  mtext(1, text = "Abundance-occupancy quantile", outer = T, font = 2, line = 1)
  mtext(2, text = "Weight of the function in quantile", outer = T, las = 0, font = 2, line = 1)
}


reformat_as_df <- function(input_list, new_var_name = NULL) {
  out <- lapply(names(input_list), function(x) {
    mutate(input_list[[x]], new_var = as.character(rep(x, nrow(input_list[[x]]))))
  }) %>% bind_rows() 
  names(out)[names(out) == "new_var"] <- new_var_name
  out
}


simpleCap <- function(x) {
  s <- strsplit(x, " ")[[1]]
  paste(toupper(substring(s, 1, 1)), substring(s, 2),
        sep = "", collapse = " ")
}


plotTokeshi <- function(x, colPch='grey', colLine='black', xlim, ylim, cexPch=1, xlab='', ylab='', main_text=NULL){
  par(las=1, mgp=c(3,0.5,0), tcl=-0.3)
  X <- x$mod$model$x
  Y <- x$mod$y
  xx <- seq(min(X), max(X), len=101)
  pre <- predict(x$mod, newdata=list(x = xx), se=TRUE)
  g <- x$mod$family$linkinv
  fv <- g(pre$fit)
  hi <- g(pre$fit + 2*pre$se)
  lo <- g(pre$fit - 2*pre$se)
  plot(X, Y, pch=21, col=colPch,xlim=xlim, ylim=ylim, bg=colPch, xlab=xlab, ylab=ylab, xaxt='n', yaxs='i', cex=cexPch, main = main_text)
  
  matlines(xx, cbind(fv, hi, lo), lty=c(1, 2, 2), lwd=c(2, 1, 1), col=colLine )
}


load_libraries <- function(x){
  eval(parse(text = paste0( 
    "if (!require(", x, ")) { 
      install.packages(", x, ") 
      require(", x, ")
    } else {
      require(", x, ")
    }"
  )))
}

# ROUND_TEXT: function to round a number (x) and get a character string with defined number of decimals (digits)

round_text<-function(x,digits=1) {
  
  round_x<-round(x,digits)
  length_round_x<-nchar(as.character(abs(round_x)))
  length_int_round_x<-nchar(as.character(trunc(abs(round_x))))
  
  # round x is an integer
  if( length_round_x==length_int_round_x) {
    rep0<-""
    for (k in 1:digits)
      rep0<-paste(rep0,"0",sep="")
    res<-paste( trunc(round_x),".", rep0,sep="") } # end of if no decimals
  
  # else
  if( length_round_x>length_int_round_x) {
    nb_dec<-length_round_x - length_int_round_x -1
    
    if (nb_dec==digits) res<-as.character(round_x) # ok
    
    if (nb_dec<digits) {
      rep0<-""
      for (k in 1:(digits-nb_dec) )
        rep0<-paste(rep0,"0",sep="")
      res<-paste( round_x, rep0 ,sep="") } # end of if 0 missing
  }# end of if decimals
  
  return(res)
  
} # end of function round_text

###########################################################################################################
###########################################################################################################
# PVALUE_STAR: function to convert a pvalue (p) in classical text notation: ns, *, **, ***
# by default if p>0.05, "ns" is returned but it can be changed with argument "ns" (e.g. to " ")
pvalue_star<-function(p,ns="ns") {
  res<-ns
  if(p<0.05) res<-"*"
  if(p<0.01) res<-"**"
  if(p<0.001) res<-"***"
  
  return(res)	
  
} # end of function pvalue_star




######################################################################################################################################
# NICEPLOT : function to build a graphical window with range, labels and title of X and Y axes as specified
# default values are for a single plot saved at resolution=150dpi, width=600pixels, height=600pixels
# (i.e. also works with res=300, width=1200, heigth=1200)
#
# "limX" and "limY": range of X and Y axes (min, max), see function "setrge" below to set limits given data
#
# "marg" : external margins width (in lines); bottom, left, top, right
#
# "tick" : length of tick marks, negative values for outside tick marks, positive for inner ones
#
# "nlab" : maximum number of ticks on each axis
# "labX" and "labY" : position of the labels on the axis
#						if "lab" argument remains empty, labels are set by default
#           			if "lab" argument is NA, no ticks and labels are added
# "nmlabX" and "nmlabY" : label names, default is nmlab=lab
# "lasX" and "lasY" : orientation of labels (1= horizontal, 3=vertical) 
# "lineX" and "lineY" : distance (in lines) between labels and axis
# "cexX" and "cexY" : size of labels characters
#
# "nmX" and "nmY": title of X and Y axes 
# lineXt" and "lineYt" : distance between an axis and its title
# "cexXt" and "cexYt" : size of axis title characters
#
#######################################################################################################################################

niceplot<-function(limX=c(0,1),limY=c(0,1), marg=c(4.5,4.5,2,2),
                      tick=-0.4,   nlab=7, labX=c(),labY=c(),    nmlabX=c(),nmlabY=c(),
                      lasX=1,lasY=1,    lineX=-0.1, lineY=0.1,   cexX=0.9, cexY=0.9,
                      nmX="X",nmY="Y",   lineXt=lineX+2, lineYt=lineY+2.5,   cexXt=1, cexYt=1   ) {
                      
par(mar=marg)   # margins
plot(limX,limY,type="n",axes=F,xaxt="n",yaxt="n",xlab="",ylab="",xlim=limX,ylim=limY) # window
rect(limX[1],limY[1],limX[2],limY[2])   # border
mtext(side=1,nmX,cex=cexXt,line=lineXt,font=2) # X title  
mtext(side=2,nmY,cex=cexYt,line=lineYt,font=2) # Y title 

# labels for X axis
if (is.na(sum(labX))==F) {
labx<-pretty(limX,n=nlab) ; labx<-labx[which(labx>=min(limX) & labx<=max(limX))] ; nmlabx<-labx  # default
if (length(labX)>0) { labx<-labX ; nmlabx<-nmlabX }                                                   # customized
axis(side=1, at=labx, labels=F, tcl=tick, pos=limY[1])  # ticks
mtext(side=1, nmlabx, at=labx, line=lineX, cex=cexX, las=lasX) # labels
                        } # end of if labels
# labels for Y axis
if (is.na(sum(labY))==F) {
laby<-pretty(limY,n=nlab) ; laby<-laby[which(laby>=min(limY) & laby<=max(limY))] ; nmlaby<-laby  # default
if (length(labY)>0) { laby<-labY ; nmlaby<-nmlabY }                                                   # customized
axis(side=2, at=laby, labels=F, tcl=tick, pos=limX[1]) # ticks 
mtext(side=2, nmlaby, at=laby, line=lineY, cex=cexY, las=lasY) # labels
 
                          } # end of if labels
 
} # end of function

##########################################################################################################################
##########################################################################################################################
# SETRGE : function to compute limits of an axis given data
# inputs: "x" a continuous variable (i.e. at least two values) and "p" a coefficent of extension (in %)
# output: limits of the axis = observed range extended by p (default 5%)
 
setrge<-function(x,p=5) {

# additional range 
plus<-(p/100)*(max(x,na.rm=T)-min(x,na.rm=T))

# lower limit
low<-min(x,na.rm=T)- plus

# upper limit
up<-max(x,na.rm=T)+ plus

return(c(low,up))
} # end of function

###########################################################################################################
###########################################################################################################
# SE : function to compute standard error of the mean for a continuous variable X
se<-function(x) { sd(x,na.rm=T)/sqrt(length(na.omit(x))) }

################################################################################################################################
################################################################################################################################
# MEANSEY : function to add points and bars representing mean and associated error (e.g. standard deviation or standard error)
#           for several modalities (X axis) and one variable (Y axis)
# inputs : 
#   x: positions on X axis at which points representing mean have to be plotted
#   meany and sey : vector of same length than X with mean and associated error values to be plotted 
#   pchp, colp, bgp, cexp: shape, border color, background color and size of points representing mean (single value or vector)
#   colb : color for error bars (single value or vector)
#   lgt : width of error bars ticks, in percentage of X range
meansey<-function(x,meany,sey,pchp=22,colp="black",bgp="black",cexp=1.5,colb="black",lgt=0.01) {
lgx<-lgt*(max(x,na.rm=T)-min(x,na.rm=T)) # tick width
segments(x,meany-sey,x,meany+sey,col=colb) # error bar
segments(x-lgx,meany-sey,x+lgx,meany-sey,col=colb) # top tick
segments(x-lgx,meany+sey,x+lgx,meany+sey,col=colb) # bottom tick
points(x,meany,pch=pchp,col=colp,bg=bgp,cex=cexp) # mean
} # end of meansey

################################################################################################################################
################################################################################################################################
# POINTBAR : generic function to add points and vertical bars for several modalities (x axis) and one variable (y axis)
#				e.g. median and 1st and 3rd quartile or mean and 95% confidence-interval      
# inputs : 
#   x: positions on X axis at which points have to be plotted
#   pointy : vector of same length than X with point position on Y axis (e.g. median)
#   boty : vector of same length than X with bottom limit of bar on Y axis, e.g. first quartile values 
#   topy : vector of same length than X with top limit of bar on Y axis, e.g. third quartlie values
#   pchp, colp, bgp, cexp: shape, border color, background color and size of points (single value or vector)
#   colb : color for error bars (single value or vector)
#   lgt : width of error bars ticks, in percentage of X range
pointbar<-function(x,pointy,boty,topy,pchp=22,colp="black",bgp="black",cexp=1.5,colb="black",lgt=0.01) {
lgx<-lgt*(max(x,na.rm=T)-min(x,na.rm=T)) # tick width
segments(x,boty,x,topy,col=colb) # error bar
segments(x-lgx,boty,x+lgx,boty,col=colb) # top tick
segments(x-lgx,topy,x+lgx,topy,col=colb) # bottom tick
points(x,pointy,pch=pchp,col=colp,bg=bgp,cex=cexp) # mean
} # end of pointbar

################################################################################################################################
################################################################################################################################
# MEANSEXY : function to add points and bars representing mean and associated error for two continuous variables for several objects
# inputs : 
#   meanxy : a matrix (k,2) with mean values for the k objects for the two variables (X=column 1, Y=column 2)
#   sexy : a 2-columns matrix (same rows than meanxy) with associated error values for the two variables (X=column 1, Y=column 2)
#   pchp, colp, bgp, cexp: shape, border color, background color and size of points representing means (single value or vector)
#   colb : color for error bars (single value or vector)
#   lgt : width of error bars ticks, in percentage of mean+-se ranges on the two axes

meansexy<-function(meanxy,sexy,pchp=22,colp="black",bgp="black",cexp=1.5,colb="black",lgt=0.01, lwd_seg=1) {
lgx<-lgt*(max(meanxy[,1]+sexy[,1],na.rm=T)-min(meanxy[,1]-sexy[,1],na.rm=T)) # x tick width
lgy<-lgt*(max(meanxy[,2]+sexy[,2],na.rm=T)-min(meanxy[,2]-sexy[,2],na.rm=T)) # y tick width

segments(meanxy[,1]-sexy[,1],meanxy[,2],meanxy[,1]+sexy[,1],meanxy[,2],col=colb, lwd=lwd_seg) # x error bar
segments(meanxy[,1],meanxy[,2]-sexy[,2],meanxy[,1],meanxy[,2]+sexy[,2],col=colb, lwd=lwd_seg) # y error bar
segments(meanxy[,1]-sexy[,1],meanxy[,2]-lgy,meanxy[,1]-sexy[,1],meanxy[,2]+lgy,col=colb, lwd=lwd_seg)
segments(meanxy[,1]+sexy[,1],meanxy[,2]-lgy,meanxy[,1]+sexy[,1],meanxy[,2]+lgy,col=colb, lwd=lwd_seg)
segments(meanxy[,1]-lgx,meanxy[,2]-sexy[,2],meanxy[,1]+lgx,meanxy[,2]-sexy[,2],col=colb, lwd=lwd_seg)
segments(meanxy[,1]-lgx,meanxy[,2]+sexy[,2],meanxy[,1]+lgx,meanxy[,2]+sexy[,2],col=colb, lwd=lwd_seg)
points(meanxy[,1],meanxy[,2],pch=pchp,col=colp,bg=bgp,cex=cexp) # mean
} # end of meansexy

###########################################################################################################
###########################################################################################################
# ADDSEGMENT : function to add a segment defined by an equation and the range on X and/or Y axes
# inputs : 
#  a, b : equation defining the segment (y=a+b*x)
# limX, limY : limit of the segment on both axes; limY is computed by default accoring to a,b and limX
# col_seg, lty_seg=1, lwd_seg : color, line type and width of the segment
addsegment<-function(a=0,b=1,limX,limY=sort(c(a+limX[1]*b,a+limX[2]*b)),col_seg="black",lty_seg=1,lwd_seg=1.5) {

x0<-limX[1] ; y0<-a+limX[1]*b
x1<-limX[2] ; y1<-a+limX[2]*b

if (y0<limY[1] ) {x0<-(limY[1]-a)/b  ; y0<-limY[1] }
if (y0>limY[2]) {x0<-(limY[2]-a)/b  ; y0<-limY[2] }

if (y1<limY[1] ) {x1<-(limY[1]-a)/b  ; y1<-limY[1] }
if (y1>limY[2]) {x1<-(limY[2]-a)/b  ; y1<-limY[2] }

segments(x0,y0,x1,y1,lty=lty_seg,lwd=lwd_seg,col=col_seg)

}  # end of function addsegment

###########################################################################################################
###########################################################################################################
# EQUA_LM : function to write the equation of a linear regression 
# inputs : 
# - reslm: output from the linear regression done using the "lm" function
# - nmy= name of Y variable (e.g. Y axis title)
# - nmx= name of X varaible (e.g. X axis title)
# - prec= number of decimals of coefficients a and b in the character string
# output :  an expression with equation of the regression "Y=a+b*X, R²=r²"
#             which can be added to a plot using "title" or "text" function

equa_lm<-function(reslm,nmy="Y",nmx="X",prec=1) {

a<-round_text(reslm$coefficients[1],prec)
b<-round_text(abs(reslm$coefficients[2]),prec)
signb<-"+" ; if(reslm$coefficients[2]<0) signb<-"-"
r2<-round_text(summary(reslm)$r.squared,2)
pval<-pvalue_star(summary(reslm)$coefficients[2,4],ns="")

res<-substitute(y*"="*a*signb * b *" "* scriptstyle("X")*" "* x*", R²=" * r2 * pval,
		list(y=nmy,x=nmx,a=a,b=b,signb=signb,r2=r2,pval=pval) )
return(res)
} # end of equa_lm



#	-------------------------------------------------------------------------------------
	# Fonction qui enlève les niveaux des lites

flattenlist <- function(x){  
  morelists <- sapply(x, function(xprime) class(xprime)[1]=="list")
  out <- c(x[!morelists], unlist(x[morelists], recursive=FALSE))
  if(sum(morelists)){ 
    Recall(out)
  }else{
    return(out)
  }
}                     
 
 
