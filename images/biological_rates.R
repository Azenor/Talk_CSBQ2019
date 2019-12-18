##### Figure Metabolic Theory #####
###################################
### Temperature dependence of biological rates ###

## BA unimodal model
rBA=function(temp,par_br,uni=1) {
  #input: temp=temperature range
  #       par_br=parameter vector
  #parameters:
  r0=par_br[[1]]  #scaling coefficient
  tpk=par_br[[2]] #optimal temperature
  E=par_br[[3]]   #activation energy
  Ed=par_br[[4]]  #deactivaiton energy
  k=8.617*10^-5   #Boltzmann constant
  # uni : specify if unimodal or exponential BA model

  if (uni==0) {
    l=1
  } else {
    l=1/(1+exp(-1/(k*temp)*(Ed-(Ed/tpk+k*log(E/(Ed-E)))*temp))) #decline phase
  }
  return(BR=r0*exp(-E/(k*temp))*l)
}


### Interaction strength per capita, per population, net ###
########## Plots according to temperature ##################
############################################################

## Temperature

temp_seq=seq(285,315,length.out=50)

## Body mass

mH=1.34*10^-2 # herbivore body-mass


### Parameters temperature dependent

par1=c(aPH0=2*10^11,topt=298,E=0.6,E2=1.15)
b1=rBA(temp=temp_seq,par_br=par1, uni = 0)

par2=c(aPH0=2*10^11,topt=298,E=0.61,E2=1.15)
b2=rBA(temp=temp_seq,par_br=par2, uni = 0)

par3=c(aPH0=2*10^11,topt=298,E=0.63,E2=1.15)
b3=rBA(temp=temp_seq,par_br=par3, uni = 0)

par4=c(aPH0=2*10^11,topt=298,E=0.65,E2=1.15)
b4=rBA(temp=temp_seq,par_br=par4, uni = 0)

# function to draw curly braces
# x, y position where to put the braces
# range is the length of the brace
# position: 1 vertical, 2 horizontal
# direction: 1 left/down, 2 right/up
# depth controls width of the shape

CurlyBraces <- function(x0, x1, y0, y1, pos = 1, direction = 1, depth = 1) {

    a=c(1,2,3,48,50)    # set flexion point for spline
    b=c(0,.2,.28,.7,.8) # set depth for spline flexion point

    curve = spline(a, b, n = 50, method = "natural")$y * depth

    curve = c(curve,rev(curve))

    if (pos == 1){
        a_sequence = seq(x0,x1,length=100)
        b_sequence = seq(y0,y1,length=100)
    }
    if (pos == 2){
        b_sequence = seq(x0,x1,length=100)
        a_sequence = seq(y0,y1,length=100)
    }

    # direction
    if(direction==1)
        a_sequence = a_sequence+curve
    if(direction==2)
        a_sequence = a_sequence-curve

    # pos
    if(pos==1)
        lines(a_sequence,b_sequence, lwd=2,   xpd=NA) # vertical
    if(pos==2)
        lines(b_sequence,a_sequence, lwd=2, xpd=NA) # horizontal

}

png("images/BAfunction_mismatch.png")
par(xaxs='i',yaxs='i')

par(mgp = c(1.5, 0.3, 0), tck = -.015, family = 'sans', oma=c(0,0,0,6))

plot(temp_seq,b1,lwd=2,axes=F,ann=F,type='l', ylim = c(min(b1,b2,b3,b4), max(b1,b2,b3,b4)), col = "orange")
lines(temp_seq, b2, lwd=2, col = 2)
lines(temp_seq, b3, lwd=2, col = 3)
lines(temp_seq, b4, lwd=2, col = 4)
abline(v=min(temp_seq), lwd = 4)
abline(h=min(b1,b2,b3,b4), lwd = 4)

mtext(text=expression(Biological~rate~b[i]),side=2,line=1.5,font=1,cex=2.5)
mtext(text='Temperature',side=1,line=1.5,font=1,cex=2.5)

legend("topleft", c(expression(D[1], D[2], A, epsilon)), lty = c(1,1,1,1), lwd = 3, col = c("orange", 2, 3, 4), bty = 'n', cex = 3)
CurlyBraces(x0=316,  x1=316,  y0=51, y1=34, pos = 1, direction = 1, depth=1.5)
CurlyBraces(x0=316,  x1=316,  y0=36, y1=15, pos = 1, direction = 1, depth=1.5)
text(320, 42.5, expression(Delta~D), cex = 2, xpd = NA)
text(321, 25.5, expression(A-D[2]), cex = 2, xpd = NA)

dev.off()
