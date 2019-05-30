tn = 20  # years
t0 = 0    # initial time
h  = 0.001 

t = seq(from=t0, to=tn, by=h)

nt = length(t)

Tr = rep(NA, nt,)
Is = rep(NA, nt,)
Ir = rep(NA, nt,)
N = rep(NA, nt,)
FrIs = rep(NA, nt,)
FrIr = rep(NA, nt,)
Fr = rep(NA, nt,)

# Initial conditions
Tr[1] = 10
Is[1] = 2
Ir[1] = 1
N[1] = Tr[1]+Is[1]+Ir[1]
FrIs[1] = Is[1]/N[1]
FrIr[1] = Ir[1]/N[1]
Fr[1] = (Is[1]+Ir[1])/N[1]

# Value of Parameters
lambda = 5
delta = 0.5
br = 0.25
bs = 0.24
a = 1
mu = 10^(-5)
epsilon = 0.1

Tr[1] = lambda/delta

for (j in 1:(nt-1) ){
  Tr[j+1] = Tr[j] + h*( lambda - delta*Tr[j] - ( (1-epsilon)*bs*Is[j] + br*Ir[j])*Tr[j]  )
  
  Is[j+1] = Is[j] + h*( (1-epsilon)*(1-mu)*bs*Is[j]*Tr[j] - a*Is[j] )
  Ir[j+1] = Ir[j] + h*( br*Tr[j]*Ir[j] - a*Ir[j] + (1-epsilon)*mu*bs*Is[j]*Tr[j]  )

  N[j+1] = Tr[j+1]+Is[j+1]+Ir[j+1]
  
  FrIs[j+1] = Is[j+1]/N[j+1]
  FrIr[j+1] = Ir[j+1]/N[j+1]
  
  Fr[j+1] = (Is[j+1]+Ir[j+1])/N[j+1]
}

#dev.off()

#plot(t[which(t == 5):length(t)], Tr[which(t == 5):length(t)], type="l", col="blue", ylim=c(0,max(Tr, Is, Ir)) )
plot(t, Tr, type="l", col="blue", ylim=c(0,max(Tr, Is, Ir)) )
lines(t, Is, type="l", col="red")
lines(t, Ir, type="l", col="green")

plot(t, Fr, type="l", col="black", ylim=c(0,max(Fr, FrIs, FrIr)) )
lines(t, FrIs, type="l", col="red")
lines(t, FrIr, type="l", col="blue")
