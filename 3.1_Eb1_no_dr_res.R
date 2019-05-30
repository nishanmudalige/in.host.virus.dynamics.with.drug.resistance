tn = 20  # years
t0 = 0    # initial time
h  = 0.01 

t = seq(from=t0, to=tn, by=h)

nt = length(t)

Tr = rep(NA, nt,)
I  = rep(NA, nt,)
N  = rep(NA, nt,)
Fr = rep(NA, nt,)

# Initial conditions
Tr[1] = 10
I[1] = 2
N[1] = Tr[1] + I[1]
Fr[1] = I[1]/N[1]

# Value of Parameters
lambda = 5
delta = 0.5
b = 0.25
a = 1
epsilon = 0.1


for (j in 1:(nt-1) ){
  Tr[j+1] = Tr[j] + h*( lambda - delta*Tr[j] - (1-epsilon)*b*I[j]*Tr[j] )
  I[j+1]  = I[j] + h*( (1-epsilon)*b*I[j]*Tr[j] - a*I[j] )
  N[j+1]  = Tr[j+1] + I[j+1]
  Fr[j+1] = I[j+1]/N[j+1]
}
# plot(t, N, type="l", col="black", ylim=c(0,max(N,Tr, I)))
# lines(t, Tr, type="l", col="blue")

par(mfrow=c(2,1))
plot(t, Tr, type="l", col="blue", ylim=c(0,max(Tr, I)) )
lines(t, I, type="l", col="red")
plot(t, Fr, type="l")

# rm(list = ls()

