mu5<-function(thet,x,tim){

theta<-exp(thet)
time<-tim*60

lsexamp0 <- function(t, y, p) {
#lsexamp is a function describing the differential equations
out<-odes5.cpp(t=t,y=y,p=p)
            list(out)
        }
        pars <- c(theta[1], theta[2], theta[3], theta[4], 
            theta[4], theta[5], 8.3145, x[4] + 273.15, 25 + 273.15)

        yini <- c(0, 0, x[2], x[3], 0, 0,x[1])/x[5] 

        out <- x[5]*ode(func = lsexamp0, y = yini, rtol =1e-15, atol = 1e-15, 
            parms = pars, method = "vode", times = c(0, time))
        
        out2 <- matrix(out[-1, c(3,5,7)],ncol=3)

log(out2[,c(2,1,3)])}





