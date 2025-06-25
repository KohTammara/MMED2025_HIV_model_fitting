

cMax <-0.7; cRate <- 0.5; cStart <- 1998
tt <- seq(1976, 2025, by = 0.1)
control_effect <- pmin(1, cMax+(1-cMax)*exp(-(tt-cStart)*cRate))
plot(tt, control_effect, type = 'l')
abline(v=cStart)
cMax <-0.7; cRate <- 0.5; cStart <- 1990
control_effect <- pmin(1, cMax+(1-cMax)*exp(-(tt-cStart)*cRate))
abline(v=cStart, col= 'red')
lines(tt, control_effect, col = 'red')
cMax <-0.7; cRate <- 0.5; cStart <- 2002
control_effect <- pmin(1, cMax+(1-cMax)*exp(-(tt-cStart)*cRate))
lines(tt, control_effect, col = 'blue')
abline(v=cStart, col= 'blue')


cMax <-0.7; cRate <- 0.5; cStart <- 1990
tt <- seq(1976, 2025, by = 0.1)
control_effect <- pmin(1, cMax+(1-cMax)*exp(-(tt-cStart)*cRate))
plot(tt, control_effect, type = 'l')
abline(v=cStart)
cMax <-0.7; cRate <- 2; cStart <- 1990
control_effect <- pmin(1, cMax+(1-cMax)*exp(-(tt-cStart)*cRate))
abline(v=cStart, col= 'red')
lines(tt, control_effect, col = 'red')
cMax <-0.7; cRate <- 0.1; cStart <- 1990
control_effect <- pmin(1, cMax+(1-cMax)*exp(-(tt-cStart)*cRate))
lines(tt, control_effect, col = 'blue')
abline(v=cStart, col= 'blue')

cMax <-0.5; cRate <- 0.5; cStart <- 1990
tt <- seq(1976, 2025, by = 0.1)
control_effect <- pmin(1, cMax+(1-cMax)*exp(-(tt-cStart)*cRate))
plot(tt, control_effect, type = 'l', ylim = c(0,1))
abline(h=cMax)
cMax <-0.8; cRate <- 0.5; cStart <- 1990
control_effect <- pmin(1, cMax+(1-cMax)*exp(-(tt-cStart)*cRate))
abline(h=cMax, col= 'red')
lines(tt, control_effect, col = 'red')
cMax <-0.2; cRate <- 0.5; cStart <- 1990
control_effect <- pmin(1, cMax+(1-cMax)*exp(-(tt-cStart)*cRate))
lines(tt, control_effect, col = 'blue')
abline(h=cMax, col= 'blue')