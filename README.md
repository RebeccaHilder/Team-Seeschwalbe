#### Team-Seeschwalbe

## Daten einlesen
dat <- read.csv("slick2-h1-h2.csv")
## Paket laden
library(moveHMM)
## Daten anzeigen
View(dat)
head(dat)

## Anzahl der Beobachtungen pro Tern
table(dat$ID)

## neue Spalte im Datensatz erstellen mit der Frequency
data2 <- dat$ID
x <- data.frame(data2, freq=ave(seq_along(data2), data2, FUN=length))
dat$freq <- x$freq
View(dat)

which(data_frame$Freq >400) #welche Terns haben mehr als 400 Beobachtungen
sum(data_frame$Freq >400) #wie viele Terns haben mehr als 400 Beobachtungen

## Datensatz nur mit Terns, die >400 bzw. >700 Beobachtungspunkte besitzen
dat_400 <- subset(dat, dat$freq > 400)
dat_700 <- subset(dat, dat$freq > 700)

## Daten aufbereiten
data <- prepData(dat_700, type="UTM")

## 2 State HMM fitten
?fitHMM
# mit turning anlge
mod <- fitHMM(data,nbStates=2, stepPar0 = c(0.1,0.25, 0.05, 0.05), anglePar0 = c(0,0,1,1),verbose=1, angleDist="vm",stationary=T)
plot(mod)
mod
