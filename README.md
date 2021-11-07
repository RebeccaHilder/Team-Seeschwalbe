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
