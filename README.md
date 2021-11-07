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
