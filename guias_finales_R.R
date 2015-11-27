v <- numeric(3);v
v[3] <- 17; v
x <- c(2, 4, 3.1, 8, 6)
x
is.integer(x) 
is.double(x)
length(x)
x <- edit(x)
y = 1:4; y
y[2] <- 5
u <- 1:12
u
u1=u[2 * 1:5]
u1
assign("z", c(x, 0, x))
z
s1 <- seq(2, 10); s1 
s2 = seq(from=-1, to=5); s2
s3<-seq(to=2, from=-2); s3
s4=seq(from=-3, to=3, by=0.2); s4
s5 <- rep(s3, times=3); s5
1/x
v=2*x+z+1
v
e<- c(1, 2, 3, 4); e2<-c(4, 5, 6, 7); crossprod(e, e2) 
xt = t(x)
xt
u = exp(y);u
options(digits=10); u
resum <- c(length(y),sum(y), prod(y), min(y), max(y)); resum
yo <- sort(y); yo
deptos <- c("Santa Ana", "Sonsonate", "San Salvador"); deptos
deptos[4]="Ahuachapán"; deptos
codDeptos <- c(11, 12, 13, 14)
Oriente <- codDeptos [c("La Unión", "San Miguel")];Oriente
etiqs<-paste(c("X", "Y"), 1:10, sep=""); etiqs
M <- matrix(numeric(), nrow = 3, ncol=4)
M[2,3] <- 6
M
A <- matrix(c(2, 4, 6, 8, 10, 12), nrow=2, ncol=3)
A
mode(A)
dim(A)
attributes(A)
is.matrix(A)
is.array(A)
B <- matrix(1:12, nrow=3, ncol=4)
B
x1 <- seq(0, 10, 2); x1
x2 <- seq(1, 11, 2); x2
x3 <- runif(6); x3
Xcol <- cbind(x1, x2, x3); Xcol
Xfil <- rbind(x1, x2, x3); Xfil
X <- Xfil[1:3, c(2, 3)]; X
v<-c(1, 2); v %*%A
P <- A %*% B; P
2*A
length(A)
T=sqrt(B); T
t(A)
C <- matrix(c(2, 1, 10, 12), nrow=2, ncol=2); C
det(C)
InvC <- solve(C)
eigen(C)
c(length(A), sum(A), prod(A), min(A), max(A))
nombres <- matrix(c("Carlos", "José", "Caren", "René", "Mar??a", "Mario"),
                  nrow=3, ncol=2); nombres
X <- array(c(1, 3, 5, 7, 9, 11), dim=c(2, 3)); X
Z <- array(1, c(3, 3)); Z
W <- 2*Z+1
W
TX <- t(X)
TX
a <- c(2, 4, 6)
a
b <- 1:3
b
M <- a %o% b
M 
Arreglo3 <- array(c(1:8, 11:18, 111:118), dim = c(2, 4, 3))
Arreglo3
#practica 3
sexo <- c("M", "F", "F", "M", "F", "F", "M")
sexo
edad <- c(19, 20, 19, 22, 20, 21, 19)
edad
FactorSexo = factor(sexo)
FactorSexo
mediaEdad <- tapply(edad, FactorSexo, mean) 
mediaEdad
is.vector(mediaEdad)
is.matrix(mediaEdad)
is.list(mediaEdad)
is.table(mediaEdad)
is.array(mediaEdad)
factor()
lista1<-list(padre="Pedro", madre="Mar?a", no.hijos=3, edad.hijos=c(4,7,9))
lista1
is.matrix(lista1)
is.vector(lista1$edad.hijos)
lista1["madre"]
lista1[[4]][2]
lista1["padre"]
lista1$padre
lista1$edad.hijos[2]
lista1[[4]][2]
lista1[["pedro"]]
x <- "nombre"; lista1[x]
subLista <- lista1[4]; subLista
lista1[5] <- list(sexo.hijos=c("F", "M", "F")); lista1
lista1 <- edit(lista1)
S <- matrix(c(3, -sqrt(2), -sqrt(2), 2), nrow=2, ncol=2);S
autovS <- eigen(S); autovS
evals <- eigen(S)$values; evals
Notas <- matrix(c(2, 5, 7, 6, 8, 2, 4, 9, 10), ncol=3,
                dimnames=list(c("Matem?tica","?lgebra","Geometr?a"),
                              c("Juan","Jos?","Ren?"))); Notas
ncol=(3)
log <- sample(c(TRUE, FALSE), size = 20, replace = TRUE)
log

comp <- rnorm(20) + runif(20) * (1i)
comp
num <- rnorm(20, mean=0, sd=1)
num
df1 <- data.frame(log, comp, num)
df1
nombres <- c("logico", "complejo", "numerico")
names(df1) <- nombres; df1
row.names(df1) <- letters[1:20]
df1
edad <- c(18, 21, 45, 54); edad
datos <- matrix(c(150, 160, 180, 205, 65, 68, 65, 69), ncol=2, dimnames=list(c(),
                                                                             c("Estatura","Peso"))); datos
sexo <- c("F", "M", "M", "M"); sexo
hoja1 <- data.frame(Edad=edad, datos, Sexo=sexo)
hoja1
search()
attach(hoja1)
search()
Edad
hoja1$Peso <- Peso+1
hoja1
detach(hoja1)
edad
#practica 4
Entrada1 <- read.table("datos01.txt", header=TRUE)
Entrada1
Edat1 <- scan("datos01.txt", list(X1=0, X2=0), skip = 1, flush = TRUE, quiet = TRUE)
Edat1
pp <- scan("datos02.txt", skip = 1, quiet= TRUE)
pp

library(foreign)
baseproductos <-read.table("productos.csv",header=TRUE,sep = ",")
baseproductos

library(Hmisc)
Baseimportante<-spss.get("Mundo.sav",use.value.labels =TRUE)
Baseimportante
#Practica 5
library(splines)
library( RcmdrMisc)
library(car)
library(sandwich)
Entrada <- readXL("F:/GCarenpr?cticas/problema5.xlsx", rownames=FALSE, 
                  header=TRUE, na="", sheet="Hoja1", stringsAsFactors=TRUE)
library(relimp, pos=15)
Entrada

Datos <- readXL("F:/GCarenpr?cticas/problema7.xlsx", rownames=FALSE, 
                header=TRUE, na="", sheet="Hoja1", stringsAsFactors=TRUE)
Datos
#practica 6
#"CC"=Coca_Cola
#"PC"=Pepsi_Cola
#"SC"=Salva_Cola
Tipo<-c("CC","PC","SC");Tipo
Consumo<-sample(Tipo,20,replace=TRUE);Consumo
data.entry(Consumo)
write(Consumo, "Consumo.txt")

frec <- table(Consumo); frec
prop <- table(Consumo)/length(Consumo); prop

summary(Consumo)

barplot(frec, main="Gráfico de barras", xlab=" Consumo", col=c("yellow", "white", "red"),
        sub="Agosto-2012")

barplot(prop, main="Gráfico de barras", xlab=" Consumo\n", col=c("yellow", "white",
                                                                 "red"), sub="Agosto-2012")

pie(frec, main="Gráfico de pastel", xlab="Tipo de Consumo", col=c("yellow", "white",
                                                                  "cyan"), sub="Agosto-2012")


names(frec) = c("Coca Cola", "Pepsi", "Salva Cola")
pie(frec, main="Gráfico de pastel", xlab=" Consumo", radius=1, col=c("red", "gray",
                                                                     "cyan"), sub="Agosto-2012")


n <- length(frec)
hoja <- data.frame(frec); hoja
etiq <- c(paste(hoja$Var1, "-", hoja$Freq)); etiq
pie(frec, main="Gráfico de pastel", labels=etiq, col=rainbow(n), border=TRUE)
#practica 7
Hijos<-c(2,1,2,1,4,2,3,0,2,3,3,2,1,0,2,4,1,2,1,3,4,1,2,3,1,5,2,3,1,2)
data.entry(Hijos)
Hijos
length(Hijos)

write(Hijos, "Hijos.txt")
ls()
rm(list=ls(all=TRUE)); ls()

X <- scan("Hijos.txt", what = integer(0), na.strings = "NA", flush=FALSE)
ls()

stripchart(X, method="stack", vertical=FALSE, col="blue", pch=1, main="Gráfico de\n
           puntos", xlab="Número de hijos")

fab <- table(X); fab
fre <- fab/length(X); fre
Fac <- cumsum(fab); Fac
Far <- Fac/length(X); Far

options(digits=2)
tabla <- data.frame(fab=fab, fre=fre, Fac=Fac, Far=Far)
names(tabla) <- c("X", "fab", "free.X", "fre", "Fac", "Far")
tabla
tfre <- data.frame(X=tabla$X, fab=tabla$fab, fre=tabla$fre, Fac=tabla$Fac, Far=tabla$Far)
tfre

media <- mean(X, na.rm = FALSE); media

for(i in 1:length(X)) if (fab[i] == max(fab)) break()
moda <- names(fab[i]); moda # R no tiene incorporada una función para la moda

mediana <- median(X); mediana

range(X)

cuasivar <- var(X); cuasivar
s <- sd(X); s

quantile(X,c(0.25, 0.5, 0.75))

quantile(X, 0.6)

resumen <- summary(X); resumen

fivenum(X)


barplot(tfre[[2]], main="Gráfico de barras", xlab="X = Número Hijos\n", ylab="frecuencia",
        col=c("yellow", "blue", "white", "orange", "cyan", "red"), sub="Agosto-2012")


pie(tfre[[2]], main="Gráfico de pastel", xlab="Número Hijos \n", col=c("yellow", "blue",
                                                                       "white", "orange", "cyan", "red"), sub="Agosto-2012")


names(fab) = c("Cero", "Uno", "Dos", "Tres", "Cuatro", "Cinco")
pie(fab, main="Gráfico de pastel", xlab="X = Número Hijos\n", col=c("yellow", "blue",
                                                                    "white", "orange", "cyan", "red"), sub="Agosto-2012")

boxplot(X, main="Gráfico de caja", ylab="Número de hijos\n")

boxplot(X, main="Gráfico de caja", xlab=" Número de hijos\n", plot=TRUE, border="red",col="yellow", horizontal=TRUE)
#practica 8
Notas<-c(4.47,4.47,3.48,5.0,3.42,3.78,3.1,3.57,4.2,4.5,3.6,3.75,4.5,2.85,3.7,4.2,3.2,4.05,4.9,5.1,5.3,4.16,4.56,3.54,3.5,5.2,4.71,3.7,4.78,4.14,4.14,4.8,4.1,3.83,3.6,2.98,4.32,5.1,4.3,3.9,3.96,3.54,4.8,4.3,3.39,4.47,3.19,3.75,3.1,4.7,3.69,3.3,2.85,5.25,4.68,4.04,4.44,5.43,3.04,2.95);Notas
data.entry(Notas)
Notas
length(Notas)

write(Notas, "Notas.txt")

ls()
rm(list=ls(all=TRUE))
ls()

X <- scan("Notas.txt", what = double(0), na.strings = "NA", flush=FALSE)
ls()

n <- length(X); n
k <- 1+3.322*logb(60, 10); k
k <- round(k); k
rango <- max(X)-min(X); rango
a=rango/k; a
a <- round(a, 3); a

# Calcula el ancho o amplitud a de cada intervalo a=rango/k
rango <- max(X)-min(X); rango
a=rango/k; a
a <- round(a, 3); a


# Define los l??mites y puntos medios de cada uno de los k intervalos
limites <- seq(from=min(X)-0.01/2, to=max(X)+0.01/2, by=a); limites
options(digits=4)
ci <- cbind(1:k); ci
for(i in 2:length(limites)) ci[i-1, 1] <- (limites[i] + limites[i-1])/2
ci

# Encuentra las frecuencias absolutas fi para cada intervalo

options(digits=2)
fi <- cbind(table(cut(X, breaks = limites, labels=NULL, include.lowest=FALSE,
                      right=FALSE, dig.lab=4))); fi


# Encuentra las frecuencias relativas o proporciones fri

options(digits=4)
fri <- fi/n; fri

# Encuentra las frecuencias acumuladas ascendentes Fi
options(digits=2)
Fi <- cumsum(fi); Fi

# Encuentra las frecuencias relativas acumuladas Fri
options(digits=4)
Fri <- Fi/n; Fri

# Completa la tabla de frecuencias.
tablaFrec <- data.frame(ci=ci, fi=fi, fri=fri, Fi=Fi, Fri=Fri); tablaFrec

h <- hist(X, breaks=c(limites[1]-a, limites, limites[k+1]+a), freq = TRUE, probability = FALSE,
          include.lowest = FALSE,right = TRUE, main = "Histograma de frecuencias",
          col="lightyellow", lty=1, border="purple", xlab=" Notas de aspirantes", ylab="Frecuencia (fi)",
          axes=TRUE, labels=FALSE)
text(h$mids, h$density, h$counts, adj=c(0.5, -0.5), col="red")
rug(jitter(X)) # adiciona marcas de los datos

# h es un objeto del tipo lista que contiene atributos del histograma
is.list(h); h

h <- hist(X, breaks=c(limites[1]-a, limites, limites[k+1]+a), freq = FALSE,
          probability = TRUE, include.lowest = FALSE, right = TRUE,
          main="Aproximación a una Normal\n", col="lightyellow",lty=1,border="purple",
          xlab="Notas de aspirantes\n", ylab="Frecuencia relativa (fri)",
          axes=TRUE, labels=FALSE)
text(h$mids, h$density, h$counts, adj=c(0.5, 0.2), col="red")
rug(jitter(X)) # adiciona marcas de los datos
curve(dnorm(x, mean=mean(X), sd=sd(X)), col = 2, lty = 2,lwd = 2, add = TRUE)


#Crea el pol??gono de frecuencias
h <- hist(X, breaks=c(limites[1]-a, limites, limites[k+1]+a), freq = TRUE,
          probability=FALSE, include.lowest=FALSE,right=TRUE,
          main = "Pol??gono de frecuencias",col="lightyellow", lty=1, border="purple", xlab="
          Notas de aspirantes", ylab="Frecuencia (fi)", axes=TRUE, labels=FALSE)
text(h$mids, h$density, h$counts, adj=c(0.5, -0.5), col="red")
rug(jitter(X)) # adiciona marcas de los datos
vCi <- c(h$mids[1]-a, h$mids, h$mids[k+1]+a); vCi
vfi <- c(0, h$counts, 0); vfi
lines(vCi, vfi, col="blue", type="l")

#Crea la Ojiva ascendente o pol??gono de frecuencias acumuladas ascendentes
Fia <- c(0, Fi); Fia
plot(limites, Fia, type = "p", pch=1, col = "blue", main="Ojiva ascendente",
     xlab="Notas de aspirantes", ylab="Frecuencia acumulada (Fi)")
text(limites, h$density, Fia, adj=c(0.5, -0.5), col="red")
lines(limites, Fia, col="black", type="l")

#Calcula los principales estad??sticos descriptivos de la variable
# Calcula la moda, ya que el R no proporciona una función para eso.
options(digits=9)
for(i in 1:k){ 
  if (fi[i] == max(fi))
    break()
  if(i > 1)
    moda <- limites[i]+((fi[i]-fi[i-1])/((fi[i]-fi[i-1])+(fi[i]-fi[i+1]) ))*a
  else moda <- limites[i]+(fi[i]/(fi[i]+(fi[i]-fi[i+1])))*a
}
moda
# Calcula los cuartiles: Q1, Q2, Q3
Q <- 1:3
for(v in 1:3) for(i in 1:k) if (Fi[i] > (v*25*n)/100)
{
  Q[v] <- limites[i]+(((25*v*n/100)-Fi[i-1])/fi[i])*a
  break
}
Q

#Calcula los principales estad??sticos.
estadisticos <- rbind(media=sum(tabEstad$cifi)/n, moda=moda, Q1=Q[1], Q2=Q[2], Q3=Q[3],
                      rango=max(X)-min(X), varianza=sum(tabEstad$ciMedia2fi)/n,
                      Desviacion=sqrt(sum(tabEstad$ciMedia2fi)/n),
                      CoeficienteVariacion=sqrt(sum(tabEstad$ciMedia2fi)/n)/(sum(tabEstad$cifi)/n),
                      CAfisher=(sum(tabEstad$ciMedia3fi)/n)/sqrt(sum(tabEstad$ciMedia2fi)/n)^3,
                      CoeficienteCurtosis=((sum(tabEstad$ciMedia4fi)/n)/sqrt(sum(tabEstad$ciMedia2fi)/n)^4)-3)
estadisticos

# Gráfico de cajas
boxplot(X, main="Gráfico de caja", xlab="Notas", notch=FALSE,
        data=parent.frame(), plot=TRUE, border="red", col="yellow",horizontal=TRUE)
#Observación: en la función boxplot(), s?? plot es FALSE se produce un resumen de los valores (los cinco números).


windows()
boxplot(X, main="Gráfico de caja", xlab="X = Notas", notch=TRUE,
        data=parent.frame(), plot=TRUE, border="red", col="yellow",horizontal=TRUE)

par(mfrow=c(1,2)) # Divide la ventana gráfica en dos partes (1 fila, 2 columnas)
mtext(side=3, line=0, cex=2, outer=T, "Titulo para Toda la Página")
hist(X); boxplot(X)

#Calcula los principales estad?sticos descriptivos de la variable
# Calcula la moda, ya que el R no proporciona una funci?n para eso.
options(digits=4)
for(i in 1:k) if (fi[i] == max(fi)) break()
if(i > 1) moda <- limites[i]+((fi[i]-fi[i-1])/((fi[i]-fi[i-1])+(fi[i]-fi[i+1]) ))*a
moda <- limites[i]+(fi[i]/(fi[i]+(fi[i]-fi[i+1])))*a
moda

#Varios gr?ficos en una misma ventana
par(mfrow=c(1,2)) # Divide la ventana gr?fica en dos partes (1 fila, 2 columnas)
mtext(side=3, line=0, cex=2, outer=T, "Titulo para Toda la P?gina")
hist(X); boxplot(X)
#practica 9
library(foreign)
HojaCat <- read.table("HojaCat.txt", header=TRUE)
HojaCat


#Conecta la hoja de datos a la segunda ruta o lista de búsqueda.
attach(HojaCat, pos=2) # pos especifica la posición donde buscar la conexión
search()

#Crea una tabla de contigencia o de doble entrada
tablaCont <- table(HojaCat); tablaCont
length(HojaCat)


# Encuentra la suma de cada fila de la tabla de contingencia
# Distribución marginal de X=Estado civil
suma.filas <- apply(tablaCont, 1, sum); suma.filas
# El 1 indica que son totales por fila

# Gráficos de barras para tabla de contingencia.
# Barras apiladas
barplot(t(tablaCont), main="Gráfico de barras (Estado, Ocupación)", xlab="Estado civil", ylab="Ocupación", legend.text=TRUE)

barplot(t(tablaCont), main="Gráfico de barras (Estado, Ocupación)", xlab="Estado civil", ylab="Ocupación", beside=TRUE, legend.text=TRUE)

# Guardar las todas las opciones iniciales y modificar número de decimales
op <- options()
options(digits=3) # sólo imprime 3 lugares decimales
options('digits')


# Proporciones basadas en el total de la muestra, la suma de filas y columnas suman 1

propTotal <- prop.table(tablaCont); propTotal

barplot(t(propTotal), main="Gráfico de barras (Estado, Ocupación)", xlab="Estado civil\n",ylab="Ocupación", beside=TRUE, legend.text=TRUE)


# Proporciones basadas en el total por fila, cada fila suma 1.

propFila <- prop.table(tablaCont, 1); propFila
# Total por fila se indica en 1
barplot(t(propFila), main="Gráfico de barras (Estado, Ocupación)", xlab="Estado civil\n",
        ylab="Ocupación", beside=TRUE, legend.text=TRUE)


propFila <- prop.table(tablaCont, 1); propFila
# Total por fila se indica en 1
barplot(t(propFila), main="Gráfico de barras (Estado, Ocupación)", xlab="Estado civil\n",
        ylab="Ocupación", beside=TRUE, legend.text=TRUE)

#Realizar la prueba o contraste Chi-cuadrado de independencia
prueba <- chisq.test(tablaCont); prueba

# Frecuencias absolutas esperadas para la prueba Chi-cuadrada
prueba$expected # fij = fi./No. column
#practica 10
A <- c(100,96,92,96,92); A
B <- c(76,80,75,84,82); B
C <- c(108,100,96,98,100); C
Baterias <- data.frame(procesoA=A, procesoB=B, procesoC=C); Baterias
# Para editar los datos puede utilizar la funci?n fix()
fix(Baterias)
write.table(Baterias, file="Baterias.txt", append=FALSE, quote=TRUE, sep=" ", na="NA",
            col.names=TRUE)
ls(); rm(list=ls(all=TRUE)); ls()
Baterias <- read.table("Baterias.txt", header=TRUE); Baterias
attach(Baterias, pos=2)
search()
stripchart(Baterias, main="Gr?fico de puntos para los tres procesos", method = "stack", vertical =
             FALSE, col="blue", pch=1, xlab="Duraci?n (semanas)", ylab="Proceso")
#Muestra un resumen estad?stico para los tres procesos.
summary(Baterias)

# Horizontal
boxplot(Baterias, width=NULL, varwidth=TRUE, names, add= FALSE, horizontal = TRUE,
        main="Gr?fico de caja por proceso", border=par("fg"), col=c("yellow", "cyan", "red"), xlab =
          "Duraci?n (semanas)", ylab="Proceso")

# Vertical
boxplot(Baterias, width=NULL, varwidth=TRUE, names, add= FALSE, horizontal = FALSE,
        main="Gr?fico de caja por proceso", border=par("fg"), col=c("yellow", "cyan", "red"), xlab =
          "Duraci?n (semanas)", ylab="Proceso")

#Presenta la matriz de covarianzas muestral.
options(digits=3) # s?lo imprime 3 lugares decimales
S <- var(Baterias); S

Baterias <- stack(Baterias); Baterias
names(Baterias) # Muestra los encabezados de los vectores


#Desconecta la hoja de datos de la segunda ruta o lista de b?squeda.
detach(Baterias, pos=2); search()

#An?lisis de una variable bidimensional

Fuma = c("Si","No","No","Si","No","Si","Si","Si","No","Si"); Fuma
Cantidad = c(1,2,2,3,3,1,2,1,3,2); Cantidad
Estudia <- data.frame(Fuma=Fuma, Cantidad=Cantidad); Estudia
fix(Estudia)
write.table(Estudia, file="Estudia.txt", append=FALSE, quote=TRUE, sep=" ", na="NA",
            col.names=TRUE)
write.table
ls()
rm(list=ls(all=TRUE))
ls()
Estudia <- read.table("Estudia.txt", header=TRUE)
Estudia
tablaCont <- table(Estudia)
tablaCont
options(digits=3) # s?lo imprime 3 lugares decimales
propTotal <- prop.table(tablaCont); propTotal
propFila <- prop.table(tablaCont, 1)
propFila
propCol <- prop.table(tablaCont, 2)
propCol
barplot(table(Estudia$Cantidad, Estudia$Fuma), beside = FALSE, horizontal=FALSE, main="Gr?fico de barras (Fuma, Cantidad de horas de estudio)", legend.text =T, xlab="Fuma", ylab="Cantidad de horas-estudio")

barplot(table(Estudia$Fuma, Estudia$Cantidad), beside = FALSE, horizontal=FALSE,main="Gr?fico
        de barras (Cantidad de horas de estudio,Fuma)", legend.text =T, xlab="Cantidad de horas-estudio",
        ylab="Fuma")
Fuma=factor(Estudia$Fuma); Fuma
barplot(table(Estudia$Cantidad, Estudia$Fuma), main="Gr?fico de barras (Fuma, Cantidad de horas
        de estudio)", xlab="Fuma", ylab="Cantidad de horas-estudio", beside=TRUE, legend.text=T)
barplot(table(Estudia$Cantidad, Estudia$Fuma), main="Gr?fico de barras (Fuma, Cantidad de horas
        de estudio)", xlab="Fuma", ylab="Cantidad de horas-estudio", beside=TRUE, legend.text=c("menor
                                                                                                que 5", "5-10", "mayor que 10"))
#practica 11
usuarios <- c(10, 15, 20, 20, 25, 30, 30); usuarios
tiempo = c(1.0, 1.2, 2.0, 2.1, 2.2, 2.0, 1.9); tiempo
Sistema <- data.frame(Usuarios=usuarios, Tiempo=tiempo);Sistema
fix(Sistema)

write.table(Sistema, file="Sistema.txt", append=FALSE, quote=TRUE, sep=" ", na="NA",
            col.names = TRUE)

ls(); rm(list=ls(all=TRUE)); ls()

#Recupera la hoja de datos.
Sistema <- read.table("Sistema.txt", header=TRUE); Sistema

#Conecta la hoja de datos a la segunda ruta o lista de búsqueda.
attach(Sistema, pos=2); search()

#Muestra un resumen de principales estad??sticos de las variables.
summary(Sistema)
cov(Sistema) # Matriz de covarianzas
cor(Sistema, use = "all.obs", method="pearson") # Matriz de correlaciones

#Elabora un gráfico de dispersión para analizar alguna relación entre las variables.
plot(Usuarios, Tiempo, xlim= c(5, 35), ylim= c(0.0, 2.5), type = "p", pch=1, col = "blue", main =
       "Gráfico de dispersión (Usuarios, Tiempo)", xlab="Número de usuarios", ylab="Tiempo de
     ejecución")

#Sin cerrar la ventana del gráfico anterior, ejecuta la siguiente instrucción
identify(Usuarios, Tiempo, n=1) # n=1 indica que solamente será un punto seleccionado

reg.Y.X <- lm(Tiempo ~ -1 + Usuarios, Sistema, na.action=NULL, method="qr", model=TRUE)
#-1 indica que no se toma en cuenta la constante en el modelo.
summary(reg.Y.X)

lines(Usuarios, 0.079437*Usuarios)

reg.anova <- anova(reg.Y.X); reg.anova
#practica 13
moneda <- c("C", "+"); moneda
n <- 10; n
espacio <- 1:54;espacio
# se define el tamaño de la muestra
n <- 6; n
muestra <- sample(espacio, n); muestra

# genera el espacio muestral del lanzamiendo de los dos dados
espacio = as.vector(outer(1:6, 1:6, paste)); espacio

# se define el tamaño de la muestra
n <- 4; n

# finalmente se selecciona la muestra
muestra <- sample(espacio, n, replace=TRUE); muestra

#genera el espacio muestral de las 52 cartas
naipe = paste(rep(c("A", 2:10, "J", "Q", "K"), 4), c("OROS","COPAS", "BASTOS",
                                                     "ESPADAS"));naipe

# se define el tamaño de la muestra
n <- 5; n

# se obtiene la muestra sin reemplazo (aunque no se especifique con replace=FALSE)
cartas <- sample(naipe, n) ; cartas

espacio <- function(num)
{
  numDiv7 <- numeric(0)
  ind <- 0
  for(i in 1:length(num))
    if ((num[i] %% 7)==0)
    {
      ind <- ind+1
      numDiv7[ind]=num[i]
    }
  return(numDiv7)
}
numeros <- 1:500

# generando el espacio muestral
s <- espacio(numeros); s

# seleccionando la muestra
muestra <- sample(s, 12, replace=TRUE); muestra
#practica 14
dbinom(4,8,0.5)
x <- 2; n=8; p=1/2

pbinom(x, size = n, prob = p, lower.tail=TRUE)

x <- 4; n=8; p=1/2

#primera forma
F <- 1 - pbinom(x, n, p, lower.tail=TRUE); F

#segunda forma
pbinom(4, size=8, prob=0.5, lower.tail=FALSE)

x <- 3; mu <- 6


ppois(x, lambda = mu, lower.tail=TRUE)

#primera forma

sum(dpois(c(6,7,8),lambda = 6))

# segunda forma
F8 <- ppois(8, lambda = 6, lower.tail=TRUE)
F5 <- ppois(5,lambda = 6, lower.tail=TRUE)
F8 - F5

n <- 30

#genera 30 valores de una distribución de Poisson con λ = 6
x <- rpois(n, lambda=mu)

#calcula las probabilidades para cada valor generado


y <- dpois(x, lambda=mu)

#genera el gráfico de distribución
plot(x, y, xlab="x", ylab="Función de probalidad", main="Distribución de Poisson: lambda = 6",
     type="h")

#une los puntos a las l??neas
points(x, y, pch=21)

x <- 0:2 
m = 11
n <- 4; k=2
# x define el número de globos con premio

# se construye la distribución de frecuencias del número de premios
Tabla <- data.frame(Probabilidad=dhyper(x, m, n, k))
rownames(Tabla) <- c(??Ningun premio?,”Solamenteuno?,”Dospremios”)
Tabla

# x define el número de intentos fallidos
x <- 0:5; p=0.1

Tabla <- data.frame(Probabilidad=dgeom(x, prob=p))

# nombrando las filas de la distribución de frecuencias
rownames(Tabla) <- c(“Ventaenelprimerintento”, “Venta en el segundointento?, “Ventaeneltercerintento?, “Ventaenelcuartointento”, “Ventaenelquintointento”, “Ventaen elsextointento?)

x=0; n=7; p=0.1
dbinom(x, n, p, log = FALSE)
#Practica 15. Distribucion de probabilidad continua
#2.Calculo de probabilidad 
#ejemplo 1
x <-55
a=0
b<-90
punif(x,min=a,max=b,lower.tail=TRUE)
F55=punif(55,min=a,max=b,lower.tail=TRUE)
F15=punif(15,min=a,max=b,lower.tail=TRUE)
F55-F15
(1-F55)*( F55-F15)
#ejemplo 2
p<-c(0.80)
media=5; 
d.t=1
qnorm(p, mean=media,sd=d.t,lower.tail=TRUE)
p<-c(0.80)
g.l<-10
qt(p,df=g.l,lower.tail=TRUE)
n<-16
x<-4.5 
mu=5
sigma=1 
d.t=sigma/sqrt(n)
pnorm(x,mean=mu,sd=d.t,lower.tail=FALSE)
#ejemplo 3
x<-5
teta=7
pexp(x,rate=1/teta,lower.tail=FALSE)
x<-3
teta=7
pexp(x,rate=1/teta,lower.tail=TRUE)
pexp(4,rate=1/teta,lower.tail=FALSE)
p<-0.9 
teta<-7
qexp(p,rate=1/teta,lower.tail=TRUE)
qexp(0.5,rate=1/teta,lower.tail=TRUE)
qexp(0.68,rate=1/teta,lower.tail=TRUE)
qexp(0.32,rate=1/teta,lower.tail=FALSE)
#3. gneracion de muestras aleatorias de las distribuciones
#ejemplo 1
min<--2
max<-4
x=runif(100,min,max)
x
hist(x,main="X ~Uniforme(min=-2,max=4",xlab="X",ylab="densidad de probabilidad",
     probability=TRUE,col="cyan")
curve(dunif(x,min,max),col="blue",add=TRUE)
#ejemplo 2
x.norm<-rnorm(n=200,mean=10,sd=2)
x.norm
hist(x.norm,breaks="Sturges",freq=TRUE,probability=FALSE,include.lowest=TRUE,right
     =TRUE,density=NULL,angle=45,col="steelblue1",border = NULL, main = "Histograma de
     datos observados",axes=TRUE,plot=TRUE,labels=FALSE)
plot(ecdf(x.norm),main="Función de distribución acumulada teórica")
#ejemplo 3
media<-4.5
desviacion<-0.75
x=rnorm(100,media,desviacion)
x
hist(x,main=expression(paste("X ~ N(",mu,"= 4.5,",sigma,"=0.75)")),xlab="X",ylab="densidad
     de probabilidad",probability=TRUE,col=gray(0.9))
hist
curve(dnorm(x,media,desviacion),col="red",lwd=2,add=TRUE)
curve
#Ejemplo 4
media<-2500
razon<-1/media
n=100 
x=rexp(n,razon)
x
hist(x,main="X ~ Exponencial( media = 2500 )", xlab="X", ylab="densidad de probabilidad", 
     probability=TRUE,col="cyan") 
curve(dexp(x,razon),col="blue",lwd=2,add=TRUE)
# 4 Funciones de distribucion y su inversa (los cuantiles)
#ejemplo 1
x <- 0.7 
p<-pnorm(x,mean=1,sd=1,lower.tail=TRUE)
p 
#ejemplo 2
z <- 0.7 
p1<-pnorm(z,mean=0,sd=1)
p1 
p2<-pnorm(z,mean=0,sd=1,lower.tail=FALSE)
p2
p3<-1-pnorm(z,mean=0,sd=1)
p3 
#ejemplo 3
p <- 0.75 
z<-qnorm(p,mean=0,sd=1,lower.tail=TRUE)
z
#ejemplo 4
x<-18.55
gl<-12 
p<-pchisq(x,gl,lower.tail=FALSE)
p 
#practica 16
#Simulacion del teorema del limite central
#ejemplo 1
tm=100
n<-10
p<-0.25
S=rbinom(tm,n,p)
Z=(S-n*p)/sqrt(n*p*(1-p))
Z
hist(Z,main="Histograma de Z ~ N(0, 1)", xlab="z = n?mero binomiales est?ndarizados", 
     ylab="f(z)", prob=TRUE, col="khaki") 
curve(dnorm(x,0,1),col = "deepskyblue",lty=2,lwd=2,add=TRUE)
#ejemplo 2
simulNorm <- function(mu,sigma, m=5, n=100) 
{ 
  vectMedias <<- numeric(0) 
  MediasEstand <<- numeric(0) 
  for (i in 1:m) 
  { 
    X = rnorm(n, mu, sigma) 
    vectMedias[i] <<- mean(X) 
    MediasEstand[i] <<- (vectMedias[i] - mu)/(sigma/sqrt(n)) 
  } 
} 
mu=5
sigma=5 
m <- 200
simulNorm(mu, sigma, m) 
hist(MediasEstand, main="Histograma de medias est?ndarizadas", xlab="Valores de m 
     medias normales est?ndarizadas", prob=TRUE, col="darkolivegreen3") 
curve(dnorm(x, 0, 1), col = "deeppink3", lty=2, lwd=2, add=TRUE)
qqnorm(MediasEstand, main="X ~ N(0, 1)") 
qqline(MediasEstand, lty=1, lwd=2, col="red")
simulExp <- function(mu, m=5, n=100) 
{ 
  razon <- 1/mu 
  vectMedias <<- numeric(0) 
  MediasEstand <<- numeric(0) 
  for (i in 1:m) 
  { 
    X = rexp(n, razon) 
    vectMedias[i] <<- mean(X) 
    MediasEstand[i] <<- (vectMedias[i] - mu)/(mu/sqrt(n)) 
  } 
}
par(mfrow=c(2,2)) 
mu=10 
m <- 100; n <- 1 
simulExp(mu, m, n) 
hist(MediasEstand, main="Medias Exp(10); n=1", xlab="m medias exp est?ndarizadas", 
     prob=TRUE, col="darkolivegreen3") 
xvals = seq(from=-3, to=3, by=0.01) 
points(xvals, dnorm(xvals, 0, 1), col = "red", type="l", lty=1, lwd=2) 
n <- 5 
simulExp(mu, m, n) 
hist(MediasEstand, main="Medias Exp(10); n=5", xlab="m medias exp est?ndarizadas", 
     prob=TRUE, col="darkolivegreen3") 
xvals = seq(from=-3, to=3, by=0.01) 
points(xvals, dnorm(xvals, 0, 1), col = "red", type="l", lty=1, lwd=2) 
n <- 15
simulExp(mu, m, n) 
hist(MediasEstand, main="Medias Exp(10); n=15", xlab="m medias exp est?ndarizadas", 
     prob=TRUE, col="darkolivegreen3") 
xvals = seq(from=-3, to=3, by=0.01) 
points(xvals, dnorm(xvals, 0, 1), col = "red", type="l", lty=1, lwd=2) 
n <- 50
simulExp(mu, m, n) 
hist(MediasEstand, main="Medias Exp(10); n=50", xlab="m medias exp est?ndarizadas", 
     prob=TRUE, col="darkolivegreen3") 
xvals = seq(from=-3, to=3, by=0.01) 
points(xvals, dnorm(xvals, 0, 1), col = "red", type="l", lty=1, lwd=2)
#practica 17
simulIntProp <- function(m=5, n=1, p, nivel.conf=0.95)
{
  X <- rbinom(m, n, p)
  # Matriz con 1000 valores aleatorios binomial(n,p), 50 muestras cada una de tama?o 20
  pe <<- X/n
  # Calcula la proporci?n estimada en cada una de las muestras.
  SE <<- sqrt(pe*(1-pe)/n)
  # Calcula la desviaci?n est?ndar estimada en cada una de las muestras.
  alfa <- 1-nivel.conf
  z <<- qnorm(1-alfa/2)
  Intervalo <<- cbind(pe - z*SE, pe + z*SE)
  # genera los extremos del intervalo de confianza
  nInter <<- 0
  # un contador para conocer en cu?ntos intervalos se encuentra la verdadera proporci?n.
  for(i in 1:m)
    if ((p >= Intervalo[i, 1]) && (p <= Intervalo[i, 2]))
      nInter <<- nInter + 1
  # funci?n que cuenta cu?ntos intervalos contienen el verdadero valor del par?metro.
  return(nInter)
}
n=20; m= 50; p=0.5; nivel.conf=0.95
simulIntProp(m, n, p, nivel.conf)
Intervalo # para visualizar cada uno de los intervalos generados
nInter
# para visualizar en cu?ntos de estos intervalos se encuentra la verdadera proporci?n.
#Gr?fico que muestra los intervalos de confianza de 95% que contienen y no contienen el verdadero valor del par?metro p.
matplot(rbind(pe - z*SE, pe + z*SE), rbind(1:m, 1:m), type="l", lty=1)
abline(v=p)
#practica 18
intervaloProp <- function(x, n, nivel.conf=0.95)
{
  pe <- x/n
  alfa <- 1-nivel.conf
  z <- qnorm(1-alfa/2)
  SE <- sqrt(pe*(1-pe)/n)
  print(rbind(pe, alfa, z, SE))
  LInf <- pe-z*SE
  LSup <- pe+z*SE
  print(" ")
  print(paste("Intervalo para p es: [", round(LInf, 2), ",", round(LSup, 2), "]"))
}
x=360; n=1200; nivel.conf=0.95
intervaloProp(x, n, nivel.conf)
#practica 21
IMC_Control <- c(23.6, 22.7, 21.2, 21.7, 20.7, 22.0, 21.8, 24.2, 20.1, 21.3, 
                 20.5, 21.1, 21.4, 22.2, 22.6, 20.4, 23.3, 24.8)
par(mfrow=c(1,2)) 



#Se genera el histograma de la variables de interés


hist(IMC_Control,main="A",xlab="IMC (kg/m2)",ylab="Frecuencia")



#Se genera el diagrama de caja de la variable de interés y se muestra en la misma ventana


boxplot(IMC_Control,main="B", lab="IMC (kg/m2)",ylim=c(20,25)) 




#Los commandos para contrastar normalidad son los siguientes


sw <- shapiro.test(IMC_Control)
sw 

ks <- ks.test(IMC_Control,"pnorm",mean=mean(IMC_Control),sd=sd(IMC_Control))
ks


#Luego se digitan los datos para pacientes y se ejecutan las mismas instrucciones

IMC_Pacientes <- c(25.6, 22.7, 25.9, 24.3, 25.2, 29.6, 21.3, 25.5, 27.4, 
                   22.3, 24.4, 23.7, 20.6, 22.8) 

par(mfrow=c(1,2))
hist(IMC_Pacientes,main="A",xlab="IMC (kg/m2)",ylab="Frecuencia")
boxplot(IMC_Pacientes,main="B", lab="IMC (kg/m2)",ylim=c(20,30))

sw <- shapiro.test(IMC_Pacientes)
sw 

ks <- ks.test(IMC_Pacientes,"pnorm",mean=mean(IMC_Pacientes),sd=sd(IMC_Pacientes))
ks
#practica 22
Prueba.prop <- function(x, n, po, H1="Distinto", alfa=0.05)
{
  op <- options();
  options(digits=2)
  pe=x/n 
  SE <- sqrt((po * (1-po))/n) 
  Zo <- (pe-po)/SE 
  
  if (H1 == "Menor" || H1 == "Mayor")
  {
    Z <- qnorm(alfa, mean=0, sd=1, lower.tail = FALSE, log.p = FALSE)
    
    valores <- rbind(Prop_Estimada=pe, Prop_Hipotetica=po, Z_critico=Z,Estadistico= Zo)
  }
  else
  {
    Z <- qnorm(alfa/2, mean=0, sd=1, lower.tail = FALSE, log.p = FALSE)
    
    valores <- rbind(Prop_Estimada=pe, Prop_Hipotetica =po, Z_critico_menor=-Z,
                     Z_critico_mayor =Z, Zo)
  } 
  if (H1 == "Menor")
  {
    if (Zo < -Z) decision <- paste("Como Estadistico <", round(-Z,3), 
                                   ", entonces rechazamos Ho")
    else decision <- paste("Como Estadistico>=", round(-Z,3), 
                           ", entonces aceptamos Ho")
  }
  if (H1 == "Mayor")
  {
    if (Zo > Z) decision <- paste("Como Estadistico >", round(Z,3), 
                                  ", entonces rechazamos Ho")
    else decision <- paste("Como Estadistico <=", round(Z,3), 
                           ", entonces aceptamos Ho")
  }
  if (H1 == "Distinto")
  {
    if (Zo < -Z) decision <- paste("Como Estadistico <", round(-Z,3), 
                                   ", entonces rechazamos Ho")
    if (Zo > Z) decision <- paste("Como Estadistico >", round(Z,3), 
                                  ", entonces rechazamos Ho")
    else decision <- paste("Como Estadistico pertenece a [", round(-Z,3), 
                           ",", round(Z,3), "], entonces aceptamos Ho") 
  } 
  print(valores)
  print(decision)
  options(op) 
} 

Prueba.prop(23, 100, 0.15, H1="Menor", alfa=0.05)
Prueba.prop(23, 100, 0.15, H1="Mayor", alfa=0.05)
Prueba.prop(23, 100, 0.15, H1="Distinto", alfa=0.05)

prop.test(x=23, n=100, p=0.15, alternative="less", conf.level=0.95)
prop.test(x=23, n=100, p=0.15, alternative="greater", conf.level=0.95)
prop.test(x=23, n=100, p=0.15, alternative="two.sided", conf.level=0.95) 
#PRUEBA DE HIPÓTESIS SOBRE UNA MEDIA, VARIANZA CONOCIDA.


X <- c(9.0, 3.41, 6.13, 1.99, 6.92, 3.12, 7.86, 2.01, 5.98, 
       4.15, 6.87, 1.97, 4.01, 3.56, 8.04, 3.24, 5.05, 7.37)


Prueba.param <- function(x, des, VEC, H1="Distinto", alfa=0.05)
{
  op <- options();
  options(digits=2)
  miu<- mean(VEC)
  RC<- 1.645
  L<- length(VEC)
  Zo<-(miu-x)/(((des^2)/L)^(0.5))
  Z<- 1.645
  
  if (H1 == "Menor")
  {
    if (Zo < -Z) decision <- paste("Como Estadistico <", round(-Z,3), 
                                   ", entonces rechazamos Ho")
    else decision <- paste("Como Estadistico>=", round(-Z,3), 
                           ", entonces aceptamos Ho")
  }
  if (H1 == "Mayor")
  {
    if (Zo > Z) decision <- paste("Como Estadistico >", round(Z,3), 
                                  ", entonces rechazamos Ho")
    else decision <- paste("Como Estadistico <=", round(Z,3), 
                           ", entonces aceptamos Ho")
  }
  if (H1 == "Distinto")
  {
    if (Zo < -Z) decision <- paste("Como Estadistico <", round(-Z,3), 
                                   ", entonces rechazamos Ho")
    if (Zo > Z) decision <- paste("Como Estadistico >", round(Z,3), 
                                  ", entonces rechazamos Ho")
    else decision <- paste("Como Estadistico pertenece a [", round(-Z,3), 
                           ",", round(Z,3), "], entonces aceptamos Ho") 
  } 
  
  print(decision)
  options(op) 
} 

Prueba.param(4, 2.45, X, H1="Mayor", alfa=0.05)

t.test(X,mu=4,alternative="greater") 
# GUIA 23
#PRUEBAS SOBRE DOS MUESTRAS INDEPENDIENTES
IMC_Control <- c(23.6, 22.7, 21.2, 21.7, 20.7, 22.0, 21.8, 24.2, 20.1,
                 21.3, 20.5, 21.1, 21.4, 22.2, 22.6,
                 20.4, 23.3, 24.8)
IMC_Pacientes <- c(25.6, 22.7, 25.9, 24.3, 25.2, 29.6, 21.3, 25.5, 27.4, 22.3,
                   24.4, 23.7, 20.6, 22.8)
t.test(IMC_Control, IMC_Pacientes, var.equal=TRUE, mu=0)
#Se concluye entonces que existe diferencia significativa en el IMC para ambos grupos
#de pacientes,#pues el p valor de la prueba resulta ser muy peque?o.

#PRUEBAS SOBRE DOS MUESTRAS PAREADAS
PAS.antes <- c(160,155,180,140,150,130,190,192,170,165)
PAS.despues <- c(139,135,175,120,145,140,170,180,149,146)
shapiro.test(PAS.antes)
shapiro.test(PAS.despues)
ks.test(PAS.antes,"pnorm",mean=mean(PAS.antes),sd=sd(PAS.antes))
ks.test(PAS.despues,"pnorm",mean=mean(PAS.despues),sd=sd(PAS.despues))
t.test(PAS.antes, PAS.despues, paired=TRUE, mu=0)
#El valor del estad?stico t es 4.0552, con gl = 9, P = 0.0029. Con estos resultados
#se rechaza 0 H y por lo tanto se concluye que la PAS antes y despu?s del 
#tratamiento es distinta, es decir, el tratamiento ha sido efectivo.

#PRUEBA DE HIP?TESIS ACERCA DE LA VARIANZA DE DOS POBLACIONES
Agente_A <- c(12, 11, 18, 16, 13)
Agente_B <- c(14, 18, 18, 17, 16)
var.test(Agente_A, Agente_B)
#Como el p valor es alto se concluye que las varianzas pueden considerarse iguales.


#EJERCICIO
#contraste de la igualdad de varianzas
Tabla_A <- c(2098,2082,2246,2340,2714,2777,2625,2388,2766,3112,3030,3375,3038,
             3017,3136,3204,3174,3220,3464,3870,3689,3783,3457,4151,4230,3707,
             415,4315,4790,4464,4499,4819,4739,4912,4494,5698,6349,6630,7585,8183)
Tabla_B<- c(1209,1115,1151,1208,1170,1198,1390,1480,1359,1337,1415,1530,1453,
            1324,1477,1501,1661,1562,1764,1796,1976,1802,2000,1923,2097,2110,
            2214,2069,2324,2309,2353,2091,2187,2399,2630,2722,2998,3392,3379,3627)
var.test(Tabla_A, Tabla_B)
#Como el p valor es bajo se concluye que las varianzas pueden considerarse distintas.

#contraste de igualdad de medias.
t.test(Tabla_A, Tabla_B, var.equal=TRUE, mu=0)
#Se concluye que existe diferencia significativa en 
#la densidad espectral para ambos grupos de pacientes,
#pues el p valor de la prueba resulta ser muy peque?o.
#GUIA 24
#AN?LISIS DE VARIANZA

#Ejemplo 1

notas <- c(20,18,18,23,22,17,15,13,21,15,20,13,12,16,17,21,15,13,12,15,18,20,
           18,17,10,24,16)
programas <- gl(n=3, k=9, labels=c("P1", "P2", "P3"))
datos <- data.frame(notas = notas, programas = programas);datos
mod1 <- aov(notas ~ programas, data = datos)
plot(mod1)

#Ejemplo 2
#(No esta la base de datos)
#-------------------página 178 a 182---------------------
#******************DISEÑOS POR BLOQUES*******************
bloques<-gl(n=4,k=1,length = 20)
bloques # Vector de bloque del experimento
tratamientos<-gl(n=5,k=4)
tratamientos # Vector de tratamientos de los novillos
peso<-c(0.9,1.4,1.4,2.3,3.6,3.2,4.5,4.1,0.5,0.9,0.5,0.9,3.6,3.6,3.2,3.6,1.8,1.8,0.9,1.4)
peso # Se han registrado los pesos de Novillos
datos2<-data.frame(bloques=bloques,tratamientos=tratamientos,peso=peso)
datos2 # Se ha registrado en una hoja de datos los resultados del experimento
mod2<-aov(peso ~ tratamientos + bloques,data = datos2)  # Se aplica el análisis de varianza
summary(mod2) # Muestra la tabla ANOVA del tratamiento

#--------------------página 183 a 188----------------------
#******************DISEÑOS BIFACTORIALES*******************
FactorA<-gl(n=4,k=8,length = 32)
FactorA # Definiendo el vector que contiene el Factor A
FactorB<-gl(n=4,k=2,length = 32)
FactorB # Definiendo el vector que contiene los tratamientos de los novillos
Porcentaje<-c(1.8,2.1,2.0,2.1,4.6,5.0,7.5,7.9,2.2,2.4,4.2,4.0,5.4,5.6,9.8,9.2,2.8,3.2,4.4,4.8,8.7,8.4,13.2,13.0,3.2,3.6,3.3,3.5,5.7,5.8,10.9,11.1)
Porcentaje # Definiendo los pesos de los novillos
datos3<-data.frame(FactorA=FactorA,FactorB=FactorB,Porcentaje=Porcentaje)
datos3 # Registrando en una hoja los datos del resultado del experimento
mod3<-aov(Porcentaje~FactorA*FactorB,data = datos3)
summary(mod3) # Se muestra la tabla ANOVA del experimento

