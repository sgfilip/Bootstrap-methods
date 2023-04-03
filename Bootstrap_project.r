#-----------------------------------------------------------------------------
# 1. ESTYMACJA PARAMETROW MODELU METODAMI KLASYCZNYMI
#-----------------------------------------------------------------------------

# usuwanie niepotrzebnej kolumny 1
dane<- subset(dane, select=-c(obs))

# zaladowanie bibliotek
install.packages("dplyr")
install.packages("nlme")
install.packages("car")
install.packages("randtests")
install.packages("lmtest")
install.packages("strucchange")

library(dplyr)
library(nlme)
library(car)
library(randtests)
library(lmtest)
library(strucchange)

# PRZYJETY POZIOM ISTOTNOSCI - 5%

# MODEL LINIOWY Z WSZYSTKIMI ZMIENNYMI
model_1 <-lm(nocars~., data=dane) #wszystkie zmienne
summary(model_1)
# normalnosc rozkladu reszt
shapiro.test(model_1$residuals)
# WNIOSEK: p-value>5% => reszty maja rozklad normalny
# istotnosc zmiennych: nie wszystkie sa statystycznie istotne, nalezy
# dokonac wyboru postaci modelu metoda krokowa wstecz

# METODA KROKOWA WSTECZ
step(model_1,direction = "backward")
model_2 <- lm(formula = nocars ~ pop + income + primert, data = dane)
summary(model_2)

# normalnosc rozkladu reszt
shapiro.test(model_2$residuals)
# WNIOSEK: p-value>5% => reszty maja rozklad normalny
# istotnosc zmiennych: wszystkie zmienne sa statystycznie istotne

#-----------------------------------------------------------------------------
# 2. BADANIE AUTOKORELACJI MODELU 2
#-----------------------------------------------------------------------------

# autokorelacja - test Durbina-Watsona
dwtest(model_2)
# WNIOSEK: p-value>5% => brak autokorelacji

#-----------------------------------------------------------------------------
# 3. BADANIE HOMOSKEDASTYCZNOSCI MODELU 2
#-----------------------------------------------------------------------------

# homoskedastycznosc - test Breuscha-Pagana
bptest(model_2)
# WNIOSEK: p-value>5% => reszty sÄ… homoskedastyczne

#-----------------------------------------------------------------------------
# 4. PRZEDZIALY UFNOSCI - METODY BOOTSTRAP
#-----------------------------------------------------------------------------
# ramka danych tylko z Y i 3 zm. objasniajacymi
dane2<-subset(dane, select=-c(price, unemp))

# pusta ramka danych na parametry do I schemtu losowania
liczba_wierszy=1000
dane_boot<-data.frame(stala=numeric(liczba_wierszy), pop=numeric(liczba_wierszy),income=numeric(liczba_wierszy), 
                      primert=numeric(liczba_wierszy))

# I SCHEMAT LOSOWANIA - OBSERWACJE
for(i in 1:liczba_wierszy){
  proba_boot <- dane2[sample(row(dane2), size = nrow(dane2), replace = T), ]
  model_boot <- lm(nocars~pop+income+primert,data=proba_boot)
  dane_boot[i,] <- model_boot$coefficients
}

# pusta ramka danych do II schematu losowania
liczba_wierszy=60 # w tym wypadku musi byc tyle wierszy, ile mamy obserwacji, bo dodajemy reszty do wartosci tepretycznych y, ktorych mamy 60
dane_boot2<-data.frame(stala=numeric(liczba_wierszy), pop=numeric(liczba_wierszy),
                       income=numeric(liczba_wierszy), primert=numeric(liczba_wierszy))

# II SCHEMAT LOSOWANIA - RESZTY
proba2=model_2$residuals

for(i in 1:1000){
  nowe_y<-c()
  nowe_reszty=sample(proba2,liczba_wierszy,replace = T)
  nowe_y=model_2$fitted.values+nowe_reszty
  nowe_dane<-data.frame(nowe_y,dane$pop,dane$income,dane$primert)
  nowy_model<-lm(nowe_y~.,data=nowe_dane)
  dane_boot2[i,]<-nowy_model$coefficients # po 1000 wartosci
}

# ---------------------
#       NORMALNY
# ---------------------

# funkcja na przedziak‚ normalny, zeby nie liczyc tego 20 razy
przedzialy_norm <- function (proba)
{
  poczatek = mean(proba) - 1.96*sd(proba)
  koniec = mean(proba) + 1.96*sd(proba)
  return (sprintf("Przedzial ufnosci oparty na przyblizeniu normalnym: (%s, %s)", round(poczatek,3), round(koniec,3)))
}

# I schemat losowania
przedzialy_norm(dane_boot$stala)
przedzialy_norm(dane_boot$pop)
przedzialy_norm(dane_boot$income)
przedzialy_norm(dane_boot$primert)

# II schemat losowania
przedzialy_norm(dane_boot2$stala)
przedzialy_norm(dane_boot2$pop)
przedzialy_norm(dane_boot2$income)
przedzialy_norm(dane_boot2$primert)

# ---------------------
#        BASIC
# ---------------------

przedzialy_basic <- function(proba){
  srednia0 <- mean(proba)
  roznica <- proba-srednia0
  kwantyl_1 <- quantile(roznica, 0.975)
  kwantyl_2 <- quantile(roznica, 0.025)
  
  poczatek <- srednia0 - kwantyl_1
  koniec <- srednia0 - kwantyl_2
  return (sprintf("Przedzial ufnosci typu basic: (%s, %s)", round(poczatek,3), round(koniec,3)))
}

# I schemat losowania
przedzialy_basic(dane_boot$stala)
przedzialy_basic(dane_boot$pop)
przedzialy_basic(dane_boot$income)
przedzialy_basic(dane_boot$primert)

# II schemat losowania
przedzialy_norm(dane_boot2$stala)
przedzialy_norm(dane_boot2$pop)
przedzialy_norm(dane_boot2$income)
przedzialy_norm(dane_boot2$primert)

# ---------------------
#    PERCENTYLOWY
# ---------------------
przedzialy_perc <- function(proba){
  kwantyl_1 <- quantile(proba, 0.025)
  kwantyl_2 <- quantile(proba, 0.975)
  
  return (sprintf("Przedzial ufnosci percentylowy: (%s, %s)", round(kwantyl_1,3), round(kwantyl_2,3)))
}

# I schemat losowania
przedzialy_perc(dane_boot$stala)
przedzialy_perc(dane_boot$pop)
przedzialy_perc(dane_boot$income)
przedzialy_perc(dane_boot$primert)

# II schemat losowania
przedzialy_norm(dane_boot2$stala)
przedzialy_norm(dane_boot2$pop)
przedzialy_norm(dane_boot2$income)
przedzialy_norm(dane_boot2$primert)

# ---------------------
#    STUDENTYZOWANY
# ---------------------
srednia <- c()
t<-c()
parametry<-data.frame(stala=numeric(), pop=numeric(),income=numeric(), primert=numeric())
parametry2<-data.frame(stala=numeric(), pop=numeric(),income=numeric(), primert=numeric())

# przedzial studentyzowany dla I schematu losowania
# funkcja na przedzial studentyowany, parametr n to numer zmiennej: 1-STALA, 2-POP, 3-INCOME, 4-PRIMERT

przedzialy_student <- function(n){
  wart_est <- as.data.frame(model_2$coefficients)
  for (i in 1:1000){
    boot <- dane2[sample(row(dane2), size=60,replace = T),]
    model_boot <- lm(nocars~pop+income+primert,data=boot)
    parametry[i,] <- model_boot$coefficients
    for (j in 1:100){
      boot2 <- boot[sample(row(boot), size=60,replace = T),]
      model_boot2 <- lm(nocars~pop+income+primert,data=boot2)
      parametry2[j,] <- model_boot2$coefficients
    }
    srednia[i] <- mean(parametry2[,n])
    t[i] <- (parametry[i,n]-wart_est[n,1])/sd(parametry2[,n])
  }
  
  kwantyl_1 <- quantile(t, 0.975)
  kwantyl_2 <- quantile(t, 0.025)
  poczatek <- wart_est[n,1] - kwantyl_1*sd(parametry[,n])
  koniec <- wart_est[n,1] - kwantyl_2*sd(parametry[,n])
  return(sprintf("Przedzial ufnosci studentyzowany: (%s, %s)", round(poczatek,3), round(koniec,3)))
}

przedzialy_student(1) # stala
przedzialy_student(2) # pop
przedzialy_student(3) # income
przedzialy_student(4) # primert

# przedzial studentyzowany dla II schematu losowania
proba2=model_2$residuals

przedzialy_student_2 <- function(n){
  for(i in 1:1000){
    nowe_y<-c()
    nowe_reszty=sample(proba2,60,replace = T)
    nowe_y=model_2$fitted.values+nowe_reszty
    nowe_dane<-data.frame(nowe_y,dane$pop,dane$income,dane$primert)
    nowy_model<-lm(nowe_y~.,data=nowe_dane)
    parametry[i,]<-nowy_model$coefficients # po 1000 wartosci
    for(j in 1:100){
      nowe_y2<-c()
      nowe_reszty2=sample(nowy_model$residuals,60,replace = T)
      nowe_y2=nowy_model$fitted.values+nowe_reszty2
      nowe_dane2<-data.frame(nowe_y2,dane$pop,dane$income,dane$primert)
      nowy_model2<-lm(nowe_y2~.,data=nowe_dane2)
      parametry2[i,]<-nowy_model2$coefficients
    }
    srednia[i] <- mean(parametry2[,n])
    t[i] <- (parametry[i,n]-wart_est[n,1])/sd(parametry2[,n])
  }
  kwantyl_1 <- quantile(t, 0.975)
  kwantyl_2 <- quantile(t, 0.025)
  poczatek <- wart_est[n,1] - kwantyl_1*sd(parametry[,n])
  koniec <- wart_est[n,1] - kwantyl_2*sd(parametry[,n])
  return(sprintf("Przedzial ufnosci studentyzowany: (%s, %s)", round(poczatek,3), round(koniec,3)))
}

przedzialy_student_2(1) # stala
przedzialy_student_2(2) # pop
przedzialy_student_2(3) # income
przedzialy_student_2(4) # primert


#-----------------------------------------------------------------------------
# 5. BADANIE ISTOTNOSCI KAZDEGO Z PARAMETROW - METODY BOOTSTRAP
#-----------------------------------------------------------------------------
# wg metod klasycznych, wszystkie 3 zmienne sa istotne na poziomie istotnosci alfa = 0,05

# istotnosc badamy stosujac II schemat losowania
# zbieramy do zmiennych wartości statystyk z wyjsciowego modelu z 3 ziennymi
stat <-summary(model_2)
stat_F <- stat$fstatistic #  statystyka F
stat_T<-stat$coefficients[,3] # statytsyki T dla kazdej zmiennej

stat_F
stat_T

# n - numer kolumny ze zmienna, ktorej istotnosc sprawdzmy -> 2 - POP; 3 - INCOME; 4 - PRIMERT

istotnosc<-function(n){ 
  model <- lm(formula=nocars~.,data=dane2[,-n])
  liczba_wierszy=60
  dane_boot4_stat<-data.frame(stala=numeric(liczba_wierszy), pop=numeric(liczba_wierszy),
                              income=numeric(liczba_wierszy), primert=numeric(liczba_wierszy))
  reszty<-model$residuals # zapisujemy reszty
  
  for(i in 1:1000){
    nowe_y<-c()
    nowe_reszty=sample(reszty,liczba_wierszy,replace = T)
    nowe_y=model$fitted.values+nowe_reszty
    nowe_dane<-data.frame(nowe_y,dane$pop,dane$income,dane$primert)
    nowy_model<-lm(nowe_y~.,data=nowe_dane)
    tabelka<-summary(nowy_model)
    wyniki<-list(tabelka$coefficients[,3]) # zapisuje 1000 wartosci statystyk T
    dane_boot4_stat[i,]<-wyniki 
  }
  pvalue <- mean(abs(dane_boot4_stat[,n])<stat_T[n]) # wartosc pvalue pokaze sie pod wykresem -> jesli mniejsza niz 0,05, to zmienna jest istotna
  
  return (plot(density(dane_boot4_stat[,n]),sub=pvalue, main="Wykres gestosci statystyki t",
               xlab="P-value dla testu istotsnosci parametru:"))
}

istotnosc(2) # pop ---- p-value zawsze = 0 ---- zmienna ISTOTNA
istotnosc(3) # income ---- p-value waha sie w granicach 0,6-0,8 ---- zmienna NIEISTOTNA
istotnosc(4) # promert ---- p-value zawsze = 0 ---- zmienna ISTOTNA

# jedna ze zmiennych, ktore wg metody krokowej byla istotna, okazala sie nieiostotna przy zastosowaniu metod bootstrapowych

#-----------------------------------------------------------------
# ^. BADANIE ISTOTNOSCI WSZYSTKICH PARAMETROW
#----------------------------------------------------------------- 

model_3 <- lm(formula=nocars~.,data=subset(dane2,select=c(nocars))) # bierzemy tylko zmienna objasniana, szacujemy stala, y = a
stat_F # wartsc statystyki F dla modelu wyjsciowego

liczba_wierszy=60
reszty2<-model_3$residuals # zapisujemy reszty

dane_stat_F<-c() # tu bedziemy wczytywac wartosci statystyki F w kazdej losowanej probie

for(i in 1:1000){
  nowe_y<-c()
  nowe_reszty=sample(reszty2,liczba_wierszy,replace = T)
  nowe_y=model_3$fitted.values+nowe_reszty
  nowe_dane<-data.frame(nowe_y,dane$pop,dane$income,dane$primert)
  nowy_model<-lm(nowe_y~.,data=nowe_dane)
  tabelka<-summary(nowy_model)
  dane_stat_F[i]<-tabelka$fstatistic[1]
}

pvalue2 <- mean(dane_stat_F>stat_F[1])  # p-value zawsze = 0 ---- laczenie wszystkie zmienne sa ISTOTNE

plot(density(dane_stat_F),sub=pvalue2,main="Wykres gestosci statystyki F",
     xlab="P-value dla testu istotsnosci wszystkich zmiennych:") # wartosc pierwotna statystyki F nie zawiera sie w przedziale wartosci wykresu gestosci -> ISTOTNE WSZYSTKIE





