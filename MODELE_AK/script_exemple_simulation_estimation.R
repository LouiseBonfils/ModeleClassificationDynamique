
setwd("C:/Users/bonfils/Documents/1.THESE/MODELE_DYNAMIQUE_LB/CODE_MODELE_AK")
source("00.Packages_Chemins.R",encoding = "utf-8")
source("01.Fonctions_Simulation_V3.R",encoding = "utf-8")
source("02.Fonctions_VEM_V3.R", encoding = "utf-8")
source("04.Fonctions_graphique.R", encoding = "utf-8")
source("03.Fonctions_evaluations_V3.R", encoding = "utf-8")
source("05.Modeles_ref_V3.R", encoding = "utf-8")

### SIMULATION D'UN JEU DE DONNEES 

n = 100
nsignaux = 100
d=1
p=2

# Données météorologique : Facteurs exogènes à partir de données centrée et réduites
METEO <-  read.csv2(file="C:/Users/bonfils/Documents/1.THESE/0.DONNEES/METEO/table_meteo.csv")%>%
  mutate(Heure = substr(Date2, 12,16), Date = substr(Date2, 1, 10)) %>%
  mutate(Date = as.Date(Date, format="%d/%m/%Y")) %>%
  dplyr::select(-Date2) %>%
  dplyr::select(Heure, Date, TemperatureC2, Humidity2) %>%
  rename('Temperature' = "TemperatureC2", "Humidite"="Humidity2")%>%
  arrange(Date)%>%
  dplyr::select(-Date, -Heure) %>%
  mutate(Temperature = (Temperature - mean(Temperature))/sd(Temperature), 
         Humidite = (Humidite - mean(Humidite))/sd(Humidite))

# Construction des vecteurs u_t : 

d=1
p=2
U <- list()
for(t in 1:(nsignaux)){
  ut <- matrix(nrow = d, ncol=(d*(p+1)))
  ut[,1:d] <- rep(1,d)
  for(j in 2:(p+1)){
    ut[,((j-1)*d + 1) : (j*d)] <- diag(METEO[t,(j-1)], d)
  }  
  
  U[[t]] <- ut
}


# Paramêtres pour la simulation: 

K=2

Sigma_0 <- list()
Sigma_0[[1]] <- diag(1, d)
Sigma_0[[2]] <- diag(1, d)

Phi_k <- list()
Phi_k[[1]] <- diag(0.9, d)
Phi_k[[2]] <- diag(0.9, d)


K <- 2
d <- 1 
facteurs <- U

V_k <- list(diag(2,d), diag(2,d))
W_k <- list(diag(0.1,d), diag(0.1,d))
pi_simu <- c(0.5, 0.5)
mu_0 <- matrix(c(6.5, 6.5), nrow=d, ncol=K)
Sigma_0 <- Sigma_0
coefficients_Ak <- list(c(0, 3, 1), c(0, 1.1, 2.4))
coef_ecartement_centres <- 0.4 # Permet de générer des centres de classes plus ou moins éloignés
coef_melange <- 0.1 # Permet de générer 

BDD <- simulation_data_VEM_delta_melange_V3(V_k, W_k,
                                           Sigma_0, Phi_k,
                                           mu_0, pi_simu,
                                           facteurs,
                                           coefficients_Ak, nsignaux, K=2, 
                                           coef_ecartement_centres, coef_melange )






### ESTIMATION DES PARAMETRES DU MODELE 
res <- algorithme_VEM_AK(BDD, K=2, TT=nsignaux) 


### Résultats 


resultat_algorithme_VEM(BDD, res$best_solution)

graphique_matrice_classes(res$best_solution, BDD)

