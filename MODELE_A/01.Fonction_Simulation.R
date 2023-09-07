

simulation_data_VEM_delta_2_melange <- function(V_k, W_k, Sigma_0, Phi_k, mu_0, pi_simu, 
                                                facteur_exo, coefficients_A, nb_obs,
                                                K, delta, degres_melange ){
  ## Récupération des différentes valeurs : 
  d <- dim(facteur_exo[[1]])[1]  # Nombre de dimension 
  p <- (dim(facteur_exo[[1]])[2]/d)  # Nombre de facteur exo  
  n <- nb_obs # Nombre d'observations
  nsignaux <- length(facteur_exo)  # Nombre de pas de temps 
  
  if(length(pi_simu)!= K | length(V_k) != K | length(W_k) != K |
     length(Sigma_0) != K){
    print("At least one element has the wrong size ! ")
    return(NULL)
  }
  
  
  #Verifier que les proportions de classe sont ok : 
  
  if(sum(pi_simu) != 1){
    print("The sum of the given class proportions is not equal to 1 !")
    print(" Default proportion has been set : 1/K")
    pi_simu = rep(1/K, K)
  }
  ## Processus b 
  
  series_b <- list()
  
  b_0 <- sapply(seq(1,K), function(k,d) if(d!=1){
    mvrnorm(1,mu_0[,k], Sigma_0[[k]])
  }else{   rnorm(1,mu_0[k], Sigma_0[[k]])},d)
  
  # vAR
  
  bruit <- lapply(seq(1,K), function(k, d, nsignaux, var) mvrnorm(n = nsignaux,
                                                                  mu=rep(0,d), 
                                                                  Sigma = var[[k]]), d, nsignaux, W_k)
  
  
  construction_processus <- function(k, bruit, Phi_k, b_0){
    
    b_k <- matrix(nrow=nsignaux+1, ncol=d)
    
    b_k[1, ] <- matrix(b_0,nrow=d, ncol=K)[,k]
    for(t in 2:(nsignaux+1)){
      b_k[t, ] <-  Phi_k[[k]]%*%b_k[(t-1),] + bruit[[k]][t-1,]
    }
    return(b_k)
  }
  
  series_b <-  lapply(seq(1,K), construction_processus, bruit, Phi_k, b_0)
  bkt1 <- series_b[[1]]
  bkt2 <- series_b[[2]]
  
  bkt2_bis <- bkt2 + delta*(bkt2-bkt1)
  
  series_b <- list(bkt1, bkt2_bis)
  
  series_b <- centrage_bkt_simulation(series_b, K=K, nsignaux = nsignaux)
  
  ## Données observées X 
  non_empty = F
  while(non_empty == F){
    z <- numeric(n)
    for(i in 1:n){
      z[i] <- which(rmultinom(1,1,pi_simu)==1)
    }
    if(length(unique(z))==K){
      non_empty = TRUE
    }else{
      print("At least one class is empty")
    }
  }
  

  
  #3- Bruitage et normalisation  
  
  t_ik <- matrix(nrow=n, ncol= K) %>%
    as.data.frame %>%
    mutate(classes= z) 
  for(k in 1:K){
    t_ik[which(t_ik$classes == k),k] <- 1
    t_ik[which(t_ik$classes != k),k] <- 0
  }
  
  
  if(degres_melange>=0.5){
    print("attention le dégré de mélange est trop élevé --> Choisir une valent < 0.5")
    degres_melange=0.45
  }
  
  min_random <- degres_melange- 2*(degres_melange/3)
  
  t_ik<-  t_ik %>%
    dplyr::select(-classes)%>%
    apply( c(1,2), function(x) x +  runif(1, min=min_random, max=degres_melange)) %>%
    normalisation_poids()
  
  X <- matrix(nrow=(nsignaux*d), ncol=n)
  REFERENCE <- data.frame(DIM = rep(seq(1,d), nsignaux), TEMPS = 
                            unlist(lapply(seq(1,nsignaux), function(x) rep(x, d)))) 
  f3 <- function(i,t,d){
    #i: individus,
    #d : nombrre de dimension
    #t: temps
    if(d==1){
      bruit <- rnorm(n=1)
    }else{
      bruit <- mvrnorm(n=1, mu=rep(0,d), Sigma=diag(d))
    }
    
    x_it <- numeric(d)
    
    x_it <- facteur_exo[[t]]%*%coefficients_A + t_ik[i,1]*series_b[[1]][t+1,] +  t_ik[i,2]*series_b[[2]][t+1,] + V_k[[z[i]]]%*%bruit
    
    X[( (t-1)*d+1 ) : (t*d), i] <- x_it
    
    return(x_it)
  }
  
  for(t in 1:nsignaux){
    X[( (t-1)*d+1 ) : (t*d), ] <-sapply(seq(1,n), f3, t=t, d=d)
  }
  
  #####
  
  res_simu <- list(data_X = X, 
                   REFERENCE = REFERENCE, 
                   data_facteurs = facteur_exo,
                   data_processus =series_b, 
                   variances_X = V_k, 
                   variances_processus = W_k, 
                   classes= z , 
                   melange = t_ik,
                   coefficient= coefficients_A, 
                   Sigma_0 = Sigma_0, 
                   mu_0 = mu_0, 
                   Phi_k = Phi_k, 
                   pi_k = pi_simu)
  
  return(res_simu)
  
}



simulation_data_VEM_delta_4_melange <- function(V_k, W_k, Sigma_0, Phi_k, mu_0, pi_simu, 
                                                facteur_exo, coefficients_A,
                                                nb_obs, K, delta, degres_melange ){
  ## Récupération des différentes valeurs : 
  d <- dim(facteur_exo[[1]])[1]  # Nombre de dimension 
  p <- (dim(facteur_exo[[1]])[2]/d)  # Nombre de facteur exo  
  n <- nb_obs # Nombre d'observations
  nsignaux <- length(facteur_exo)  # Nombre de pas de temps 
  
  if(length(pi_simu)!= K | length(V_k) != K | length(W_k) != K |
     length(Sigma_0) != K){
    print("At least one element has the wrong size ! ")
    return(NULL)
  }
  
  
  #Verifier que les proportions de classe sont ok : 
  
  if(sum(pi_simu) != 1){
    print("The sum of the given class proportions is not equal to 1 !")
    print(" Default proportion has been set : 1/K")
    pi_simu = rep(1/K, K)
  }
  ## Processus b 
  
  series_b <- list()
  
  b_0 <- sapply(seq(1,K), function(k,d) if(d!=1){
    mvrnorm(1,mu_0[,k], Sigma_0[[k]])
  }else{   rnorm(1,mu_0[k], Sigma_0[[k]])},d)
  
  # vAR
  
  bruit <- lapply(seq(1,K), function(k, d, nsignaux, var) mvrnorm(n = nsignaux,
                                                                  mu=rep(0,d), 
                                                                  Sigma = var[[k]]), d, nsignaux, W_k)
  
  
  construction_processus <- function(k, bruit, Phi_k, b_0){
    
    b_k <- matrix(nrow=nsignaux+1, ncol=d)
    
    b_k[1, ] <- matrix(b_0,nrow=d, ncol=K)[,k]
    for(t in 2:(nsignaux+1)){
      b_k[t, ] <-  Phi_k[[k]]%*%b_k[(t-1),] + bruit[[k]][t-1,]
    }
    return(b_k)
  }
  series_b <-  lapply(seq(1,K), construction_processus, bruit, Phi_k, b_0)
  
  bkt1 <- series_b[[1]]
  bkt2 <- series_b[[2]]
  bkt3 <- series_b[[3]]
  bkt4 <- series_b[[4]]
  
  centres_b <- (1/4)*(bkt1 + bkt2 + bkt3 + bkt4)
  bkt1_bis <- bkt1 + delta*(bkt1 - centres_b)
  bkt2_bis <- bkt2 +delta*(bkt2-centres_b)
  bkt3_bis <- bkt3 + delta*(bkt3-centres_b) 
  bkt4_bis <- bkt4 + delta*(bkt4-centres_b) 
  
  
  series_b <- list(bkt1_bis, bkt2_bis, bkt3_bis, bkt4_bis)
  series_b <- centrage_bkt_simulation(series_b, K=K, nsignaux = nsignaux)
  
  ## Données observées X 
  non_empty = F
  while(non_empty == F){
    z <- numeric(n)
    for(i in 1:n){
      z[i] <- which(rmultinom(1,1,pi_simu)==1)
    }
    if(length(unique(z))==K){
      non_empty = TRUE
    }else{
      print("At least one class is empty")
    }
  }

  
  #3- Bruitage et normalisation  
  
  t_ik <- matrix(nrow=n, ncol= K) %>%
    as.data.frame %>%
    mutate(classes= z) 
  for(k in 1:K){
    t_ik[which(t_ik$classes == k),k] <- 1
    t_ik[which(t_ik$classes != k),k] <- 0
  }
  
  
  if(degres_melange>=0.5){
    print("attention le dégré de mélange est trop élevé --> Choisir une valeur < 0.5")
    degres_melange=0.45
  }
  
  min_random <- degres_melange- 2*(degres_melange/3)
  
  t_ik<-  t_ik %>%
    dplyr::select(-classes)%>%
    apply( c(1,2), function(x) x +  runif(1, min=min_random, max=degres_melange)) %>%
    normalisation_poids()
  
  X <- matrix(nrow=(nsignaux*d), ncol=n)
  REFERENCE <- data.frame(DIM = rep(seq(1,d), nsignaux), TEMPS = 
                            unlist(lapply(seq(1,nsignaux), function(x) rep(x, d)))) 
  f3 <- function(i,t,d){
    #i: individus,
    #d : nombrre de dimension
    #t: temps
    if(d==1){
      bruit <- rnorm(n=1)
    }else{
      bruit <- mvrnorm(n=1, mu=rep(0,d), Sigma=diag(d))
    }
    
    x_it <- numeric(d)
    
    x_it <- facteur_exo[[t]]%*%coefficients_A + t_ik[i,1]*series_b[[1]][t+1,] + t_ik[i,2]*series_b[[2]][t+1,] + t_ik[i,3]*series_b[[3]][t+1,] + t_ik[i,4]*series_b[[4]][t+1,] + V_k[[z[i]]]%*%bruit
    
    X[( (t-1)*d+1 ) : (t*d), i] <- x_it
    
    return(x_it)
  }
  
  for(t in 1:nsignaux){
    X[( (t-1)*d+1 ) : (t*d), ] <-sapply(seq(1,n), f3, t=t, d=d)
  }
  
  #####
  
  res_simu <- list(data_X = X, 
                   REFERENCE = REFERENCE, 
                   data_facteurs = facteur_exo,
                   data_processus =series_b, 
                   variances_X = V_k, 
                   variances_processus = W_k, 
                   classes= z , 
                   melange = t_ik,
                   coefficient= coefficients_A, 
                   Sigma_0 = Sigma_0, 
                   mu_0 = mu_0, 
                   Phi_k = Phi_k, 
                   pi_k = pi_simu)
  
  return(res_simu)
  
}
centrage_bkt_simulation <- function(series_b, nsignaux, K){
  
  # Centrage de l'ensemble des centres de classes
  mb=matrix(0,K,d)
  for (k in 1:K)
  {
    mb[k,] =colMeans(series_b[[k]])
  }
  
  mmb=colMeans(mb)
  
  
  for (k in 1:K)
  {
    series_b[[k]] = series_b[[k]]- matrix(1,nsignaux+1,1)%*%mmb
  }
  
  print("centered")
  
  
  return(series_b)
  
}







