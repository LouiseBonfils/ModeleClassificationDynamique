###### MODELE REF 1 : Algorithme EM : 
## x_it = a.ut + bk.zkt 
library(Rmixmod)



etape_E <- function(X, centres_init, mat_V_k, pi_k){
  n <- dim(X)[1]
  d <- dim(X)[2]
  K <- dim(centres_init)[1]
  
  #----------
  
  # A: Mettre à jours les probabilités d'appartenances : 
  t_ik <- matrix(nrow=n, ncol=K) # Matrice des probabilités d'appartenances qui sont recalculées au cours de l'étape E 
  logfikl <- matrix(0, nrow=n, ncol=K)
  
  for(k in 1:K){
    Y <- ((X-matrix(rep(centres_init[k,], n), nrow=n, ncol=d, byrow = T))^2)/diag(mat_V_k[[k]])
    logfikl[,k] <- - 0.5*as.matrix(rowSums(Y), nrow=n, ncol=1) + matrix(log(pi_k[k]), nrow=n, ncol=1) - (d/2)*matrix(log(2*pi), nrow=n, ncol=1) -0.5*matrix(sum(log(diag(mat_V_k[[k]]))), nrow=n, ncol=1)
  }
  
  #print(logfikl)
  
  t_ik <- round(exp(logfikl - matrix(rep(rowMax(logfikl),K), nrow=n, ncol=K)),4)
  t_ik <- normalisation_poids(t_ik) 
  #print(t_ik)
  
  #---- Verifier que la somme des deux fait bien 1 
  # if(sum(rowSums(t_ik))==n){
  #   # print("----------------------------------------------------------")
  #   # print("Les probabilités d'appartenances ont correctement été mises à jours")
  #   # print("----------------------------------------------------------")
  # }
  # 
  # B: Calculer la log-vraisemblance : 
  L= sum(matrix(rowMax(logfikl), ncol=1) + matrix(log(rowSums(exp(logfikl - matrix(rep(rowMax(logfikl),K), nrow=n, ncol=K)))),ncol=1))
  
  return(list(tik=t_ik, LogL=L, pi_k=pi_k, mat_V_k=mat_V_k, centres= centres_init))
}

etape_m <- function(X, t_ik, pi_k, mat_V_k, centres_init){
  
  n= dim(X)[1]
  d = dim(X)[2]
  K= dim(t_ik)[2]
  
  pi_k_new <- numeric(K)
  mat_V_k_new <- list()
  centres_init_new <- matrix(nrow=K, ncol=d)
  
  for(k in 1:K){
    poids_k <- t_ik[,k]
    
    # mise à jours des pi_k :
    pi_k_new[k] <- 1/K
    
    # mise à jours des centres de classes :  
    centres_init_new[k,] <- colSums(t_ik[,k]*X)/sum(t_ik[,k])
    
    # Mise à, jours des matrices de variances covariances
    A <- t(X-matrix(rep(centres_init_new[k,], n), nrow=n, ncol=d, byrow = TRUE))
    #%*%diag(t_ik[,k])%*%t(X-matrix(rep(centres_init_new[k,], n), nrow=n, ncol=d, byrow = T))
    B <- diag(t_ik[,k])
    C <- as.matrix((X-rep(centres_init_new[k,], n)))
    #, nrow=n, ncol=d, byrow = T)
    
    M <- A%*%B%*%C
    if(d!=1){
      mat_V_k_new[[k]] <- diag(diag(M/sum(poids_k)))
    }else{
      mat_V_k_new[[k]] <- M/sum(poids_k)
    }
    
  }
  return(list(pi_k=pi_k_new, centres_init = centres_init_new, mat_V_k= mat_V_k_new)) 
}



# Algorithme EM 

iteration_EM <- function(X, res_E){
  
  iter=0 
  fin_iter = F
  lik = res_E$LogL
  epsilon = 0.00001 
  l_past = res_E$LogL
  
  
  while((!fin_iter)&&(iter<200)){ 
    print(iter)
    iter= iter+1 
    
    param_new <- etape_m( X, t_ik=res_E$tik, mat_V_k = res_E$mat_V_k, pi_k=res_E$pi_k, centres_init = res_E$centres)
    
    # print("Les poids estimés :")
    # print(param_new$pi_k)
    # print('---------')
    # print("Les matrices de variances covariances : ")
    # print(param_new$mat_V_k)
    # print('---------')
    # print("Les centres de classes estimés :")
    # print(param_new$centres_init)
    # print('---------')
    # print('---------')
    
    res_E <- etape_E(X, param_new$centres_init, param_new$mat_V_k, param_new$pi_k)
    
    lik <- cbind(lik, res_E$LogL)
    
    #print(res_E$LogL)
    
    l_current = res_E$LogL
    
    if(is.na(res_E$LogL)){
      fin_iter <- T
      
    }
    
    if((!is.na(res_E$LogL)) && (abs(l_past-l_current) <= epsilon)){
      fin_iter <- T
    }
    print(abs(l_past-l_current))
    
    l_past <- res_E$LogL
    #print(fin_iter)
  }
  # classes <- numeric(0)
  # for(i in 1:n){
  #   classes[i] <- which.max(resultat_e$tik[i,])
  #   
  # }
  # 
  # classes <- resultat_e$tik 
  # 
  solution_EM <-list(parametre=param_new, Max_likelihood= l_past, vec_lik = lik, t_ik = res_E$tik) 
  
  return(solution_EM)
}


initialisation_EM <- function(X, K){
  n <- dim(X)[1]
  d <- dim(X)[2]
  
  best_Max_Likelihood <- -1000000000000
  nb_try = 0 
  
  mat_V_k_0 <- list() # Liste des matrices de variances covariances des observations x_it, V_k
  pi_k_0 <- numeric(K)
  
  centres_0 <- X[c(sample(seq(1,n),K)),] # Matrice avec les centres de classes initiaux avec les k vecteurs en colones 
  
  centres_0 <- matrix(centres_0, nrow=K, ncol=d)
  
  for(k in 1:K){
    
    mat_V_k_0[[k]] <- diag(1,d)
    #diag(apply(X, 2, var))
    pi_k_0[k] <- 1/K
  }
  
  first_iter <- etape_E(X, centres_init = centres_0, mat_V_k = mat_V_k_0, pi_k = pi_k_0)
  
  return(first_iter)
}


algorithm_EM <- function(X, K, try_max){
  
  n= dim(X)[1]
  d=dim(X)[2]
  
  nb_try <- 0
  nb_good_try <- 0 
  best_likelihood <- -100000000000
  best_solution <- NA 
  
  good_solution <- list()
  
  while(nb_good_try<= try_max & nb_try <2){
    #print("algoe em !!")
    nb_try <- nb_try + 1 
    
    # Initialisation 
    
    #print('Initialisation...')
    initialisation <- initialisation_EM(X,K)
    
    # Itération de l'algorithme 
    #print("Iteration...")
    solution <- iteration_EM(X, initialisation)
    
    # Comparaison des solutions 
    if((!is.na(solution$Max_likelihood)) && (solution$Max_likelihood > best_likelihood)){
      best_likelihood <- solution$Max_likelihood
      best_solution <- solution
      good_solution <- append(good_solution, solution)
      
    }else if((!is.na(solution$Max_likelihood)) && (solution$Max_likelihood <= best_likelihood)){
      nb_good_try <- nb_good_try +1
      good_solution <- append(good_solution, solution)
    }
  }
  
  # Attribution des classes : 
  
  best_solution$cluster <- apply(best_solution$t_ik, 1, which.max)
  
  # Calcul du BIC 
  
  best_solution$BIC <- 2*best_solution$Max_likelihood - (2*K*d + K -1)*log(n)
  
  return(best_solution)
}

estimation_coef_A <- function(BDD){
  
  d <- dim(BDD$data_facteurs[[1]])[1]
  K <- length(BDD$data_processus)
  nsignaux <- dim(BDD$data_X)[1]
  n <- dim(BDD$data_X)[2]
  dim_U <- dim(BDD$data_facteurs[[1]])[2]
  
  A <- numeric(dim(BDD$data_facteurs[[1]])[2])
  
  
  X <- BDD$data_X
  REFERENCE <- BDD$REFERENCE
  
  s2 <- 0
  s1 <- 0
  V_k <- diag(1,d)
  
  for(t in 1:nsignaux){
    a2 <- t(BDD$data_facteurs[[t]]) %*% solve(V_k) 
    a1 <- t(BDD$data_facteurs[[t]]) %*% solve(V_k) %*% BDD$data_facteurs[[t]]  
    
    for(i in 1:n){
      b2 <- (X[which(REFERENCE$TEMPS==t),i])
      s2 <- s2 + (a2%*%b2)
      s1 <- s1 + a1
    } 
  } 
  
  coefficient_A <- solve(s1)%*%s2
  
  return(coefficient_A)
  
}



estimation_A_retraitement <- function(simulation){
  coeff_A <- estimation_coef_A(simulation)
  
  X <- simulation$data_X
  REFERENCE <- simulation$REFERENCE
  nsignaux <- length(simulation$data_facteurs)
  X_bis <- matrix(nrow = dim(X)[1], ncol= dim(X)[2])
  
  for(i in 1:dim(X)[2]){
    for(t in 1:nsignaux){
      
      X_bis[which(REFERENCE$TEMPS==t), i] <-  X[which(REFERENCE$TEMPS==t), i] -  simulation$data_facteurs[[t]]%*%coeff_A
    }
  }
  
  simulation$data_X <- X_bis
  simulation$coeff_estim <- coeff_A
  
  
  return(simulation)
}




algorithme_EM_CST <- function(simulation, K, nsignaux){
  # 1 Estimation A 
  
  simulation_EM <- estimation_A_retraitement(simulation)
  
  nb_obs <- dim(simulation$data_X)[2]
  facteur_exo <- simulation$data_facteurs
  
  # Mélange 
  
  simulation_EM_X <- simulation_EM$data_X %>%
    as.data.frame %>%
    mutate(temps = seq(1, nsignaux)) %>%
    gather(key='indiv', value="X", -temps) %>%
    mutate(indice = seq(1, nb_obs*nsignaux))
  
  X <- simulation_EM_X %>%
    dplyr::select(X)
  
  new_ref <- simulation_EM_X
  
  res_EM <-  mixmodCluster(X, nbCluster = K)
  
  
  # Return les bonnes informations 
  
  solution <- list()
  n <- dim(simulation$data_X)[2]
  d <- 1
  nb_var <- dim(facteur_exo[[1]])[2]/d
  nbClasse <- K
  solution$size <- list()
  solution$size$d <- d
  solution$size$n <- n
  solution$size$nsignaux <- nsignaux
  solution$size$nbClasse <- K
  
  #classe 
  solution$parametre <- list()
  solution$parametre_var <- list()
  
  solution$parametre$sigma <- res_EM@bestResult@parameters@variance
  
  table_indiv_classe <- merge(new_ref, res_EM@results[[1]]@proba%>%
                                as.data.frame %>%
                                mutate(indice = seq(1, dim(X)[1])))
  
  classes <- table_indiv_classe %>%
    mutate(classes = res_EM@results[[1]]@partition)
  
  classes_X <- table(indiv = classes$indiv, classes = classes$classes) %>%
    as.data.frame %>%
    pivot_wider(values_from = Freq, names_from = classes) %>%
    column_to_rownames('indiv')
  
  classes_X$classes <- apply(classes_X, 1, function(x) which.max(x)) 
  
  
  classes_XX <-  merge(table_indiv_classe, classes_X%>%
          dplyr::select(classes) %>%
          rownames_to_column('indiv')) %>%
    mutate(indiv = as.numeric(str_replace(indiv, 'V', '')) )%>%
    dplyr::select(indiv, classes) %>%
    group_by(indiv) %>%
    summarise(classes = mean(classes)) %>%
    arrange(indiv)
  
  
  t_ik <- matrix(nrow=n, ncol= K) %>%
    as.data.frame %>%
    mutate(classes= classes_XX$classes) 
  
  for(k in 1:K){
    t_ik[which(t_ik$classes == k),k] <- 1
    t_ik[which(t_ik$classes != k),k] <- 0
  }
  
  
  
  
  solution$parametre_var$t_ik <- t_ik%>%
    dplyr::select(-classes) %>%
    as.matrix
    
    # table_indiv_classe %>%
    # gather(key='classe', value='proba', -X, -temps, -indice, -indiv)%>%
    # group_by( indiv, classe) %>%
    # summarise(proba = mean(proba)) %>%
    # as.data.frame %>%
    # pivot_wider(values_from = proba, names_from = classe)%>%
    # dplyr::select(-indiv) %>%
    # as.matrix
    # 
  
  
  
  #Centres 
  
  solution$filtre_Kalman$mu_kt_B <- list()
  for(k in 1:K){
    solution$filtre_Kalman$mu_kt_B[[k]]<- matrix(res_EM@results[[1]]@parameters@mean[k], nrow=(nsignaux+1), ncol=d)
    
  }
  
  #coefficient 
  solution$coefficient <- simulation_EM$coeff_estim
  
  res <- list()
  res$best_solution <- solution
  return(res)
}





algorithme_REG_KMEANS <- function(simulation, K, nsignaux){
  # 1 Estimation A 
  
  simulation_KMEANS <- estimation_A_retraitement(simulation)
  
  facteur_exo <- simulation$data_facteurs
  
  ###############################################
  
  res_KMEANS <-  kmeans(t(simulation_KMEANS$data_X), centers = K, nstart = 30)
  
  
  # Return les bonnes informations 
  
  solution <- list()
  n <- dim(simulation$data_X)[2]
  d <- 1
  nb_var <- dim(facteur_exo[[1]])[2]/d
  nbClasse <- K
  solution$size <- list()
  solution$size$d <- d
  solution$size$n <- n
  solution$size$nsignaux <- nsignaux
  solution$size$nbClasse <- K
  
  #classe 
  solution$parametre <- list()
  solution$parametre_var <- list()
  
  
  #### Proba 
  
  t_ik <- matrix(nrow=n, ncol= K) %>%
    as.data.frame %>%
    mutate(classes= res_KMEANS$cluster) 
  
  for(k in 1:K){
    t_ik[which(t_ik$classes == k),k] <- 1
    t_ik[which(t_ik$classes != k),k] <- 0
  }
  
  # set.seed(45)
  t_ik <-  t_ik %>%
    dplyr::select(-classes)
  #   apply( c(1,2), function(x) x +  runif(1, min=0, max=0.05)) %>%
  #   normalisation_poids()
  # set.seed(NULL)
  
  solution$parametre_var$t_ik <- t_ik
  
  
  #Centres 
  
  solution$filtre_Kalman$mu_kt_B <- list()
  for(k in 1:K){
    solution$filtre_Kalman$mu_kt_B[[k]]<- matrix(c(res_KMEANS$centers[k,1], res_KMEANS$centers[k,]), nrow=(nsignaux+1),ncol=1)
    
  }
  
  #coefficient 
  solution$coefficient <- simulation_KMEANS$coeff
  
  
  res <- list()
  res$best_solution <- solution
  return(res)
}
