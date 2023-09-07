initialisation_VEM_Kmeans <- function(X, K, facteur_exo, REFERENCE, nsignaux){
  
  deb_init <- tic()
  
  solution <- list()
  n <- dim(X)[2]
  d <- 1
  nb_var <- dim(facteur_exo[[1]])[2]/d
  nbClasse <- K
  solution$size <- list()
  solution$size$d <- d
  solution$size$n <- n
  solution$size$nsignaux <- nsignaux
  solution$size$nbClasse <- K
  
  
  
  ## INITIALISATION DES PARAMETRE
  
  variances_Vk<- list() 
  variances_Wk <- list() 
  
  for(k in 1:K){
    variances_Wk[[k]] <- diag(1,d)
    variances_Vk[[k]] <- diag(1,d)
  }
  
  solution$parametre$variances_W_k <- variances_Wk 
  solution$parametre$variances_V_k <- variances_Vk 
  
  pi_k <- rep(1/K, K)
  
  ## Filtre de Kalman 
  mu_kt_F <- list()
  mu_kt_B <- list()
  
  P_kt_F <- list()
  P_kt_B <- list()
  
  for(k in 1:K){
    P_kt_F[[k]] <- list()
    P_kt_B[[k]] <- list()
  }
  
  for(t in 1:(nsignaux+1)){
    for(k in 1:K){
      P_kt_F[[k]][[t]] <- diag(1, d)
      P_kt_B[[k]][[t]] <- diag(1, d)
    }
  }
  
  
  data <- X 
  
  ## MISE A JOUR DES PROBABILITES D'APPARTENANCE : 
  
  data_k_means <- t(X)
  
  #2- Classification k-means
  
  res_k_means <- kmeans(data_k_means, K, iter.max=20, nstart=15)
  classes_initiales <- res_k_means$cluster
 
  
  solution$parametre$pi_k <- rep(1/K, K)
  
  #3- Bruitage et normalisation  
  
  t_ik <- matrix(nrow=n, ncol= K) %>%
    as.data.frame %>%
    mutate(classes= classes_initiales) 
  for(k in 1:K){
    t_ik[which(t_ik$classes == k),k] <- 1
    t_ik[which(t_ik$classes != k),k] <- 0
  }
  set.seed(45)
  t_ik <-  t_ik %>%
    dplyr::select(-classes)%>%
    apply( c(1,2), function(x) x +  runif(1, min=0, max=0.3)) %>%
    normalisation_poids()
  set.seed(NULL)
  
  mu_k <- list() 
  
  classes_initiales <- res_k_means$cluster
  
  
  for(k in 1:K){
    
    mu_k[[k]] <- matrix(0, nrow=nsignaux, ncol= d)
  }
  
  
  iter= 0 
  
  coefficient_A <- list()  
  
  
  
  while(iter<2){
    
    for(k in 1:K){
      INDICE <- merge(data.frame(individu=seq(1,n)), data.frame(temps=seq(1,nsignaux)))
      
      
      if(d*nb_var==1){
        
        t2 <- sum(mapply(function(i,t)
          t_ik[i,k]* t(facteur_exo[[t]]) %*% solve(variances_Vk[[k]])%*%(X[which(REFERENCE$TEMPS==t),i] - mu_k[[k]][t,]), 
          INDICE$individu, INDICE$temps, SIMPLIFY = "array"))
        
        t1 <- sum(mapply(function(i,t) t_ik[i,k]* t(facteur_exo[[t]]) %*% solve(variances_Vk[[k]]) %*% facteur_exo[[t]], 
                         INDICE$individu, INDICE$temps, SIMPLIFY = "array"))
        
        
        
      }else{
        t2 <- rowSums(mapply(function(i,t)
          t_ik[i,k]* t(facteur_exo[[t]]) %*% solve(variances_Vk[[k]])%*%(X[which(REFERENCE$TEMPS==t),i] - mu_k[[k]][t,]), 
          INDICE$individu, INDICE$temps, SIMPLIFY = "array"), dims=2)
        
        t1 <- rowSums(mapply(function(i,t) t_ik[i,k]* t(facteur_exo[[t]]) %*% solve(variances_Vk[[k]]) %*% facteur_exo[[t]], 
                             INDICE$individu, INDICE$temps,  SIMPLIFY = "array"), dims=2)
        
        
      }
      
      coefficient_A[[k]] <- solve(t1)%*%t2
    }
    
    #### Initialisation des centres de classes 
    
    
    X_bis <- matrix(nrow = dim(X)[1], ncol= dim(X)[2])
    for(i in 1:n){
      for(t in 1:nsignaux){
        if(d==1){
          X_bis[which(REFERENCE$TEMPS==t), i] <-  sum(sapply(seq(1,K), function(k, i, t) 
            t_ik[i,k]*(X[which(REFERENCE$TEMPS==t), i] - facteur_exo[[t]]%*%coefficient_A[[k]]),i,t))
          
        }else{
          X_bis[which(REFERENCE$TEMPS==t), i] <-  rowSums(sapply(seq(1,K), function(k, i, t) 
            t_ik[i,k]*(X[which(REFERENCE$TEMPS==t), i] - facteur_exo[[t]]%*%coefficient_A[[k]]),i,t))
          
        }
      }
    }
    
    
    data_k_means <- t(X_bis)
    #2- Classification k-means
    
    res_k_means <- kmeans(t(X_bis), K, iter.max=20, nstart=15)
    classes_initiales <- res_k_means$cluster
    
    #3- Bruitage et normalisation  
    
    t_ik <- matrix(nrow=n, ncol= K) %>%
      as.data.frame %>%
      mutate(classes= classes_initiales) 
    
    for(k in 1:K){
      t_ik[which(t_ik$classes == k),k] <- 1
      t_ik[which(t_ik$classes != k),k] <- 0
    }
    
    set.seed(45)
    t_ik <-  t_ik %>%
      dplyr::select(-classes)%>%
      apply( c(1,2), function(x) x +  runif(1, min=0, max=0.3)) %>%
      normalisation_poids()
    set.seed(NULL)
    
    classes_initiales <- res_k_means$cluster
    
    for(k in 1:K){
      indice <- which(classes_initiales==k)
      XX <- X_bis[,indice]
      
      if(d != 1){
        mu_k[[k]] <- cbind(XX, REFERENCE) %>%
          gather(key='indiv', value="value", -DIM, -TEMPS)%>%
          group_by(DIM, TEMPS) %>%
          summarise(mm = mean(value)) %>%
          pivot_wider(names_from = TEMPS, values_from = mm) %>%
          as.data.frame%>%
          column_to_rownames(var="DIM")%>%
          t
      }else{ 
        mu_k[[k]] <- cbind(XX, REFERENCE) %>%
          gather(key='indiv', value="value", -DIM, -TEMPS)%>%
          group_by(TEMPS) %>%
          summarise(mm = mean(value)) %>%
          pivot_wider(names_from = TEMPS, values_from = mm) %>%
          as.data.frame%>%
          t
        
      }
      
    }
    
    
    # Centrage de l'ensemble des centres de classes
    if(K !=1){
      mb=matrix(0,K,d)
      for (k in 1:K)
      {
        mb[k,] =colMeans(mu_k[[k]])
      }
      
      mmb=colMeans(mb)
      
      for (k in 1:K)
      {
        mu_k[[k]] = mu_k[[k]]- matrix(1,nsignaux,1)%*%mb[k,]
      }
      
    }
    
    iter <- iter + 1
    
    
  }
  
  
  solution$coefficient <-coefficient_A
  
  solution$parametre_var$lambda_k <- rep(1,K)
  
  solution$parametre$mu_0 <- list()
  solution$parametre$Sigma_0 <- list()
  
  for(k in 1:K){
    if(d!=1){
      solution$parametre$mu_0[[k]] <- colMeans(mu_k[[k]][1:min(7, nsignaux),])
      
      solution$parametre$Sigma_0[[k]] <- diag(10,d)
    }else{
      solution$parametre$mu_0[[k]] <- mean(mu_k[[k]][1:min(7, nsignaux),])
      solution$parametre$Sigma_0[[k]] <- diag(10,d)
    }
    
  }
  
  
  for(k in 1:K){
    if(d!=1){
      mu_k[[k]] <- rbind(colMeans(mu_k[[k]][1:min(7, nsignaux),]), mu_k[[k]])
    }else{
      mu_k[[k]] <- rbind(mean(mu_k[[k]][1:min(7, nsignaux),]), mu_k[[k]])
    }
  }
  
  
  solution$parametre_var$t_ik <- t_ik
  
  solution$filtre_Kalman$mu_kt_F <- mu_k # Forward 
  solution$filtre_Kalman$mu_kt_B <- mu_k # Backward
  
  
  solution$filtre_Kalman$P_kt_F <- P_kt_F # Forward
  solution$filtre_Kalman$P_kt_B <- P_kt_B # Backward
  
  
  ## Initalisation Phi_k 
  
  solution$parametre$Phi_k <- list()
  
  for(k in 1:K){
    if(d!=1){
      M1 <- rowSums(sapply(seq(2, (nsignaux+1)), function(t,k) mu_k[[k]][t,]%*%t(mu_k[[k]][(t-1),]),k, simplify="array" ), dims=2)
      
      M2 <- rowSums(sapply(seq(2, (nsignaux+1)), function(t,k,d) mu_k[[k]][t,]%*%t(mu_k[[k]][(t),])+ diag(1,d) ,k,d , simplify="array" ), dims=2)
      
    }else{ 
      M1 <- sum(sapply(seq(2, (nsignaux+1)), function(t,k) mu_k[[k]][t,]%*%t(mu_k[[k]][(t-1),]),k, simplify="array" ))
      
      M2 <- sum(sapply(seq(2, (nsignaux+1)), function(t,k,d) mu_k[[k]][t,]%*%t(mu_k[[k]][(t),])+ diag(1,d) ,k,d , simplify="array" ))
      
    }
    
    solution$parametre$Phi_k[[k]] <- diag(diag(M1%*%solve(M2)),d)
    
  }
  
  ## Initialisation des poids t_ik à 0 pour pouvoir remplir ensuite 
  
  solution$parametre_var$t_ik <- t_ik
  
  fin_init <- toc(quiet=F)
  
  
  return(solution)
}


log_densite_normale <- function(x, mu, Sigma){

  d <- length(x)
  res <- -(d/2)*log(2*pi) - (1/2)*log(det(Sigma)) - (1/2)*(matrix((x-mu), nrow=1, ncol=d)%*%solve(Sigma)%*%matrix((x-mu), 
                                                                                                                  nrow=d, ncol=1))
  
  return(res)
}


trace <- function(M){
  if(dim(M)[1]!=dim(M)[2]){
    print('The matrix have to be a square matrix')
    return(NULL)
  }else{
    return(sum(diag(M))) 
  }
}


iterations_VEM <- function(simulation, K, TT, solution){
  
  table_res <- data.frame(erreur_centres =NA, erreur_facteur = NA)
  X <- simulation$data_X 
  REFERENCE <- simulation$REFERENCE
  facteur_exo <- simulation$data_facteurs
  
  d <- solution$size$d
  K <- solution$size$nbClasse 
  nsignaux <- solution$size$nsignaux
  n <- solution$size$n
  
  table_erreurs <- data.frame(iter=NA, erreur_centres = NA, erreur_facteur = NA, taux_classif = NA, ELBO = NA)
  
  # ITERATION 
  iter = 0
  log_vraisemblance <- NA
  ml_av <- -1000000
  fin_iter <- F
  debut <- tic()
  
  while(iter < 500 &  !fin_iter){
    dITER <- tic()
    ## MISE A JOUR DES CENTRES DE CLASSES ###
    c_kt <- list()
    m_kt <- list()
    
    for(k in 1:K){
      c_kt[[k]] <- matrix(nrow=(nsignaux+1), ncol=d)
      m_kt[[k]] <- matrix(nrow=(nsignaux+1), ncol=d)
    }
    
    P_kt_F <- solution$filtre_Kalman$P_kt_F
    P_kt_B <- solution$filtre_Kalman$P_kt_B
    
    for(k in 1:K){
      
      ## INIT
      c_kt[[k]][1,]<- solution$parametre$mu_0[[k]]
      P_kt_F[[k]][[1]] <- solution$parametre$Sigma_0[[k]]
      # FORWARD
      for(t in 2:(nsignaux+1)){
        
        if(d==1){
          ss_1 <- sum(sapply(seq(1,n), function(i,t,k)
            solution$parametre_var$t_ik[i,k]*(X[which(REFERENCE$TEMPS == (t-1)),i] - facteur_exo[[t-1]]%*%solution$coefficient[[k]] - solution$parametre$Phi_k[[k]]%*%matrix(c_kt[[k]][(t-1),],ncol=1)), t,k))
          
        }else{
          ss_1 <- rowSums(sapply(seq(1,n), function(i,t,k)
            solution$parametre_var$t_ik[i,k]*(X[which(REFERENCE$TEMPS == (t-1)),i]- facteur_exo[[t-1]]%*%solution$coefficient[[k]] - solution$parametre$Phi_k[[k]]%*%matrix(c_kt[[k]][(t-1),],
                                                                                                                                                                            ncol=1)), t,k))
        }
        
        
        P_kt_F[[k]][[t]] <- solve( solve(solution$parametre$Phi_k[[k]]%*%P_kt_F[[k]][[t-1]]%*%solution$parametre$Phi_k[[k]]+ solution$parametre$variances_W_k[[k]])+ solve(solution$parametre$variances_V_k[[k]])*(sum(solution$parametre_var$t_ik[,k])))
        
        c_kt[[k]][t,] <- solution$parametre$Phi_k[[k]]%*%c_kt[[k]][(t-1),] + P_kt_F[[k]][[t]]%*%solve(solution$parametre$variances_V_k[[k]])%*%ss_1
        
      }
    }
    
    
    ## Backward
    
    # Init
    J_kt <- list()
    
    for(k in 1:K){
      
      J_kt[[k]]<-list()
      
      J_kt[[k]][[nsignaux+1]] <- matrix(nrow=d, ncol=d)
      
      m_kt[[k]][(nsignaux+1),] <- c_kt[[k]][(nsignaux+1),]
      P_kt_B[[k]][[(nsignaux+1)]] <- P_kt_F[[k]][[(nsignaux+1)]]
      
      M <- solve(solution$parametre$Phi_k[[k]]%*%P_kt_F[[k]][[(nsignaux+1)]]%*%solution$parametre$Phi_k[[k]] + solution$parametre$variances_W_k[[k]])
      
      J_kt[[k]][[nsignaux+1]] <- P_kt_F[[k]][[(nsignaux+1)]]%*%solution$parametre$Phi_k[[k]]%*%M
      
      
      for(t in (nsignaux):1){
        
        M <- solve(solution$parametre$Phi_k[[k]]%*%P_kt_F[[k]][[(t)]]%*%solution$parametre$Phi_k[[k]] + solution$parametre$variances_W_k[[k]])
        
        
        J_kt[[k]][[t]] <- P_kt_F[[k]][[(t)]]%*%solution$parametre$Phi_k[[k]]%*%M
        
        m_kt[[k]][t,] <- c_kt[[k]][t,] + J_kt[[k]][[t]]%*%(m_kt[[k]][t+1,] - solution$parametre$Phi_k[[k]]%*%c_kt[[k]][t,])
        
        p_kt_moins <- solution$parametre$Phi_k[[k]]%*%P_kt_F[[k]][[t]]%*%t(solution$parametre$Phi_k[[k]]) + solution$parametre$variances_W_k[[k]]
        
        P_kt_B[[k]][[t]] <- P_kt_F[[k]][[t]] - J_kt[[k]][[t]] %*%(p_kt_moins - P_kt_B[[k]][[(t+1)]])%*%t(J_kt[[k]][[t]] )
        
        
        
      }
      
      
    }
    
    
    
    
    
    
    # Centrage de l'ensemble des centres de classes
    if(K !=1){
      mb=matrix(0,K,d)
      for (k in 1:K)
      {
        mb[k,] =colMeans(m_kt[[k]])
      }
      
      mmb=colMeans(mb)
      
      for (k in 1:K)
      {
        m_kt[[k]] = m_kt[[k]]- matrix(1,nsignaux+1,1)%*%mb[k,]
      }
      
    }
    
    
    solution$filtre_Kalman$P_kt_F <- P_kt_F
    solution$filtre_Kalman$P_kt_B <- P_kt_B
    
    solution$filtre_Kalman$mu_kt_F <- c_kt
    solution$filtre_Kalman$mu_kt_B <- m_kt 
    
    
    ### MISE A JOUR DES POIDS ###
    if(d==1){
      
      numerateur<- matrix(0,nrow=n,ncol=K) 
      
      for(k in 1:K){
        
        effet_ext <- numeric(nsignaux)
        
        for(t in 1:nsignaux){
          effet_ext[t] <- facteur_exo[[t]]%*%solution$coefficient[[k]]
        }
        
        
        a <- nsignaux*(log(solution$parametre$pi_k[k]) -0.5*solution$parametre_var$lambda_k[k]*trace(solve(solution$parametre$variances_V_k[[k]])))
        
        
        ##########
        
        mu <- solution$filtre_Kalman$mu_kt_B[[k]][-1]*matrix(1,nrow=solution$size$nsignaux,ncol=solution$size$n) + effet_ext
        
        
        Sigma <- solution$parametre$variances_V_k[[k]][1,1]
        
        XM <- X-mu
        XM_s <- XM^2
        RES <- log((1/sqrt(2*pi*Sigma))*exp((-1/(2*Sigma))*XM_s)) 
        
        numerateur[,k] <- colSums(RES) + a
        
      }
      
      
      M <-exp(numerateur-matrix(rep(rowMax(numerateur),K),nrow=n, ncol=K))
      MM <- apply(M, c(1, 2), function(x) max(x,exp(-50)))
      tt <- normalisation_poids(MM)
      
      if(sum(is.na(tt))==0){
        
        solution$parametre_var$t_ik <- tt
        
      }
      
      
      
      
    }else{

      solution <- mise_a_jours_poids_2(solution, facteur_exo, X, REFERENCE)
    }
    
    
    ### MISE A JOURS DES LAMBDA_K ###
    
    lambda_k<- numeric(K)
    
    for(k in 1:K){
      
      numerateur <- (nsignaux+1)*d
      a <- nsignaux*trace(solve(solution$parametre$variances_V_k[[k]]))*sum(solution$parametre_var$t_ik[,k])
      b <- nsignaux*trace(solve(solution$parametre$variances_W_k[[k]]))
      c <- nsignaux*trace(solve(solution$parametre$variances_W_k[[k]])%*%solution$parametre$Phi_k[[k]]%*%solution$parametre$Phi_k[[k]])
      e <- trace(solve(solution$parametre$Sigma_0[[k]]))
      denominateur <-  a+b+c+e
      
      lambda_k[k] <- numerateur/denominateur
      
      
    }
    solution$parametre_var$lambda_k <- lambda_k
    
    
    ### MISE A JOUR DES PI ###
    K <- solution$size$nbClasse 
    n <- solution$size$n
    
    solution$parametre$pi_k <- rep(1/K, K)
    
    ### MISE A JOUR DES MU_0 ET SIGMA_0 ###
    
    for(k in 1:K){
      solution$parametre$mu_0[[k]]<- solution$filtre_Kalman$mu_kt_B[[k]][1,]  
      solution$parametre$Sigma_0[[k]] <- diag(solution$parametre_var$lambda_k[k], d)
    }
    
    
    
    ### MISE A JOUR DE V_k ###
    Vk <- list()
    
    
    for(k in 1:K){
      
      effet_ext <- numeric(nsignaux)
      for(t in 1:nsignaux){
        effet_ext[t] <- facteur_exo[[t]]%*%solution$coefficient[[k]]
      }
      
      mu <- solution$filtre_Kalman$mu_kt_B[[k]][-1]*matrix(1,nrow=solution$size$nsignaux,ncol=solution$size$n) + effet_ext
      
      XM <- X - mu
      XM_s <- XM^2
      
      XM_s[t,1]
      solution$parametre_var$t_ik[1,k]*(XM_s[1,1] + solution$parametre_var$lambda_k[k])
      
      
      RES <- solution$parametre_var$t_ik[,k]*t(XM_s + solution$parametre_var$lambda_k[k])
      
      Vk[[k]] <- (1/(nsignaux*sum(solution$parametre_var$t_ik[,k])))*matrix(sum(RES), nrow=d, ncol=d)
      
    }
    
    solution$parametre$variances_V_k <- Vk
    
    
    
    ### MISE A JOUR DE W_K ###
    
    for(k in 1:K){
      
      mu <- solution$filtre_Kalman$mu_kt_B[[k]][-1,] - solution$parametre$Phi_k[[k]]%*%solution$filtre_Kalman$mu_kt_B[[k]][-(nsignaux+1),]
      
      M <- (1/nsignaux)*(sum(mu^2)) + solution$parametre_var$lambda_k[k]*(1+ solution$parametre$Phi_k[[k]]^2)
      
      solution$parametre$variances_W_k[[k]] <- matrix(M, nrow=d, ncol=d) 
      
    }
    
    ### MISE A JOUR DES PHI_K ###
    
    
    for(k in 1:K){
      
      somme_1 <- 0
      somme_2 <- 0 
      
      epsilon <- 0.0000165 
      
      for(t in 2:nsignaux+1){
        
        z1 <- solution$filtre_Kalman$mu_kt_B[[k]][t,] %*% t(solution$filtre_Kalman$mu_kt_B[[k]][t-1,])
        z2 <- solution$filtre_Kalman$mu_kt_B[[k]][t-1,] %*% t(solution$filtre_Kalman$mu_kt_B[[k]][t-1,]) + diag(solution$parametre_var$lambda_k[k],d)
        
        somme_1 <- somme_1 + z1
        somme_2 <- somme_2 + z2
        
      }
      
      M <- somme_1%*%solve(somme_2)
      
      phi_k <- diag(diag(M),d)
      
      if(sum(diag(phi_k)>=1)!=0){
        diag(phi_k)[abs(diag(phi_k))>=1] <- 1-epsilon
      }
      
      solution$parametre$Phi_k[[k]]<- phi_k
    }
    
    
    ### MISE A JOUR DES COEFFICIENTS A ###
    X <- simulation$data_X 
    REFERENCE <- simulation$REFERENCE
    facteur_exo <- simulation$data_facteurs
    d <- solution$size$d
    K <- solution$size$nbClasse 
    nsignaux <- solution$size$nsignaux
    n <- solution$size$n
    dim_U <- length(solution$coefficient[[1]])
    
    for(k in 1:K){
      
      X_M <- X-solution$filtre_Kalman$mu_kt_B[[k]][-1]*matrix(1,nrow=solution$size$nsignaux,
                                                              ncol=solution$size$n)
      
      somme_effect_square <- 0
      somme_effect <- rep(0, dim_U)
      for(t in 1:nsignaux){
        somme_effect <- somme_effect + rowSums(t(facteur_exo[[t]])%*%(solution$parametre_var$t_ik[,k]*t(X_M[t,])))
      }
      
      for(t in 1:nsignaux){
        somme_effect_square <- somme_effect_square + t(facteur_exo[[t]])%*%facteur_exo[[t]]
      }
      
      
      argument_1 <- matrix(0, nrow=dim_U, ncol=dim_U)
      argument_2 <- rep(0,dim_U)
      
      argument_1 <- sum(solution$parametre_var$t_ik[,k]*(1/solution$parametre$variances_V_k[[k]][1,1]))*somme_effect_square
      argument_2 <- (1/solution$parametre$variances_V_k[[k]][1,1])%*%somme_effect
      
      solution$coefficient[[k]] <-  solve(argument_1)%*%t(argument_2)
      
    }
    
   
    if(iter==0){
      classe_av = solution$filtre_Kalman$mu_kt_B
    }else{
      classe_new = solution$filtre_Kalman$mu_kt_B
      
      ss = 0
      for(k in 1:solution$size$nbClasse){
        
        ss= ss+ (mean((classe_av[[k]]- classe_new[[k]])^2))
        
      }
      if(ss <0.00001){
        fin_iter <- TRUE
      }
    }
    
    
    classe_av <-solution$filtre_Kalman$mu_kt_B   
    
    iter <- iter + 1

    
  }
  
  
  elbo <- compute_ELBO(solution, facteur_exo = facteur_exo,X = X, REFERENCE = REFERENCE)
  bb <- calcul_BIC(solution, elbo)
  
  fin <- toc(quiet = TRUE)
  
  
  temps_calcul <- fin$toc-debut
  
  
  return(list(BIC=bb, best_solution = solution, nbr_iteration = iter,
              temps_calcul= temps_calcul, elbo=elbo))
  
  
}

compute_ELBO <- function(solution, facteur_exo, X, REFERENCE){
  
  
  d <- solution$size$d
  K <- solution$size$nbClasse 
  nsignaux <- solution$size$nsignaux
  n <- solution$size$n
  
  ## Premier élément, t=0 
  log_mod_1 <- 0 
  log_mod_2 <- 0 
  log_mod_3 <- 0 
  log_mod_4 <- 0 
  
  
  for(k in 1:K){
    
    ## 1 
    z1 <- log_densite_normale(solution$filtre_Kalman$mu_kt_B[[k]][1,],
                              solution$parametre$mu_0[[k]], solution$parametre$Sigma_0[[k]])
    z2 <-  z1 -0.5*solution$parametre_var$lambda_k[k]*trace(solve(solution$parametre$Sigma_0[[k]]))
    log_mod_1 <- log_mod_1 + z2
    
    for(t in 1:nsignaux){
      
      ## 3 
      
      mu_3 <- solution$parametre$Phi_k[[k]]%*%solution$filtre_Kalman$mu_kt_B[[k]][(t),]
      Sigma_3 <- solution$parametre$variances_W_k[[k]]
      
      z1 <- log_densite_normale(solution$filtre_Kalman$mu_kt_B[[k]][(t+1),],mu_3, Sigma_3)
      z2 <- -0.5*solution$parametre_var$lambda_k[k]*(trace(solve(Sigma_3)) + trace(solve(Sigma_3)%*% solution$parametre$Phi_k[[k]]%*% t(solution$parametre$Phi_k[[k]])))
      z3 <- z1+z2
      log_mod_3 <- log_mod_3 + z3
      
      ## 2 
      
      for(i in 1:n){
        
        mu_2 <- solution$filtre_Kalman$mu_kt_B[[k]][t+1,] + facteur_exo[[t]]%*%solution$coefficient[[k]]
        Sigma_2 <- matrix(solution$parametre$variances_V_k[[k]], nrow=d, ncol=d)
        
        z1 <- log_densite_normale(X[which(REFERENCE$TEMPS==t),i], mu_2, Sigma_2)
        z2 <- log(solution$parametre$pi_k[k])
        z3 <- -0.5*solution$parametre_var$lambda_k[k]*trace(solve(Sigma_2))
        z4 <- solution$parametre_var$t_ik[i,k]*(z1+z3)
        
        log_mod_2 <- log_mod_2+z4
        
      }
      
      
    }
    
    prod_tau_ik_pi_k <- 0 
    entropie_1 <- 0
    for(i in 1:n){
      ent_1 <- (-1)*solution$parametre_var$t_ik[i,k]*log(solution$parametre_var$t_ik[i,k])
      entropie_1 <- entropie_1 + ent_1
      prod_tau_ik_pi_k <- prod_tau_ik_pi_k + log(solution$parametre$pi_k[k])*solution$parametre_var$t_ik[i,k]
      
    }
    
    entropie_2 <- ((nsignaux+1)*d/2)*(log(2*pi*exp(1)) + log(solution$parametre_var$lambda_k[k]))
    
    log_mod_4 <- log_mod_4 + entropie_2 + entropie_1
    
    log_mod_2_b <- log_mod_2 + prod_tau_ik_pi_k
  }
  
  ELBO <- log_mod_1 + log_mod_2_b + log_mod_3 + log_mod_4
  
  
  return(ELBO)
}


calcul_BIC <- function(solution, ELBO){
  d <- solution$size$d
  K <- solution$size$nbClasse 
  nsignaux <- solution$size$nsignaux
  n <- solution$size$n
  
  q <- length(solution$coefficient[[1]])
  
  nbr_parametre <- (5*K*d + K-1 + K*q)
  
  
  L <- ELBO
  
  BIC <- -2*L + nbr_parametre*log(n*nsignaux)
  
  return(BIC)
}

normalisation_poids <- function(matrice_poids){
  Z = as.matrix(rowSums(matrice_poids), ncol=1)
  S <- Z
  S[S==0]<-1  
  return(matrice_poids/matrix(rep(S, dim(matrice_poids)[2]), ncol=dim(matrice_poids)[2]))
}

resultat_algorithme_VEM <- function(simulation_VEM, solution, chemin_sauvegarde = "./"){
  
  
  graphique_matrice_classes(solution, simulation_VEM, chemin_sauvegarde)
  
  graphique_proccesus(solution, simulation_VEM,chemin_sauvegarde)
  
  graphique_effet_facteur_AK(solution, simulation_VEM, chemin_sauvegarde)
  
}



mise_a_jours_poids_2 <- function(solution, facteur_exo, X, REFERENCE){
  d <- solution$size$d
  K <- solution$size$nbClasse 
  nsignaux <- solution$size$nsignaux
  n <- solution$size$n
  
  ## t_ik
  numerateur <- matrix(0,nrow=n,ncol=K) 
  
  for(i in 1:n){ 
    for(k in 1:K){
      for(t in 1:nsignaux){
        mu <- solution$filtre_Kalman$mu_kt_B[[k]][(t+1),] + facteur_exo[[t]]%*%solution$coefficient[[k]]
        Sigma <- matrix(solution$parametre$variances_V_k[[k]], nrow=d, ncol=d)
        
        z <- matrix(X[which(REFERENCE$TEMPS == t),i] - mu, nrow=1) %*% solve(Sigma) %*% matrix(X[which(REFERENCE$TEMPS == t),i] - mu, nrow=d)  
        
        
        z_lambda <- solution$parametre_var$lambda_k[k]*trace(solve(Sigma))
        
        a <- -0.5*z 
        b <- log(solution$parametre$pi_k[k]) 
        c <- -(d/2)*log(2*pi) 
        e <- -0.5*log(max(det(Sigma), exp(-323)))
        f <- -0.5*z_lambda
        
        
        
        numerateur[i,k] <- numerateur[i,k] + (a+b+c+e+f)
      }
      
    }
  }
  
  M <- exp(numerateur-matrix(rep(rowMax(numerateur),K),nrow=n, ncol=K))
  MM <- apply(M, c(1, 2), function(x) max(x,exp(-50)))
  
  solution$parametre_var$t_ik <- normalisation_poids(MM)
  return(solution)
}



retraitement_biais <- function(solution, simulation){
  
  
  K <- solution$size$nbClasse
  
  
  if(length(simulation$data_processus)==0){
    K_true <- K
    vec_classes_estimees <- apply(solution$parametre_var$t_ik, 1, which.max)
    table_corresp_classes <- data.frame(estimation = seq(1,K), simulation = seq(1,K))
  }else{
    K_true <- length(simulation$data_processus)
    
    vec_classes_estimees <- apply(solution$parametre_var$t_ik, 1, which.max)
    vec_classes_simulees <- simulation$classes
    
    table_corresp_classes <- create_table_corresp_classes(solution, simulation)
    
  }
  
  d <- solution$size$d
  nsignaux <- solution$size$nsignaux
  
  b_simulation <- simulation$data_processus
  b_estimation <- solution$filtre_Kalman$mu_kt_B
  
  table_processus <- data.frame(simulation=NA, estimation=NA, classe=NA, temps=NA, dim=NA)
  
  
  if(K > K_true){
    for(j in 1:d){
      for(k in 1:K){
        
        indice <- as.numeric(as.character(table_corresp_classes$simulation[table_corresp_classes$estimation==k]))
        table_processus <- rbind(table_processus, data.frame(simulation = as.matrix(b_simulation[[indice]], nrow=nsignaux)[-1,j], 
                                                             estimation = b_estimation[[k]][-1,j] , 
                                                             classe=k, temps = seq(1,nsignaux), dim=j ))   
      }
      
    }  
  }else{
    
    for(j in 1:d){
      for(k in 1:K){
        indice <- k 
        table_processus <- rbind(table_processus, data.frame(simulation = as.matrix(b_simulation[[k]], nrow=nsignaux)[-1,j], 
                                                             estimation = b_estimation[[indice]][-1,j] , 
                                                             classe=k, temps = seq(1,nsignaux), dim=j ))   
      }
      
    } 
  }
  
  
  table_difference <- table_processus[-1,]%>%
    group_by(classe)%>%
    summarise(erreur_moyenne = mean(simulation-estimation)) %>%
    as.data.frame %>%
    rename("classe_process" = 'classe') %>%
    arrange(erreur_moyenne)
  
  termes_biais <- data.frame(classe_coeff = seq(1,K), biais = lapply(solution$coefficient, function(x) x[1]) %>%
                               unlist()) %>%
    arrange(biais)
  
  table_corresp_BIAIS <- cbind(table_difference, termes_biais)
  
  solution$centre_not_retraited <- solution$filtre_Kalman$mu_kt_B
  
  for(k in 1:K){
    indice <- table_corresp_BIAIS$classe_coeff[table_corresp_BIAIS$classe_process ==k]
    solution$filtre_Kalman$mu_kt_B[[k]] <- solution$centre_not_retraited[[k]] + solution$coefficient[[indice]][1]*matrix(1,nrow=(nsignaux +1), ncol=1)
    
  }
  
  
  
  solution$coefficient_BIAIS <- solution$coefficient
  
  
  for(k in 1:K){
    solution$coefficient[[k]][1] <- 0
  }
  
  return(list(solution_new=solution, simulation_new= simulation))
  
}



log_vrais_indiv <- function(i, solution, simulation){
  
  nsignaux <- solution$size$nsignaux
  K <- solution$size$nbClasse
  x_it <- simulation$data_X[,i]
  t_ik <- solution$parametre_var$t_ik[i,]
  
  b_kt <- solution$filtre_Kalman$mu_kt_B
  
  effet_ext <- rep(0,nsignaux)
  for(t in 1:nsignaux){
    for(k in 1:K){
      effet_ext[t] <-  effet_ext[t] +  as.numeric(simulation$data_facteurs[[t]]%*%( t_ik[k]*solution$coefficient[[k]]))
    }
  }
  
  ll <- rep(0, nsignaux)
  
  for(k in 1:K){
    
    ll <- ll +  t_ik[k]*(as.numeric((1/solution$parametre$variances_V_k[[k]]))*(-(x_it - b_kt[[k]][-1,]- effet_ext)^2) - log(as.numeric(solution$parametre$variances_V_k[[k]])))
  }
  
  return(ll)
  
}

compute_log_vrais_VEM <- function(solution, simulation){
  
  table_vrais <- mapply(i= seq(1, solution$size$n), log_vrais_indiv, 
                        MoreArgs= list(solution , simulation))
  
  mean_vraisemblance <- (1/(solution$size$n*solution$size$nsignaux))*sum(table_vrais)
  
  return(mean_vraisemblance)
}


create_table_corresp_classes <- function(solution, simulation){
  
  
  K <- solution$size$nbClasse
  vec_classes_estimees <- apply(solution$parametre_var$t_ik, 1, which.max)
  vec_classes_simulees <- simulation$classes
  
  table_corresp_classes <- table(vec_classes_estimees, vec_classes_simulees) %>%
    as.data.frame %>%
    rename(simulation = vec_classes_simulees, 
           estimation = vec_classes_estimees)
  
  K_true <- length(unique(simulation$classes))
  
  
  if(K < K_true){
    
    
    table_corresp_classes <- table(vec_classes_estimees, vec_classes_simulees) %>%
      as.data.frame %>%
      rename(simulation = vec_classes_simulees, 
             estimation = vec_classes_estimees)
    
    if(length(unique(table_corresp_classes$estimation))== 1){
      
      
      table_corresp_classes <- table(vec_classes_estimees, vec_classes_simulees) %>%
        as.data.frame %>%
        rename(simulation = vec_classes_simulees, 
               estimation = vec_classes_estimees) %>%
        dplyr::select(-Freq)
      
      combi_correspondante <- table_corresp_classes
    }else{
      
      
      combi_possible = permutations(n=max(vec_classes_estimees), r=K_true, v=1:max(vec_classes_estimees), repeats.allowed=T)
      
      cc <- combi_possible %>%
        as.data.frame 
      
      cc$nb <- apply(cc, 1, function(x) length(unique(x)))
      
      combi_possible <- combi_possible[which(cc$nb==max(vec_classes_estimees)),]
      
      base <- seq(1,K_true)
      
      table_trace <- numeric(dim(combi_possible)[1])
      
      for(i in 1:dim(combi_possible)[1] ){
        
        combi <- combi_possible[i,]
        
        corresp <- cbind(base, combi) %>%
          as.data.frame
        
        
        t <- merge(table_corresp_classes, corresp,
                   by.x=c("simulation", "estimation"),
                   by.y=c("base", "combi"))
        
        table_trace[i] <- sum(t$Freq)
        
      }
      
      combi_correspondante <-cbind(combi_possible[which.max(table_trace),], base) %>%
        as.data.frame() %>%
        rename("simulation" = "base", 
               "estimation" = "V1")
    }
  }else if(K == K_true){
    
    
    
    combi_possible = permutations(n=K_true, r=K_true, v=1:K)
    base <- seq(1,K_true)
    
    table_trace <- numeric(dim(combi_possible)[1])
    
    for(i in 1:dim(combi_possible)[1] ){
      
      combi <- combi_possible[i,]
      
      corresp <- cbind(base, combi) %>%
        as.data.frame
      
      
      t <- merge(table_corresp_classes, corresp,
                 by.x=c("simulation", "estimation"),
                 by.y=c("base", "combi"))
      
      table_trace[i] <- sum(t$Freq)
      
    }
    
    
    combi_correspondante <-cbind(combi_possible[which.max(table_trace),], base) %>%
      as.data.frame() %>%
      rename("simulation" = "base", 
             "estimation" = "V1")
    
  }else if (K > K_true){
    
    
    
    combi_possible = permutations(n=K_true, r=K, repeats.allowed = T)
    
    cc <- combi_possible %>%
      as.data.frame 
    
    cc$nb <- apply(cc, 1, function(x) length(unique(x)))
    
    combi_possible <- combi_possible[which(cc$nb==max(vec_classes_simulees)),]
    
    base <- seq(1,K)
    
    table_trace <- numeric(dim(combi_possible)[1])
    
    for(i in 1:dim(combi_possible)[1] ){
      
      combi <- combi_possible[i,]
      
      corresp <- cbind(base, combi) %>%
        as.data.frame
      
      
      t <- merge(table_corresp_classes, corresp,
                 by.x=c("simulation", "estimation"),
                 by.y=c("base", "combi"))
      
      table_trace[i] <- sum(t$Freq)
    }
    
    
    combi_correspondante <-cbind(combi_possible[which.max(table_trace),], base) %>%
      as.data.frame() %>%
      rename("simulation" = "V1", 
             "estimation" = "base")
  }
  
  
  
  return(combi_correspondante)
}

initialisation_VEM_Kmeans_AK <- function(X, K, facteur_exo, REFERENCE, nsignaux){
  
  deb_init <- tic()
  
  solution <- list()
  n <- dim(X)[2]
  d <- 1
  nb_var <- dim(facteur_exo[[1]])[2]/d
  nbClasse <- K
  solution$size <- list()
  solution$size$d <- d
  solution$size$n <- n
  solution$size$nsignaux <- nsignaux
  solution$size$nbClasse <- K
  
  
  
  ## INITIALISATION DES PARAMETRE
  
  variances_Vk<- list() 
  variances_Wk <- list() 
  
  for(k in 1:K){
    variances_Wk[[k]] <- diag(1,d)
    variances_Vk[[k]] <- diag(1,d)
  }
  
  solution$parametre$variances_W_k <- variances_Wk 
  solution$parametre$variances_V_k <- variances_Vk 
  
  pi_k <- rep(1/K, K)
  
  ## Filtre de Kalman 
  mu_kt_F <- list()
  mu_kt_B <- list()
  
  P_kt_F <- list()
  P_kt_B <- list()
  
  for(k in 1:K){
    P_kt_F[[k]] <- list()
    P_kt_B[[k]] <- list()
  }
  
  for(t in 1:(nsignaux+1)){
    for(k in 1:K){
      P_kt_F[[k]][[t]] <- diag(1, d)
      P_kt_B[[k]][[t]] <- diag(1, d)
    }
  }
  
  
  data <- X 
  
  ## MISE A JOUR DES PROBABILITES D'APPARTENANCE : 
  
  data_k_means <- t(X)
  
  #2- Classification k-means
  
  res_k_means <- kmeans(data_k_means, K, iter.max=20, nstart=15)
  classes_initiales <- res_k_means$cluster
  
  
  solution$parametre$pi_k <- rep(1/K, K)
  
  #3- Bruitage et normalisation  
  
  t_ik <- matrix(nrow=n, ncol= K) %>%
    as.data.frame %>%
    mutate(classes= classes_initiales) 
  for(k in 1:K){
    t_ik[which(t_ik$classes == k),k] <- 1
    t_ik[which(t_ik$classes != k),k] <- 0
  }
  set.seed(45)
  t_ik <-  t_ik %>%
    dplyr::select(-classes)%>%
    apply( c(1,2), function(x) x +  runif(1, min=0, max=0.3)) %>%
    normalisation_poids()
  set.seed(NULL)
  
  mu_k <- list() 
  
  for(k in 1:K){
    
    mu_k[[k]] <- matrix(0, nrow=nsignaux, ncol= d)
  }
  
  iter= 0 
  
  coefficient_A <- list()  
  
  
  
  while(iter<2){
    
    for(k in 1:K){
      INDICE <- merge(data.frame(individu=seq(1,n)), data.frame(temps=seq(1,nsignaux)))
      
      
      if(d*nb_var==1){
        
        t2 <- sum(mapply(function(i,t)
          t_ik[i,k]* t(facteur_exo[[t]]) %*% solve(variances_Vk[[k]])%*%(X[which(REFERENCE$TEMPS==t),i] - mu_k[[k]][t,]), 
          INDICE$individu, INDICE$temps, SIMPLIFY = "array"))
        
        t1 <- sum(mapply(function(i,t) t_ik[i,k]* t(facteur_exo[[t]]) %*% solve(variances_Vk[[k]]) %*% facteur_exo[[t]], 
                         INDICE$individu, INDICE$temps, SIMPLIFY = "array"))
        
        
        
      }else{
        t2 <- rowSums(mapply(function(i,t)
          t_ik[i,k]* t(facteur_exo[[t]]) %*% solve(variances_Vk[[k]])%*%(X[which(REFERENCE$TEMPS==t),i] - mu_k[[k]][t,]), 
          INDICE$individu, INDICE$temps, SIMPLIFY = "array"), dims=2)
        
        t1 <- rowSums(mapply(function(i,t) t_ik[i,k]* t(facteur_exo[[t]]) %*% solve(variances_Vk[[k]]) %*% facteur_exo[[t]], 
                             INDICE$individu, INDICE$temps,  SIMPLIFY = "array"), dims=2)
        
        
      }
      
      coefficient_A[[k]] <- solve(t1)%*%t2
    }
    
    #### Initialisation des centres de classes 
    
    
    X_bis <- matrix(nrow = dim(X)[1], ncol= dim(X)[2])
    for(i in 1:n){
      for(t in 1:nsignaux){
        if(d==1){
          X_bis[which(REFERENCE$TEMPS==t), i] <-  sum(sapply(seq(1,K), function(k, i, t) 
            t_ik[i,k]*(X[which(REFERENCE$TEMPS==t), i] - facteur_exo[[t]]%*%coefficient_A[[k]]),i,t))
          
        }else{
          X_bis[which(REFERENCE$TEMPS==t), i] <-  rowSums(sapply(seq(1,K), function(k, i, t) 
            t_ik[i,k]*(X[which(REFERENCE$TEMPS==t), i] - facteur_exo[[t]]%*%coefficient_A[[k]]),i,t))
          
        }
      }
    }
    
    
    data_k_means <- t(X_bis)
    #2- Classification k-means
    
    res_k_means <- kmeans(t(X_bis), K, iter.max=20, nstart=15)
    classes_initiales <- res_k_means$cluster
    
    #3- Bruitage et normalisation  
    
    t_ik <- matrix(nrow=n, ncol= K) %>%
      as.data.frame %>%
      mutate(classes= classes_initiales) 
    
    for(k in 1:K){
      t_ik[which(t_ik$classes == k),k] <- 1
      t_ik[which(t_ik$classes != k),k] <- 0
    }
    
    set.seed(45)
    t_ik <-  t_ik %>%
      dplyr::select(-classes)%>%
      apply( c(1,2), function(x) x +  runif(1, min=0, max=0.3)) %>%
      normalisation_poids()
    set.seed(NULL)
    
    classes_initiales <- res_k_means$cluster
    
    for(k in 1:K){
      indice <- which(classes_initiales==k)
      XX <- X_bis[,indice]
      
      if(d != 1){
        mu_k[[k]] <- cbind(XX, REFERENCE) %>%
          gather(key='indiv', value="value", -DIM, -TEMPS)%>%
          group_by(DIM, TEMPS) %>%
          summarise(mm = mean(value)) %>%
          pivot_wider(names_from = TEMPS, values_from = mm) %>%
          as.data.frame%>%
          column_to_rownames(var="DIM")%>%
          t
      }else{ 
        mu_k[[k]] <- cbind(XX, REFERENCE) %>%
          gather(key='indiv', value="value", -DIM, -TEMPS)%>%
          group_by(TEMPS) %>%
          summarise(mm = mean(value)) %>%
          pivot_wider(names_from = TEMPS, values_from = mm) %>%
          as.data.frame%>%
          t
        
      }
      
    }
    
    iter <- iter + 1
    
    
  }
  
  
  solution$coefficient <-coefficient_A
  
  solution$parametre_var$lambda_k <- rep(1,K)
  
  solution$parametre$mu_0 <- list()
  solution$parametre$Sigma_0 <- list()
  
  for(k in 1:K){
    if(d!=1){
      solution$parametre$mu_0[[k]] <- colMeans(mu_k[[k]][1:min(7, nsignaux),])
      solution$parametre$Sigma_0[[k]] <- diag(2,d)
    }else{
      solution$parametre$mu_0[[k]] <- mean(mu_k[[k]][1:min(7, nsignaux),])
      solution$parametre$Sigma_0[[k]] <- diag(2,d)
    }
    
  }
  
  
  for(k in 1:K){
    if(d!=1){
      mu_k[[k]] <- rbind(colMeans(mu_k[[k]][1:min(10, nsignaux),]), mu_k[[k]])
    }else{
      mu_k[[k]] <- rbind(mean(mu_k[[k]][1:min(10, nsignaux),]), mu_k[[k]])
    }
  }
  
  
  solution$parametre_var$t_ik <- t_ik
  
  solution$filtre_Kalman$mu_kt_F <- mu_k # Forward 
  solution$filtre_Kalman$mu_kt_B <- mu_k # Backward
  
  
  solution$filtre_Kalman$P_kt_F <- P_kt_F # Forward
  solution$filtre_Kalman$P_kt_B <- P_kt_B # Backward
  
  
  ## Initalisation Phi_k 
  
  solution$parametre$Phi_k <- list()
  
  for(k in 1:K){
    if(d!=1){
      M1 <- rowSums(sapply(seq(2, (nsignaux+1)), function(t,k) mu_k[[k]][t,]%*%t(mu_k[[k]][(t-1),]),k, simplify="array" ), dims=2)
      
      M2 <- rowSums(sapply(seq(2, (nsignaux+1)), function(t,k,d) mu_k[[k]][t,]%*%t(mu_k[[k]][(t),])+ diag(1,d) ,k,d , simplify="array" ), dims=2)
      
    }else{ 
      M1 <- sum(sapply(seq(2, (nsignaux+1)), function(t,k) mu_k[[k]][t,]%*%t(mu_k[[k]][(t-1),]),k, simplify="array" ))
      
      M2 <- sum(sapply(seq(2, (nsignaux+1)), function(t,k,d) mu_k[[k]][t,]%*%t(mu_k[[k]][(t),])+ diag(1,d) ,k,d , simplify="array" ))
      
    }
    
    solution$parametre$Phi_k[[k]] <- diag(diag(M1%*%solve(M2)),d)
    
  }
  
  ## Initialisation des poids t_ik à 0 pour pouvoir remplir ensuite 
  
  solution$parametre_var$t_ik <- t_ik
  
  fin_init <- toc(quiet=F)
  
  
  
  return(solution)
}


iterations_VEM_AK <- function(simulation, K, TT, solution){
  
  table_res <- data.frame(erreur_centres =NA, erreur_facteur = NA)
  X <- simulation$data_X 
  REFERENCE <- simulation$REFERENCE
  facteur_exo <- simulation$data_facteurs
  
  d <- solution$size$d
  K <- solution$size$nbClasse 
  nsignaux <- solution$size$nsignaux
  n <- solution$size$n
  
  table_erreurs <- data.frame(iter=NA, erreur_centres = NA, erreur_facteur = NA, taux_classif = NA, ELBO = NA)
  
  # ITERATION 
  iter = 0
  log_vraisemblance <- NA
  ml_av <- -1000000
  fin_iter <- F
  debut <- tic()
  
  while(iter < 500 &  !fin_iter){
    dITER <- tic()
    ### MISE A JOUR DES CENTRES DE CLASSES ###
    c_kt <- list()
    m_kt <- list()
    
    for(k in 1:K){
      c_kt[[k]] <- matrix(nrow=(nsignaux+1), ncol=d)
      m_kt[[k]] <- matrix(nrow=(nsignaux+1), ncol=d)
    }
    
    P_kt_F <- solution$filtre_Kalman$P_kt_F
    P_kt_B <- solution$filtre_Kalman$P_kt_B
    
    for(k in 1:K){
      
      ## INIT 
      c_kt[[k]][1,]<- solution$parametre$mu_0[[k]]
      P_kt_F[[k]][[1]] <- solution$parametre$Sigma_0[[k]]
      
      # FORWARD
      for(t in 2:(nsignaux+1)){
        
        if(d==1){
          ss_1 <- sum(sapply(seq(1,n), function(i,t,k) 
            solution$parametre_var$t_ik[i,k]*(X[which(REFERENCE$TEMPS == (t-1)),i] - facteur_exo[[t-1]]%*%solution$coefficient[[k]] - solution$parametre$Phi_k[[k]]%*%matrix(c_kt[[k]][(t-1),],ncol=1)), t,k)) 
          
        }else{
          ss_1 <- rowSums(sapply(seq(1,n), function(i,t,k) 
            solution$parametre_var$t_ik[i,k]*(X[which(REFERENCE$TEMPS == (t-1)),i]- facteur_exo[[t-1]]%*%solution$coefficient[[k]] - solution$parametre$Phi_k[[k]]%*%matrix(c_kt[[k]][(t-1),], 
                                                                                                                                                                            ncol=1)), t,k)) 
        }
        
        
        P_kt_F[[k]][[t]] <- solve( solve(solution$parametre$Phi_k[[k]]%*%P_kt_F[[k]][[t-1]]%*%solution$parametre$Phi_k[[k]]+ solution$parametre$variances_W_k[[k]])+ solve(solution$parametre$variances_V_k[[k]])*(sum(solution$parametre_var$t_ik[,k])))
        
        c_kt[[k]][t,] <- solution$parametre$Phi_k[[k]]%*%c_kt[[k]][(t-1),] + P_kt_F[[k]][[t]]%*%solve(solution$parametre$variances_V_k[[k]])%*%ss_1
        
      } 
    }
    
    
    ## Backward 
    
    # Init 
    J_kt <- list()
    
    for(k in 1:K){
      
      J_kt[[k]]<-list()
      
      J_kt[[k]][[nsignaux+1]] <- matrix(nrow=d, ncol=d) 
      
      m_kt[[k]][(nsignaux+1),] <- c_kt[[k]][(nsignaux+1),]
      P_kt_B[[k]][[(nsignaux+1)]] <- P_kt_F[[k]][[(nsignaux+1)]]
      
      M <- solve(solution$parametre$Phi_k[[k]]%*%P_kt_F[[k]][[(nsignaux+1)]]%*%solution$parametre$Phi_k[[k]] + solution$parametre$variances_W_k[[k]])
      
      J_kt[[k]][[nsignaux+1]] <- P_kt_F[[k]][[(nsignaux+1)]]%*%solution$parametre$Phi_k[[k]]%*%M
      
      
      for(t in (nsignaux):1){
        
        M <- solve(solution$parametre$Phi_k[[k]]%*%P_kt_F[[k]][[(t)]]%*%solution$parametre$Phi_k[[k]] + solution$parametre$variances_W_k[[k]])
        
        
        J_kt[[k]][[t]] <- P_kt_F[[k]][[(t)]]%*%solution$parametre$Phi_k[[k]]%*%M
        
        m_kt[[k]][t,] <- c_kt[[k]][t,] + J_kt[[k]][[t]]%*%(m_kt[[k]][t+1,] - solution$parametre$Phi_k[[k]]%*%c_kt[[k]][t,])
        
        p_kt_moins <- solution$parametre$Phi_k[[k]]%*%P_kt_F[[k]][[t]]%*%t(solution$parametre$Phi_k[[k]]) + solution$parametre$variances_W_k[[k]]
        
        P_kt_B[[k]][[t]] <- P_kt_F[[k]][[t]] - J_kt[[k]][[t]] %*%(p_kt_moins - P_kt_B[[k]][[(t+1)]])%*%t(J_kt[[k]][[t]] )
        
        
        
      }
      
    }
    
    
    
    
    
    # Centrage de l'ensemble des centres de classes
    if(K !=1){
      mb=matrix(0,K,d)
      for (k in 1:K)
      {
        mb[k,] =colMeans(m_kt[[k]])
      }
      
      mmb=colMeans(mb)
      
      for (k in 1:K)
      {
        m_kt[[k]] = m_kt[[k]]- matrix(1,nsignaux+1,1)%*%mb[k,]
      }
      
    }
    
    
    
    solution$filtre_Kalman$P_kt_F <- P_kt_F
    solution$filtre_Kalman$P_kt_B <- P_kt_B
    
    solution$filtre_Kalman$mu_kt_F <- c_kt
    solution$filtre_Kalman$mu_kt_B <- m_kt 
    
    
    ### MISE A JOUR DES POIDS ###
    if(d==1){
      
      numerateur<- matrix(0,nrow=n,ncol=K) 
      
      for(k in 1:K){
        
        effet_ext <- numeric(nsignaux)
        
        for(t in 1:nsignaux){
          effet_ext[t] <- facteur_exo[[t]]%*%solution$coefficient[[k]]
        }
        
        
        a <- nsignaux*(log(solution$parametre$pi_k[k]) -0.5*solution$parametre_var$lambda_k[k]*trace(solve(solution$parametre$variances_V_k[[k]])))
        
        
        ##########
        
        mu <- solution$filtre_Kalman$mu_kt_B[[k]][-1]*matrix(1,nrow=solution$size$nsignaux,ncol=solution$size$n) + effet_ext
        
        
        Sigma <- solution$parametre$variances_V_k[[k]][1,1]
        
        XM <- X-mu
        XM_s <- XM^2
        RES <- log((1/sqrt(2*pi*Sigma))*exp((-1/(2*Sigma))*XM_s)) 
        
        numerateur[,k] <- colSums(RES) + a
        
      }
      
      
      M <-exp(numerateur-matrix(rep(rowMax(numerateur),K),nrow=n, ncol=K))
      MM <- apply(M, c(1, 2), function(x) max(x,exp(-50)))
      tt <- normalisation_poids(MM)
      
      if(sum(is.na(tt))==0){
        
        solution$parametre_var$t_ik <- tt
        
      }
      
      
      
      
    }else{
      ## Pour le moment, le code précédent est adapté au cas d=1 
      solution <- mise_a_jours_poids_2(solution, facteur_exo, X, REFERENCE)
    }
    
    
    ### MISE A JOURS DES LAMBDA_K ###
    
    lambda_k<- numeric(K)
    
    for(k in 1:K){
      
      numerateur <- (nsignaux+1)*d
      a <- nsignaux*trace(solve(solution$parametre$variances_V_k[[k]]))*sum(solution$parametre_var$t_ik[,k])
      b <- nsignaux*trace(solve(solution$parametre$variances_W_k[[k]]))
      c <- nsignaux*trace(solve(solution$parametre$variances_W_k[[k]])%*%solution$parametre$Phi_k[[k]]%*%solution$parametre$Phi_k[[k]])
      e <- trace(solve(solution$parametre$Sigma_0[[k]]))
      denominateur <-  a+b+c+e
      
      lambda_k[k] <- numerateur/denominateur
      
      
    }
    solution$parametre_var$lambda_k <- lambda_k
    
    
    ### MISE A JOUR DES PI ###
    K <- solution$size$nbClasse 
    n <- solution$size$n
    
    solution$parametre$pi_k <- rep(1/K, K)
    
    ### MISE A JOUR DES MU_0 ET SIGMA_0 ###
    
    for(k in 1:K){
      solution$parametre$mu_0[[k]]<- solution$filtre_Kalman$mu_kt_B[[k]][1,]  
      solution$parametre$Sigma_0[[k]] <- diag(solution$parametre_var$lambda_k[k], d)
    }
    
    
    
    ### MISE A JOUR DE V_k ###
    Vk <- list()
    
    
    for(k in 1:K){
      
      effet_ext <- numeric(nsignaux)
      for(t in 1:nsignaux){
        effet_ext[t] <- facteur_exo[[t]]%*%solution$coefficient[[k]]
      }
      
      mu <- solution$filtre_Kalman$mu_kt_B[[k]][-1]*matrix(1,nrow=solution$size$nsignaux,ncol=solution$size$n) + effet_ext
      
      XM <- X - mu
      XM_s <- XM^2
      
      XM_s[t,1]
      solution$parametre_var$t_ik[1,k]*(XM_s[1,1] + solution$parametre_var$lambda_k[k])
      
      
      RES <- solution$parametre_var$t_ik[,k]*t(XM_s + solution$parametre_var$lambda_k[k])
      
      Vk[[k]] <- (1/(nsignaux*sum(solution$parametre_var$t_ik[,k])))*matrix(sum(RES), nrow=d, ncol=d)
      Vk[[k]] <- matrix(1, nrow=d, ncol=d)
    }
    
    solution$parametre$variances_V_k <- Vk
    
    
    
    ### MISE A JOUR DE W_K ###
    
    for(k in 1:K){
      
      mu <- solution$filtre_Kalman$mu_kt_B[[k]][-1,] - solution$parametre$Phi_k[[k]]%*%solution$filtre_Kalman$mu_kt_B[[k]][-(nsignaux+1),]
      
      M <- (1/nsignaux)*(sum(mu^2)) + solution$parametre_var$lambda_k[k]*(1+ solution$parametre$Phi_k[[k]]^2)
      
      solution$parametre$variances_W_k[[k]] <- matrix(M, nrow=d, ncol=d) 
      solution$parametre$variances_W_k[[k]] <- matrix(0.2, nrow=d, ncol=d) 
      
    }
    
    ### MISE A JOUR DES PHI_K ###
    
    
    for(k in 1:K){
      
      somme_1 <- 0
      somme_2 <- 0 
      
      epsilon <- 0.0000165 
      
      for(t in 2:nsignaux+1){
        
        z1 <- solution$filtre_Kalman$mu_kt_B[[k]][t,] %*% t(solution$filtre_Kalman$mu_kt_B[[k]][t-1,])
        z2 <- solution$filtre_Kalman$mu_kt_B[[k]][t-1,] %*% t(solution$filtre_Kalman$mu_kt_B[[k]][t-1,]) + diag(solution$parametre_var$lambda_k[k],d)
        
        somme_1 <- somme_1 + z1
        somme_2 <- somme_2 + z2
        
      }
      
      M <- somme_1%*%solve(somme_2)
      
      phi_k <- diag(diag(M),d)
      
      if(sum(diag(phi_k)>=1)!=0){
        diag(phi_k)[abs(diag(phi_k))>=1] <- 1-epsilon
      }
      
      solution$parametre$Phi_k[[k]]<- phi_k
    }
    
    
    ### MISE A JOUR DES COEFFICIENTS A ###
    X <- simulation$data_X 
    REFERENCE <- simulation$REFERENCE
    facteur_exo <- simulation$data_facteurs
    d <- solution$size$d
    K <- solution$size$nbClasse 
    nsignaux <- solution$size$nsignaux
    n <- solution$size$n
    dim_U <- length(solution$coefficient[[1]])
    
    for(k in 1:K){
      
      X_M <- X-solution$filtre_Kalman$mu_kt_B[[k]][-1]*matrix(1,nrow=solution$size$nsignaux,
                                                              ncol=solution$size$n)
      
      somme_effect_square <- 0
      somme_effect <- rep(0, dim_U)
      for(t in 1:nsignaux){
        somme_effect <- somme_effect + rowSums(t(facteur_exo[[t]])%*%(solution$parametre_var$t_ik[,k]*t(X_M[t,])))
      }
      
      for(t in 1:nsignaux){
        somme_effect_square <- somme_effect_square + t(facteur_exo[[t]])%*%facteur_exo[[t]]
      }
      
      
      argument_1 <- matrix(0, nrow=dim_U, ncol=dim_U)
      argument_2 <- rep(0,dim_U)
      
      argument_1 <- sum(solution$parametre_var$t_ik[,k]*(1/solution$parametre$variances_V_k[[k]][1,1]))*somme_effect_square
      argument_2 <- (1/solution$parametre$variances_V_k[[k]][1,1])%*%somme_effect
      
      solution$coefficient[[k]] <-  solve(argument_1)%*%t(argument_2)
      
    }
    
    if(iter==0){
      classe_av = solution$filtre_Kalman$mu_kt_B
    }else{
      classe_new = solution$filtre_Kalman$mu_kt_B
      
      ss = 0
      for(k in 1:solution$size$nbClasse){
        
        ss= ss+ (mean((classe_av[[k]]- classe_new[[k]])^2))
        
      }
      if(ss <0.00001){
        fin_iter <- TRUE
      }
    }
    
    
    classe_av <-solution$filtre_Kalman$mu_kt_B   
    
    iter <- iter + 1
    
    table_res<- rbind(table_res, c(erreur_moyenne_centres(simulation = simulation ,
                                                          solution = solution), 
                                   erreur_effet_facteurs(simulation = simulation,
                                                         solution = solution)))
    
    
  }
  
  
  elbo <- compute_ELBO(solution, facteur_exo = facteur_exo,X = X, REFERENCE = REFERENCE)
  bb <- calcul_BIC(solution, elbo)
  
  fin <- toc(quiet = TRUE)
  
  
  temps_calcul <- fin$toc-debut
  
  
  return(list(BIC=bb, best_solution = solution, nbr_iteration = iter,
              temps_calcul= temps_calcul, elbo=elbo))
  
  
}

algorithme_VEM_AK <- function(simulation, K, TT){
  
  
  solution <- initialisation_VEM_Kmeans_AK(X= simulation$data_X,
                                           K=K, facteur_exo=simulation$data_facteurs,
                                           REFERENCE=simulation$REFERENCE, nsignaux = TT)
  
  
  res <- iterations_VEM_AK(simulation = simulation, K, TT, solution =solution)
  
  
  return(res)
}
