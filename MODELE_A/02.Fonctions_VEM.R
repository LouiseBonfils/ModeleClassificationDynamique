initialisation_VEM_Kmeans <- function(X, K, facteur_exo, REFERENCE, nsignaux){
  

  solution <- list()
  n <- dim(X)[2]
  d <- length(unique(REFERENCE$DIM))
  if(d == 1){
    nb_var <- length(facteur_exo[[1]])
    
  }else{
    nb_var <- dim(facteur_exo[[1]])[2]/d
  }
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
  
  if(K == 1){
    classes_initiales <- rep(1, n)
  }else if(K == nrow(data_k_means) & K>1){
    res_k_means <- kmeans(data_k_means, K-1, iter.max=40, nstart=5)
    classes_initiales <- seq(1, K)
  }else{
    res_k_means <- kmeans(data_k_means, K, iter.max=40, nstart=5)
    classes_initiales <- res_k_means$cluster
    
  }
  
  
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
    apply( c(1,2), function(x) x +  runif(1, min=0, max=0.05)) %>%
    normalisation_poids()
  set.seed(NULL)
  
  mu_k <- list() 
  
  for(k in 1:K){
    
    mu_k[[k]] <- matrix(0, nrow=nsignaux, ncol= d)
  }
  
  iter= 0 
  
  while(iter<2){
    
    
    somme_effect_square <- 0
    somme_effect <- list()
    if(d==1){
      
      for(k in 1:K){
        
        X_M <- X-mu_k[[k]][,1]*matrix(1,nrow=nsignaux, ncol=n)
        
        somme_effect[[k]] <- rep(0, nb_var)
        for(t in 1:nsignaux){

          somme_effect[[k]] <- somme_effect[[k]] + rowSums(t(facteur_exo[[t]])%*%(t_ik[,k]*t(X_M[t,])))
        }
      }
      for(t in 1:nsignaux){
        somme_effect_square <- somme_effect_square + t(facteur_exo[[t]])%*%facteur_exo[[t]]
      }
      
      
      argument_1 <- matrix(0, nrow=nb_var, ncol=nb_var)
      argument_2 <- rep(0,nb_var)
      for(k in 1:K){
        argument_1 <- argument_1 + sum(t_ik[,k]*variances_Vk[[k]][1,1])*somme_effect_square
        argument_2 <- argument_2 + variances_Vk[[k]]%*%somme_effect[[k]]
      }
      
      
      coefficient_A <- solve(argument_1)%*%t(argument_2)
      
    }else{
      INDICE <- merge(data.frame(classe=seq(1,K)), 
                      merge(data.frame(individu=seq(1,n)), data.frame(temps=seq(1,nsignaux))))
      
      s2 <- 0
      s1 <- 0
      for(k in 1:K){
        for(t in 1:nsignaux){
          a2 <- t(facteur_exo[[t]]) %*% solve(variances_Vk[[k]]) 
          a1 <- t(facteur_exo[[t]]) %*% solve(variances_Vk[[k]]) %*% facteur_exo[[t]]  
          
          for(i in 1:n){
            b2 <- (X[which(REFERENCE$TEMPS==t),i] - mu_k[[k]][t,])
            s2 <- s2 +  t_ik[i,k]*(a2%*%b2)
            s1 <- s1 +  t_ik[i,k]*a1
          } 
        } 
        #}
        
      }
      
      coefficient_A <- solve(s1)%*%s2
      
    }
    
    
    
    
    #### Initialisation des centres de classes 
    
    
    X_bis <- matrix(nrow = dim(X)[1], ncol= dim(X)[2])
    
    for(i in 1:n){
      for(t in 1:nsignaux){
    
        X_bis[which(REFERENCE$TEMPS==t), i] <-  X[which(REFERENCE$TEMPS==t), i] - facteur_exo[[t]]%*%coefficient_A
      }
    }
    
    #2- Classification k-means
    if(K == 1){
      classes_initiales <- rep(1, n)
    }else if(K == nrow(data_k_means) & K>1){
      res_k_means <- kmeans(t(X_bis), K-1, iter.max=10, nstart=15)
      classes_initiales <- seq(1, K)
      res_k_means$cluster <- seq(1, K)
      
    }else{
      res_k_means <- kmeans(t(X_bis), K, iter.max=10, nstart=15)
      classes_initiales <- res_k_means$cluster
      
    }
    
    
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
      apply( c(1,2), function(x) x +  runif(1, min=0, max=0.05)) %>%
      normalisation_poids()
    set.seed(NULL)
    
    
    for(k in 1:(K)){
      
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
 
    # Les centres de classes sont centrés
    
    if(K!=1){
      mb=matrix(0,K,d)
      for (k in 1:K){
        mb[k,] =colMeans(mu_k[[k]])
      }
      
      mmb=colMeans(mb)
      
      
      for (k in 1:K){
        mu_k[[k]] = mu_k[[k]]-matrix(1,nsignaux,1)%*%mmb
      }
    }
  
    iter <- iter + 1
  }
  
  
  solution$coefficient <-coefficient_A
  
  solution$parametre_var$lambda_k <- rep(1, K)
  
  solution$parametre$mu_0 <- list()
  solution$parametre$Sigma_0 <- list()
  
  for(k in 1:K){
    if(d!=1){
      solution$parametre$mu_0[[k]] <- colMeans(mu_k[[k]][1:min(nsignaux, 7),])
      solution$parametre$Sigma_0[[k]] <- diag(2,d)
    }else{
      solution$parametre$mu_0[[k]] <- mean(mu_k[[k]][1:min(nsignaux, 7),])
      solution$parametre$Sigma_0[[k]] <- diag(2,d)
    }
    
  }
  
  
  for(k in 1:K){
    if(d!=1){
      mu_k[[k]] <- rbind(colMeans(mu_k[[k]][1:min(nsignaux, 7),]), mu_k[[k]])
    }else{
      mu_k[[k]] <- rbind(mean(mu_k[[k]][1:min(nsignaux, 7),]), mu_k[[k]])
    }
  }
  
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
  
  return(solution)
}


log_densite_normale <- function(x, mu, Sigma){

  d <- length(x)
  res <- -(d/2)*log(2*pi) - (1/2)*log(det(Sigma)) -
    (1/2)*(matrix((x-mu), nrow=1, ncol=d)%*%solve(Sigma)%*%matrix((x-mu), 
                                                                  nrow=d, ncol=1))
  
  return(res)
}



densite_normale <- function(x, mu, sigma){
  
  a <- 1/sqrt(2*pi*sigma)
  b <- (-1/2)*( ((x-mu)^2)/sigma )
  
  res <- a*exp(b)
  
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

  X <- simulation$data_X 
  REFERENCE <- simulation$REFERENCE
  facteur_exo <- simulation$data_facteurs
  
  d <- solution$size$d
  K <- solution$size$nbClasse 
  nsignaux <- solution$size$nsignaux
  n <- solution$size$n
  
  table_erreurs <- data.frame(iter=NA, erreur_centres = NA, erreur_facteur = NA, erreur_totale = NA,  elbo = NA)
  
  # ITERATION 
  iter = 0
  log_vraisemblance <- NA
  ml_av <- -1000000
  fin_iter <- F
  debut <- tic()
  
  while(iter < 750 && !fin_iter){
    
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
            solution$parametre_var$t_ik[i,k]*(X[which(REFERENCE$TEMPS == (t-1)),i] - facteur_exo[[t-1]]%*%solution$coefficient - solution$parametre$Phi_k[[k]]%*%matrix(c_kt[[k]][(t-1),],
                                                                                                                                                                        ncol=1)), t,k)) 
          
        }else{
          ss_1 <- rowSums(sapply(seq(1,n), function(i,t,k) 
            solution$parametre_var$t_ik[i,k]*(X[which(REFERENCE$TEMPS == (t-1)),i]- facteur_exo[[t-1]]%*%solution$coefficient- solution$parametre$Phi_k[[k]]%*%matrix(c_kt[[k]][(t-1),], 
                                                                                                                                                                      ncol=1)), t,k)) 
        }
        
        
        P_kt_F[[k]][[t]] <- solve( solve(solution$parametre$Phi_k[[k]]%*%P_kt_F[[k]][[t-1]]%*%t(solution$parametre$Phi_k[[k]]) + solution$parametre$variances_W_k[[k]])+ solve(solution$parametre$variances_V_k[[k]])*(sum(solution$parametre_var$t_ik[,k])))
        
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
      
      P_kt_B[[k]][[(nsignaux+1)]]
      
      for(t in (nsignaux):1){
        
        M <- solve(solution$parametre$Phi_k[[k]]%*%P_kt_F[[k]][[(t)]]%*%t(solution$parametre$Phi_k[[k]]) + solution$parametre$variances_W_k[[k]])
        J_kt[[k]][[t]] <- P_kt_F[[k]][[(t)]]%*%solution$parametre$Phi_k[[k]]%*%M
        
        m_kt[[k]][t,] <- c_kt[[k]][t,] + J_kt[[k]][[t]]%*%(m_kt[[k]][t+1,]-c_kt[[k]][t,])
       
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
        m_kt[[k]] = m_kt[[k]]- matrix(1,nsignaux+1,1)%*%mmb
      }

    }

    solution$filtre_Kalman$P_kt_F <- P_kt_F
    solution$filtre_Kalman$P_kt_B <- P_kt_B
    
    solution$filtre_Kalman$mu_kt_F <- c_kt
    solution$filtre_Kalman$mu_kt_B <- m_kt 
  
    
    ### MISE A JOUR DES POIDS ###
    if(d==1){
      
      numerateur <- matrix(0,nrow=n,ncol=K) 
      
      effet_ext <- numeric(nsignaux)
      for(t in 1:nsignaux){
        effet_ext[t] <- facteur_exo[[t]]%*%solution$coefficient
      }
      
      for(k in 1:K){
        
        ab <- nsignaux*(log(solution$parametre$pi_k[k]) -0.5*solution$parametre_var$lambda_k[k]*trace(solve(solution$parametre$variances_V_k[[k]])))
        
        a <- log(solution$parametre$pi_k[k]) -nsignaux*(0.5*solution$parametre_var$lambda_k[k]*trace(solve(solution$parametre$variances_V_k[[k]])))
        
        ##########
        
        mu <- solution$filtre_Kalman$mu_kt_B[[k]][-1]*matrix(1,nrow=solution$size$nsignaux,ncol=solution$size$n) + effet_ext
        
        Sigma <- solution$parametre$variances_V_k[[k]][1,1]
        
        XM <- X-mu
        XM_s <- XM^2
        RES <- log((1/sqrt(2*pi*Sigma))*exp((-1/(2*Sigma))*XM_s)) 
        
        
        
        numerateur[,k] <- colSums(RES) + a

        
      }
      
      
      M <- exp(numerateur-matrix(rep(rowMax(numerateur),K),nrow=n, ncol=K))
      MM <- apply(M, c(1, 2), function(x) max(x,exp(-50)))
      solution$parametre_var$t_ik <- normalisation_poids(MM)
      
    
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
      
      solution$parametre$mu_0[[k]] <- solution$filtre_Kalman$mu_kt_B[[k]][1,]  
      
      
      solution$parametre$Sigma_0[[k]] <- diag(solution$parametre_var$lambda_k[k], d)
    }

    ### MISE A JOUR DE V_k ###

    Vk <- list()
    if(d==1){
      
      
      effet_ext <-matrix(nrow = nsignaux, ncol=d)
      for(t in 1:nsignaux){
        effet_ext[t, ] <- facteur_exo[[t]]%*%solution$coefficient
      }
      
    
      for(k in 1:K){
        
        mu <- solution$filtre_Kalman$mu_kt_B[[k]][-1]*matrix(1,nrow=solution$size$nsignaux,ncol=solution$size$n) + matrix(effet_ext,nrow=solution$size$nsignaux,ncol=solution$size$n)
        
        
        XM <- X - mu
        XM_s <- XM^2
        
        RES <- solution$parametre_var$t_ik[,k]*t(XM_s + solution$parametre_var$lambda_k[k])
        
        Vk[[k]] <- (1/(nsignaux*sum(solution$parametre_var$t_ik[,k])))*matrix(sum(RES), nrow=d, ncol=d)
        
      }
      
      solution$parametre$variances_V_k <- Vk
      
      
      
      ### MISE A JOUR DE W_K ###
      for(k in 1:K){
        
        mu <- solution$filtre_Kalman$mu_kt_B[[k]][-1,] - 
          solution$parametre$Phi_k[[k]]%*%solution$filtre_Kalman$mu_kt_B[[k]][-(nsignaux+1),]
        
        M <- (1/nsignaux)*(sum(mu^2)) + solution$parametre_var$lambda_k[k]*(1+ solution$parametre$Phi_k[[k]]^2)
        
        solution$parametre$variances_W_k[[k]] <- matrix(M, nrow=d, ncol=d) 
        
      }
      
      
      
      
    }else{
      
      Vk <- mise_a_jours_V_k(solution, simulation)
      
      solution$parametre$variances_V_k <- Vk
      
      ### MISE A JOUR DE W_K ###
      
      Wk <- mise_a_jours_W_k(solution)
      solution$parametre$variances_W_k <- Wk
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
       
        print("PHI_k set to <1")
        diag(phi_k)[abs(diag(phi_k))>=1] <- 1-epsilon
      }
      
      solution$parametre$Phi_k[[k]]<- phi_k
    }
    
    
    
    ### MISE A JOUR DES COEFFICIENTS A ##
    
    d <- solution$size$d
    K <- solution$size$nbClasse 
    nsignaux <- solution$size$nsignaux
    n <- solution$size$n
    dim_U <- length(solution$coefficient)
    
    
    if(d ==1 ){
      
      
      somme_effect_square <- 0
      somme_effect <- list()
      for(k in 1:K){
        X_M <- X-solution$filtre_Kalman$mu_kt_B[[k]][-1]*matrix(1,nrow=solution$size$nsignaux, 
                                                                ncol=solution$size$n)
        
        somme_effect[[k]] <- rep(0, dim_U)
        for(t in 1:nsignaux){
          somme_effect[[k]] <- somme_effect[[k]] + rowSums(t(facteur_exo[[t]])%*%(solution$parametre_var$t_ik[,k]*t(X_M[t,])))
        }
      }
      for(t in 1:nsignaux){
        somme_effect_square <- somme_effect_square + t(facteur_exo[[t]])%*%facteur_exo[[t]]
      }
      
      
      argument_1 <- matrix(0, nrow=dim_U, ncol=dim_U)
      argument_2 <- rep(0,dim_U)
      for(k in 1:K){
        argument_1 <- argument_1 + sum(solution$parametre_var$t_ik[,k]*solution$parametre$variances_V_k[[k]][1,1])*somme_effect_square
        argument_2 <- argument_2 + solution$parametre$variances_V_k[[k]][1,1]%*%somme_effect[[k]]
      }
      
      
      solution$coefficient <-  solve(argument_1)%*%t(argument_2)
      
      
    }else{
      
      INDICE <- merge(data.frame(classe=seq(1,K)), 
                      merge(data.frame(individu=seq(1,n)), data.frame(temps=seq(1,nsignaux))))
      
      
      
      s2 <- 0
      s1 <- 0
      for(k in 1:K){
        for(t in 1:nsignaux){
          a2 <- t(facteur_exo[[t]]) %*% solve(solution$parametre$variances_V_k[[k]]) 
          a1 <- t(facteur_exo[[t]]) %*% solve(solution$parametre$variances_V_k[[k]]) %*% facteur_exo[[t]]  
          
          for(i in 1:n){
            b2 <- (X[which(REFERENCE$TEMPS==t),i] - solution$filtre_Kalman$mu_kt_B[[k]][t+1,])
            s2 <- s2 + solution$parametre_var$t_ik[i,k]*(a2%*%b2)
            s1 <- s1 + solution$parametre_var$t_ik[i,k]*a1
          } 
        } 
      }
      
      solution$coefficient <- solve(s1)%*%s2
    }
    
    
  
  
    ### CRITERE ARRET ###
    
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
  aic <- calcul_AIC(solution, elbo)

  aic_r <- calcul_AIC_r(solution, elbo)
  icl <- calcul_ICL(solution, facteur_exo = facteur_exo,X = X, REFERENCE = REFERENCE)
  
  fin <- toc(quiet = TRUE)
  
  temps_calcul <- fin$toc-debut

  
  return(list(BIC=bb, AIC = aic, AIC_R = aic_r, best_solution = solution, nbr_iteration = iter,
              temps_calcul= temps_calcul, elbo=elbo,  ICL = icl, table_erreurs = table_erreurs ))
  
  
}





compute_ELBO<- function(solution, facteur_exo, X, REFERENCE){
  
  
  d <- solution$size$d
  K <- solution$size$nbClasse
  nsignaux <- solution$size$nsignaux
  n <- solution$size$n
  
  ## Premier élément, t=0
  log_mod_1 <- 0
  log_mod_2 <- 0
  log_mod_2_b <- 0
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
        
        mu_2 <- solution$filtre_Kalman$mu_kt_B[[k]][t+1,] + facteur_exo[[t]]%*%solution$coefficient
        Sigma_2 <- matrix(solution$parametre$variances_V_k[[k]], nrow=d, ncol=d)
        z1 <- log_densite_normale(X[which(REFERENCE$TEMPS==t),i], mu_2, Sigma_2)
        z2 <- log(solution$parametre$pi_k[k])
        z3 <- -0.5*(solution$parametre_var$lambda_k[k]*trace(solve(Sigma_2)))
        z4 <- solution$parametre_var$t_ik[i,k]*(z1+ z3)
        
        log_mod_2 <- log_mod_2 + z4
        
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
  
  q <- length(solution$coefficient)

  nbr_parametre <- (5*K*d + K-1 + q)
  
  L <- ELBO
  
  BIC <- -2*L + nbr_parametre*log(nsignaux*n)
  
  return(BIC/100)
}

calcul_AIC <- function(solution, ELBO){
  
  d <- solution$size$d
  K <- solution$size$nbClasse 
  nsignaux <- solution$size$nsignaux
  n <- solution$size$n
  
  q <- length(solution$coefficient)
  nbr_parametre <- (5*K*d + K-1 + q) 
  
  L <- ELBO
  
 AIC <- -2*L + 2*nbr_parametre
  
  return(AIC/100)
}


calcul_AIC_r <- function(solution, ELBO){
  
  d <- solution$size$d
  K <- solution$size$nbClasse 
  nsignaux <- solution$size$nsignaux
  n <- solution$size$n
  
  q <- length(solution$coefficient)

  nbr_parametre <- (5*K*d + K-1 + q)
  
  L <- ELBO
  
  AIC_r <- -2*L + 2*nbr_parametre + (2*nbr_parametre*(nbr_parametre+1))/(n*nsignaux - nbr_parametre -1)
  
  return(AIC_r/100)
}




calcul_ICL <- function(solution,facteur_exo, X, REFERENCE ){
  
  d <- solution$size$d
  K <- solution$size$nbClasse 
  nsignaux <- solution$size$nsignaux
  n <- solution$size$n
  
  q <- length(solution$coefficient)
  nbr_parametre <- (5*K*d + K-1 + q) 
  ## Premier élément, t=0
  log_mod_1 <- 0
  log_mod_2 <- 0
  log_mod_3 <- 0
  
  for(k in 1:K){
    
    z1 <- log_densite_normale(solution$filtre_Kalman$mu_kt_B[[k]][1,],
                              solution$parametre$mu_0[[k]], solution$parametre$Sigma_0[[k]])
    
    log_mod_1 <- log_mod_1 + z1
    
    
    for(t in 1:nsignaux){
      ## 3
      mu_3 <- solution$parametre$Phi_k[[k]]%*%solution$filtre_Kalman$mu_kt_B[[k]][(t),]
      Sigma_3 <- solution$parametre$variances_W_k[[k]]
      
      z1 <- log_densite_normale(solution$filtre_Kalman$mu_kt_B[[k]][(t+1),],mu_3, Sigma_3)
      
      log_mod_3 <- log_mod_3 + z1
      
      ## 2
      for(i in 1:n){
        mu_2 <- solution$filtre_Kalman$mu_kt_B[[k]][t+1,] + facteur_exo[[t]]%*%solution$coefficient
        Sigma_2 <- matrix(solution$parametre$variances_V_k[[k]], nrow=d, ncol=d)
        z1 <- log_densite_normale(X[which(REFERENCE$TEMPS==t),i], mu_2, Sigma_2)
  
        z2 <- log(solution$parametre$pi_k[k])
        log_mod_2 <- log_mod_2 + z1 + z2
      }
    }
  }
  
  Log_vrais_complete <- log_mod_1 + log_mod_2 + log_mod_3 
  
  #####################################################
  
  
  z_ik_map <- round(solution$parametre_var$t_ik)
  log_tau <- log(solution$parametre_var$t_ik)
  
  ICL <- Log_vrais_complete + sum(z_ik_map*log_tau) - 0.5*log(n*nsignaux)*nbr_parametre
  
  
  return(ICL)
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
  
  graphique_facteur_processus(solution, simulation_VEM, chemin_sauvegarde)
  
  graphique_effet_facteur(solution, simulation_VEM, chemin_sauvegarde)
  
}


algorithme_VEM <- function(simulation, K, TT){
  solution <- initialisation_VEM_Kmeans(X= simulation$data_X,
                                        K=K, facteur_exo=simulation$data_facteurs,
                                        REFERENCE=simulation$REFERENCE, nsignaux = TT)

  
  res <- iterations_VEM(simulation, K, TT, solution)
  
  
  return(res)
}


mise_a_jours_poids_2 <- function(solution, facteur_exo, X, REFERENCE){
  
  d <- solution$size$d
  K <- solution$size$nbClasse 
  nsignaux <- solution$size$nsignaux
  n <- solution$size$n
  
  
  numerateur <- matrix(0,nrow=n,ncol=K) 
  
  for(k in 1:K){
    
    a <- nsignaux*(log(solution$parametre$pi_k[k]) -
                     0.5*solution$parametre_var$lambda_k[k]*
                     trace(solve(solution$parametre$variances_V_k[[k]])))
    
    for(i in 1:n){
      sum_t <- 0 
      for(t in 1:nsignaux){
        zz <- log_densite_normale(X[which(REFERENCE$TEMPS == t),i],
                                  (solution$filtre_Kalman$mu_kt_B[[k]][(t+1),] +facteur_exo[[t]]%*%solution$coefficient), 
                                  matrix(solution$parametre$variances_V_k[[k]], nrow=d, ncol=d) )
        
        sum_t <- sum_t + zz 
      }
      
      numerateur[i,k] <- a + sum_t
      
    }
    
    
  }
  
  M <-exp(numerateur-matrix(rep(rowMax(numerateur),K),nrow=n, ncol=K))
  MM <- apply(M, c(1, 2), function(x) max(x,exp(-50)))
  
  
  solution$parametre_var$t_ik <- normalisation_poids(MM)
  
  return(solution)
}

mise_a_jours_V_k <- function(solution, simulation){
  
  facteur_exo <- simulation$data_facteurs
  X <- simulation$data_X
  REFERENCE <- simulation$REFERENCE
  
  
  d <- solution$size$d
  K <- solution$size$nbClasse 
  nsignaux <- solution$size$nsignaux
  n <- solution$size$n
  INDICE  <- merge(data.frame(individu = seq(1,n)), data.frame(temps = seq(1,nsignaux)))
  
  f <- function(k,d){
    
    a <- 1/(nsignaux*sum(solution$parametre_var$t_ik[,k]))
    
    if(d > 1){
      V_k <- diag(diag(a*rowSums(mapply(function(i,t,k) 
        solution$parametre_var$t_ik[i,k]*((X[which(REFERENCE$TEMPS == t),i] 
                                           - solution$filtre_Kalman$mu_kt_B[[k]][t+1,] 
                                           - facteur_exo[[t]]%*%solution$coefficient )
                                          %*%t((X[which(REFERENCE$TEMPS == t),i] 
                                                - solution$filtre_Kalman$mu_kt_B[[k]][t+1,]
                                                - facteur_exo[[t]]%*%solution$coefficient))
                                          + solution$parametre_var$lambda_k[k]*diag(d)),
        INDICE$individu, INDICE$temps, k, SIMPLIFY = "array"),dims=2)))
      
      if(sum(diag(V_k)==0)){
        diag(V_k)[diag(V_k)==0] <- exp(-100)
      } ## Pour s'assurer que la matrice soit inversible 
      
    } else {
      V_k <- a*sum(mapply(function(i,t,k) 
        solution$parametre_var$t_ik[i,k]*((X[which(REFERENCE$TEMPS == t),i] 
                                           - solution$filtre_Kalman$mu_kt_B[[k]][t+1,]
                                           - facteur_exo[[t]]%*%solution$coefficient)%*%t((X[which(REFERENCE$TEMPS == t),i]- solution$filtre_Kalman$mu_kt_B[[k]][t+1,]- facteur_exo[[t]]%*%solution$coefficient))+ solution$parametre_var$lambda_k[k]*diag(d)), 
        INDICE$individu, INDICE$temps, k, SIMPLIFY = "array"))
      
      if(V_k == 0 ){
        V_k <- exp(-100)
      } 
    }
    
    
    return(V_k)
  }
  
  V_k <- lapply(seq(1,K), f,d)
  
  return(V_k)
}


mise_a_jours_W_k <- function(solution){
  d <- 1
  K <- solution$size$nbClasse 
  nsignaux <- solution$size$nsignaux
  n <- solution$size$n
  W_k <- list()
  for(k in 1:K){
    
    somme_t <- 0
    for(t in 2:(nsignaux+1)){
      
      z <- solution$filtre_Kalman$mu_kt_B[[k]][t,] - 
        solution$parametre$Phi_k[[k]]%*%solution$filtre_Kalman$mu_kt_B[[k]][t-1,]
      zz <- z%*%t(z) 
      somme_t <- somme_t + zz
      
    }
    
    W_k[[k]] <- diag(diag((1/nsignaux)*somme_t),d)+ solution$parametre_var$lambda_k[k]*(diag(1,d)+ solution$parametre$Phi_k[[k]]%*%solution$parametre$Phi_k[[k]])
  }
  
  return(W_k)
}



library(cluster)
silhouette_score <- function(cluster, df){
  ss <- silhouette(cluster, dist(df))
  return(list(list_silhouette = ss[,3], S = mean(ss[, 3])))
}



compute_silhouette <- function(solution, simulation){

  classes <- apply(solution$parametre_var$t_ik, 1, which.max)
  
  res <- silhouette(classes, dist(t(simulation$data_X)))
  
  return(list(list_silhouette = res[,3], S = mean(res[, 3])))
}


log_vrais_indiv_KM <- function(i, solution, simulation, sigma){
  
  nsignaux <- solution$size$nsignaux
  K <- solution$size$nbClasse
  x_it <- simulation$data_X[,i]
  t_ik <- solution$parametre_var$t_ik[i,]
  b_kt <- solution$filtre_Kalman$mu_kt_B
  
  ll <- rep(0, nsignaux)
  
  effet_ext <- rep(0,nsignaux)
  for(t in 1:nsignaux){
    effet_ext[t] <- as.numeric(simulation$data_facteurs[[t]]%*%solution$coefficient)
  }
  
  for(k in 1:K){
    
    ll <- ll +  as.numeric(t_ik[k])*(as.numeric((1/sigma))*(-(x_it - b_kt[[k]][-1,]- effet_ext)^2) - log(as.numeric(sigma)))
  }

  
  return(ll) 
}




compute_log_vrais_KM <- function(solution, simulation){
  
  nsignaux <- solution$size$nsignaux
  n <- solution$size$n
  
  effet_ext <- rep(0,nsignaux)
  for(t in 1:nsignaux){
    effet_ext[t] <- as.numeric(simulation$data_facteurs[[t]]%*%solution$coefficient)
  }
  classes <- apply(solution$parametre_var$t_ik, 1, which.max)
  residus <- matrix(nrow=nsignaux, ncol=n)
  
  for(i in 1:solution$size$n){
    residus[,i] <- simulation$data_X[,i] - solution$filtre_Kalman$mu_kt_B[[classes[i]]][-1] - effet_ext
  }
  
  sigma <- (1/((n*nsignaux) - 1))*sum((residus - mean(residus))^2)    
  
  
  table_vrais <- mapply(i= seq(1, solution$size$n), log_vrais_indiv_KM, 
                        MoreArgs= list(solution , simulation, sigma))
  
  mean_vraisemblance <- (1/(solution$size$n*solution$size$nsignaux))*sum(table_vrais)
  
  print(mean_vraisemblance)
  
  return(mean_vraisemblance)
}


log_vrais_indiv <- function(i, solution, simulation){
  nsignaux <- solution$size$nsignaux
  K <- solution$size$nbClasse
  x_it <- simulation$data_X[,i]
  t_ik <- solution$parametre_var$t_ik[i,]
  
  b_kt <- solution$filtre_Kalman$mu_kt_B
  
  effet_ext <- rep(0,nsignaux)
  for(t in 1:nsignaux){
    effet_ext[t] <- as.numeric(simulation$data_facteurs[[t]]%*%solution$coefficient)
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
  
  print(mean_vraisemblance)
  return(mean_vraisemblance)
}


log_vrais_indiv_EM <- function(i, solution, simulation){
  
  nsignaux <- solution$size$nsignaux
  K <- solution$size$nbClasse
  x_it <- simulation$data_X[,i]
  t_ik <- solution$parametre_var$t_ik[i,]
  b_kt <- solution$filtre_Kalman$mu_kt_B
  sigma <- solution$parametre$sigma
  ll <- rep(0, nsignaux)
  
  effet_ext <- rep(0,nsignaux)
  for(t in 1:nsignaux){
    effet_ext[t] <- as.numeric(simulation$data_facteurs[[t]]%*%solution$coefficient)
  }
  
  for(k in 1:K){
    
    ll <- ll +  as.numeric(t_ik[k])*(as.numeric((1/sigma[[k]]))*(-(x_it - b_kt[[k]][-1,]- effet_ext)^2) - log(as.numeric(sigma[[k]])))
  }
  
  
  return(ll) 
}




compute_log_vrais_EM <- function(solution, simulation){
  
  nsignaux <- solution$size$nsignaux
  n <- solution$size$n

  table_vrais <- mapply(i= seq(1, solution$size$n), log_vrais_indiv_EM, 
                        MoreArgs= list(solution , simulation))
  
  mean_vraisemblance <- (1/(solution$size$n*solution$size$nsignaux))*sum(table_vrais)
  
  print(mean_vraisemblance)
  
  return(mean_vraisemblance)
}






