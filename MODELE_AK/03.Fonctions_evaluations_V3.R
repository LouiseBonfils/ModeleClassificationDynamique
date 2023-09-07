erreur_moyenne_centres <- function(solution, simulation){
  
  K <- solution$size$nbClasse
  
  table_corresp_classes <- create_table_corresp_classes(solution, simulation )
  d <- solution$size$d
  nsignaux <- solution$size$nsignaux
  
  b_simulation <- simulation$data_processus
  b_estimation <- solution$filtre_Kalman$mu_kt_B
  
  table_processus <- data.frame(simulation=NA, estimation=NA, classe=NA, temps=NA, dim=NA)
  
  for(j in 1:d){
    for(k in 1:K){
      indice <- as.numeric(as.character(table_corresp_classes$estimation[table_corresp_classes$simulation==k]))
      table_processus <- rbind(table_processus, data.frame(simulation = as.matrix(b_simulation[[k]], nrow=nsignaux)[-1,j], 
                                                           estimation = b_estimation[[indice]][-1,j] , 
                                                           classe=k, temps = seq(1,nsignaux), dim=j ))   
    }
    
  }  
  
  table_difference <- table_processus$simulation[-1] - table_processus$estimation[-1] 
  
  mean_error <- mean(table_difference^2)
  
  return(mean_error)  
}




create_table_corresp_classes <- function(solution, simulation){
  
  vec_classes_estimees <- apply(solution$parametre_var$t_ik, 1, which.max)
  vec_classes_simulees <- simulation$classes
  
  table_corresp_classes <- table(vec_classes_estimees, vec_classes_simulees) %>%
    as.data.frame %>%
    rename(simulation = vec_classes_simulees, 
           estimation = vec_classes_estimees)
  
  K <- length(unique(simulation$classes))
  
  
  if(length(unique(table_corresp_classes$estimation)) < K){
    
    table_corresp_classes <- table(vec_classes_estimees, vec_classes_simulees) %>%
      as.data.frame %>%
      rename(simulation = vec_classes_simulees, 
             estimation = vec_classes_estimees)%>%
      arrange(desc(Freq)) %>%
      top_n(K)%>%
      filter(Freq !=0) %>%
      dplyr::select(-Freq)
    
    if(dim(table_corresp_classes)[1]>K){
      table_corresp_classes <- table_corresp_classes[1:K, ]
    }
    
    return(table_corresp_classes)
    
  } else {
    
    combi_possible = permutations(n=K, r=K, v=1:K)
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
      rename("simulation" = "base", 
             "estimation" = "V1")
    
    
    
    
    return(combi_correspondante)
  }
}

taux_bon_classe_2 <- function(classe_bdd, proba_appartenance){
  
  classe_estime <- apply(proba_appartenance, 1, which.max)
  K = max(classe_bdd)
  
  table_corresp_classes <- create_table_corresp_classes_proba(proba_appartenance, classe_bdd)

  classe_estime_c <-classe_estime %>% 
    as.data.frame%>%
    rename('classe'='.')%>%
    mutate(classe = ifelse(classe == 1, table_corresp_classes$simulation[table_corresp_classes$estimation == 1], 
                           table_corresp_classes$simulation[table_corresp_classes$estimation != 1]))
  
  
  taux <- (sum(classe_estime_c == classe_bdd)/length(classe_bdd))
  
  
  return(taux)
  
  
}

taux_bon_classe_4 <- function(classe_bdd, proba_appartenance){
  
  classe_estime <- apply(proba_appartenance, 1, which.max)
  K = max(classe_bdd)
  
  table_corresp_classes <- create_table_corresp_classes_proba(proba_appartenance, classe_bdd)

  
  classe_estime_c <-classe_estime %>% 
    as.data.frame%>%
    rename('classe'='.')%>%
    mutate(classe = ifelse(classe == 1, table_corresp_classes$simulation[table_corresp_classes$estimation 
                                                                         == 1], 
                           ifelse(classe==2, table_corresp_classes$simulation[table_corresp_classes$estimation 
                                                                              == 2], 
                                  ifelse(classe==3, table_corresp_classes$simulation[table_corresp_classes$estimation 
                                                                                     == 3],
                                         table_corresp_classes$simulation[table_corresp_classes$estimation == 4]))))
  
  
  taux <- (sum(classe_estime_c == classe_bdd)/length(classe_bdd))
  
  
  return(taux)
  
  
}

erreur_facteurs_centres <- function(solution, simulation){
  
  K <- solution$size$nbClasse
  
  table_corresp_classes <- create_table_corresp_classes(solution, simulation )
  d <- solution$size$d
  nsignaux <- solution$size$nsignaux
  
  effet_simulation <- simulation$data_processus
  effet_estimation <- solution$filtre_Kalman$mu_kt_B
  b_simulation <- simulation$data_processus
  b_estimation <- solution$filtre_Kalman$mu_kt_B
  
  table_effet <- data.frame(simulation=NA, estimation=NA, classe=NA, temps=NA, dim=NA)
  
  for(j in 1:d){
    for(k in 1:K){
      indice <- as.numeric(as.character(table_corresp_classes$estimation[table_corresp_classes$simulation==k]))
      table_effet <- rbind(table_effet, data.frame(simulation = 
                                                     sapply(seq(1,nsignaux), 
                                                            function(t,j)
                                                              (simulation$data_facteur[[t]]%*%simulation$coefficient[[k]])[j,] +
                                                              as.matrix(b_simulation[[k]], nrow=nsignaux)[(t+1),j], j), 
                                                   estimation = 
                                                     sapply(seq(1,nsignaux), function(t,j)
                                                       (simulation$data_facteur[[t]]%*%
                                                          solution$coefficient[[indice]])[j,] + b_estimation[[indice]][(t+1),j] , j),
                                                   classe=k, temps = seq(1,nsignaux), dim=j ))   
    }
    
  }  
  
  table_error <- table_effet[-1,] %>%
    group_by(dim)%>%
    summarise(error = mean((simulation-estimation)^2))%>%
    as.data.frame
  
  return(mean(table_error$error))
  
}



erreur_effet_facteurs <- function(simulation, solution){
  K <- solution$size$nbClasse
  
  table_corresp_classes <- create_table_corresp_classes(solution, simulation )
  d <- solution$size$d
  nsignaux <- solution$size$nsignaux
  
  effet_simulation <- simulation$data_processus
  effet_estimation <- solution$filtre_Kalman$mu_kt_B
  
  
  table_effet <- data.frame(simulation=NA, estimation=NA, classe=NA, temps=NA, dim=NA)
  
  for(j in 1:d){
    for(k in 1:K){
      indice <- as.numeric(as.character(table_corresp_classes$estimation[table_corresp_classes$simulation==k]))
      table_effet <- rbind(table_effet, data.frame(simulation = 
                                                     sapply(seq(1,nsignaux), 
                                                            function(t,j)
                                                              (simulation$data_facteur[[t]]
                                                               %*%simulation$coefficient[[k]])[j,], j), 
                                                   estimation = 
                                                     sapply(seq(1,nsignaux), function(t,j)
                                                       (simulation$data_facteur[[t]]%*%
                                                          solution$coefficient[[indice]])[j,], j),
                                                   classe=k, temps = seq(1,nsignaux), dim=j ))   
    }
    
  }  
  
  table_error <- table_effet[-1,] %>%
    group_by(dim)%>%
    summarise(error = mean((simulation-estimation)^2))%>%
    as.data.frame
  
  
  return(mean(table_error$error))
}


create_table_corresp_classes_proba <- function(proba_appartenance, classe_bdd){
  
  vec_classes_estimees <- apply(proba_appartenance, 1, which.max)
  vec_classes_simulees <- classe_bdd
  
  K <- length(unique(classe_bdd))
  
  table <- table(vec_classes_estimees, vec_classes_simulees) %>%
    as.data.frame %>%
    rename(simulation = vec_classes_simulees, 
           estimation = vec_classes_estimees)
  
  
  if(length(unique(table$estimation)) < K){ ## Cas où au moins une classe est vide 
    
    table_corresp_classes <- table(vec_classes_estimees, vec_classes_simulees) %>%
      as.data.frame %>%
      rename(simulation = vec_classes_simulees, 
             estimation = vec_classes_estimees)%>%
      arrange(desc(Freq)) %>%
      top_n(K)%>%
      filter(Freq !=0) %>%
      dplyr::select(-Freq) 
    
    if(length(unique(table_corresp_classes$estimation)) != length(unique(table$estimation))){
      
      diff <-abs(length(unique(table_corresp_classes$estimation))-length(unique(table$estimation)))
      
      table_corresp_classes <- table(vec_classes_estimees, vec_classes_simulees) %>%
        as.data.frame %>%
        rename(simulation = vec_classes_simulees, 
               estimation = vec_classes_estimees)%>%
        arrange(desc(Freq)) %>%
        top_n(K+diff)
      
      use_while=F
      while(length(unique(table_corresp_classes$estimation[1:(K+diff)])) != length(unique(table$estimation))){
        diff <- diff+1
        table_corresp_classes <- table(vec_classes_estimees, vec_classes_simulees) %>%
          as.data.frame %>%
          rename(simulation = vec_classes_simulees, 
                 estimation = vec_classes_estimees)%>%
          arrange(desc(Freq)) %>%
          top_n(K+diff)
        print("WHILE")      
        use_while <- T
      }
      
      
      frequence_classe_estimation <- table(table_corresp_classes$estimation)%>%
        as.data.frame%>%
        arrange(Freq)
      
      position <- rep(NA, K)
      pp <- 1
      for(i in 1:dim(frequence_classe_estimation)[1]){
        sum_freq <- sum(frequence_classe_estimation$Freq)
        if(sum_freq == 5){
          classe <- frequence_classe_estimation$Var1[i]
          position[pp]<- which(table_corresp_classes$estimation==classe)[1]
          pp <- pp+1
        }
      }
      
      if(sum(is.na(position))!=0){
        count <- sum(is.na(position))
        classe <- frequence_classe_estimation$Var1[frequence_classe_estimation$Freq!=1]
        j=1
        for(c in 1:length(classe)){
          cc <- classe[c]
          if(j <= count){
            position[pp] <- which(table_corresp_classes$estimation==classe)[j+1]
            pp <- pp+1
            j <- j+1
          }
        }
      }
      
      table_corresp_classes <- table_corresp_classes[position, ]
      
      
      
     
      
    }
    return(table_corresp_classes)
    
  } else {
    
    combi_possible = permutations(K, K, 1:K)
    base <- seq(1,K)
    
    table_trace <- numeric(dim(combi_possible)[1])
    
    for(i in 1:dim(combi_possible)[1] ){
      
      combi <- combi_possible[i,]
      
      corresp <- cbind(base, combi) %>%
        as.data.frame
      
      
      t <- merge(table, corresp,
                 by.x=c("simulation", "estimation"),
                 by.y=c("base", "combi"))
      
      table_trace[i] <- sum(t$Freq)
      
    }
    
    combi_correspondante <-cbind(combi_possible[which.max(table_trace),], base) %>%
      as.data.frame() %>%
      rename("simulation" = "base", 
             "estimation" = "V1")
    
    return(combi_correspondante)
  }
}

erreur_moyenne_centres_mixe_reg <- function(solution, simulation){
  
  K <- solution$size$nbClasse
  
  table_corresp_classes <- create_table_corresp_classes(solution, simulation)
  
  
  
  d <- solution$size$d
  nsignaux <- solution$size$nsignaux
  
  b_simulation <- simulation$data_processus
  b_estimation <- solution$filtre_Kalman$mu_kt_B
  
  table_processus <- data.frame(simulation=NA, estimation=NA, classe=NA, temps=NA, dim=NA)
  
  for(j in 1:d){
    for(k in 1:K){
      indice <- as.numeric(as.character(table_corresp_classes$estimation[table_corresp_classes$simulation==k]))
      table_processus <- rbind(table_processus, data.frame(simulation = as.matrix(b_simulation[[k]], nrow=nsignaux)[-1,j], 
                                                           estimation = b_estimation[[indice]][,j] , 
                                                           classe=k, temps = seq(1,nsignaux), dim=j ))   
    }
    
  }  
  
  table_difference <- table_processus$simulation[-1] - table_processus$estimation[-1] 
  
  mean_error <- mean(table_difference^2)
  
  return(mean_error)  
}


erreur_facteurs_centres_mix_reg <- function(solution, simulation){
  
  K <- solution$size$nbClasse
  
  table_corresp_classes <- create_table_corresp_classes(solution, simulation )
  d <- solution$size$d
  nsignaux <- solution$size$nsignaux
  
  effet_simulation <- simulation$data_processus
  effet_estimation <- solution$filtre_Kalman$mu_kt_B
  b_simulation <- simulation$data_processus
  b_estimation <- solution$filtre_Kalman$mu_kt_B
  
  table_effet <- data.frame(simulation=NA, estimation=NA, classe=NA, temps=NA, dim=NA)
  
  for(j in 1:d){
    for(k in 1:K){
      indice <- as.numeric(as.character(table_corresp_classes$estimation[table_corresp_classes$simulation==k]))
      table_effet <- rbind(table_effet, data.frame(simulation = 
                                                     sapply(seq(1,nsignaux), 
                                                            function(t,j)
                                                              (simulation$data_facteur[[t]]%*%simulation$coefficient[[k]])[j,] +
                                                              as.matrix(b_simulation[[k]], nrow=nsignaux)[(t+1),j], j), 
                                                   estimation = 
                                                     sapply(seq(1,nsignaux), function(t,j)
                                                       (simulation$data_facteur[[t]]%*%
                                                          solution$coefficient[[indice]])[j,] + b_estimation[[indice]][(t),j] , j),
                                                   classe=k, temps = seq(1,nsignaux), dim=j ))   
    }
    
  }  
  
  table_error <- table_effet[-1,] %>%
    group_by(dim)%>%
    summarise(error = mean((simulation-estimation)^2))%>%
    as.data.frame
  
  return(mean(table_error$error))
  
}

