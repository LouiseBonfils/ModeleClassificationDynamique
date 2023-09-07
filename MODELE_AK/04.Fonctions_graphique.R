graphique_simulation_effet_facteur <- function(simulation, path_resultat = './', sub=''){
  
  K <- length(simulation$coefficient)
  
  vec_classes_simulees <- simulation$classes
  
  
  
  table_effet <- data.frame(simulation=NA, dim= NA , temps=NA, classe=NA)

  d <- 1
  for(k in 1:K){
    for(j in 1:d){
      table_effet <- rbind(table_effet, data.frame(simulation = sapply(seq(1,nsignaux), function(t,j) (simulation$data_facteur[[t]]%*%simulation$coefficient[[k]])[j,], j), 
                                                    dim=j, temps = seq(1,nsignaux), classe=k ))   
    }
  }
  
  table_effet <- table_effet[-1,]
  
  gg <- table_effet %>%
    gather(key='type', value='valeurs', -temps, -dim, -classe) %>%
    ggplot(aes(x=temps, y=valeurs, color=as.factor(classe))) + geom_line(size=0.8)+ 
    labs(title = "Effet des facteurs exogènes simulé", subtitle = sub) +
    theme_bw() + scale_color_manual(values=cols)
  
  ggsave(paste0(path_resultat, "effet_facteurs_simulation.png"),plot=gg, width = 15, height = 10, units = "cm")
  
  print(gg)
}


cols <- c("#FAAB18", "#1380A1","#990000", "#588300", "black",
          'red', 'yellow', 'green', 'pink', 'grey',
          "yellowgreen", "aquamarine1","chocolate3", "brown2" , "darkorchid3")


graphique_simulation <- function(X,REFERENCE, d, path_essais='./', sub=""){
  
  XX <- data.frame(temps = NA, dim=NA, individu=NA, x=NA)
  
  for(j in 1:d){
    XX <- rbind(XX, X[which(REFERENCE$DIM==j), ] %>%
                  as.data.frame %>%
                  mutate(temps = seq(1, (dim(X)[1]/d)), dim=j) %>%
                  gather(key="individu", value="x", -temps, -dim)
    )
  }
  
 gg <-  XX[-1,] %>%
    ggplot(aes(x=temps, y=x, group=individu)) + geom_line(size=0.8) + theme_bw() + 
    facet_grid(dim~., labeller = label_both) + 
    theme(legend.position = "none") + labs(x="Temps", y="Données simulées", title = "Jeu de données simulées, dimension", subtitle = sub)

   
  ggsave(paste0(path_essais, "simulation_data_X.png"), plot=gg, width = 15, height = 10, units = "cm")
  
}
#### FONCTION AFFICHAGE DE RESULTATS  


## 1Heatmap des classes 

graphique_matrice_classes <- function(solution, simulation, path_resultat='./', subtitle= ""){
  
  
  vec_classes_estimees <- apply(solution$parametre_var$t_ik, 1, which.max)
  vec_classes_simulees <- simulation$classes
  
  
  gg <- table(vec_classes_estimees, vec_classes_simulees) %>%
    as.data.frame %>%
    rename(effectifs = Freq) %>%
    ggplot(aes(x=vec_classes_estimees, y=vec_classes_simulees, fill=effectifs)) + 
    geom_tile(colour = "grey50") + scale_fill_gradient2(low = "white", high = "tomato3") +
    geom_text(aes(label=effectifs)) +
    labs(x="Classes estimées", y="Classes simulées", 
         title = "Effectifs croisés entre les classes estimées et simulées", 
         subtitle = subtitle)  + 
    theme_bw()  +  theme(legend.position = "none") 
  
  ggsave(paste0(path_resultat, "matrice_classes.png"),plot=gg, width = 8, height = 8, units = "cm") 
  
  
  print(gg)
}


##  Tracer les processus : 

graphique_proccesus <- function(solution, simulation , path_resultat='./', subtitle = ""){
  
  K <- solution$size$nbClasse
  vec_classes_estimees <- apply(solution$parametre_var$t_ik, 1, which.max)
  vec_classes_simulees <- simulation$classes
  
  table_corresp_classes <- create_table_corresp_classes(solution, simulation)
  
  d <- solution$size$d
  nsignaux <- solution$size$nsignaux
  K <- solution$size$nbClasse
  
  b_simulation <- simulation$data_processus
  b_estimation <- solution$filtre_Kalman$mu_kt_B
  
  
  table_processus <- data.frame(simulation=NA, estimation=NA, classe=NA, temps=NA, dim=NA)
  
  
  
  for(j in 1:d){
    for(k in 1:K){
      table_processus <- rbind(table_processus, data.frame(simulation = b_simulation[[k]][-1,j], 
                                                           estimation = b_estimation[[table_corresp_classes$estimation[table_corresp_classes$simulation==k]]][-1,j] , 
                                                           classe=k, temps = seq(1,nsignaux), dim=j ))   
    }
    
  }
  
  gg <-  table_processus[-1,] %>%
    gather(key='type', value='centres', -temps, -classe, -dim) %>%
    ggplot(aes(x=temps, y=centres, color=type)) + geom_line(size=0.8)+ 
    facet_grid(classe~dim, scale='free', labeller = label_both)+
    labs(title = "Centres de classes estimés et simulés", 
         subtitle = subtitle)+
    # , 
    # subtitle="Cas sans effet croisé") +
    theme_bw()+ scale_color_manual(values=c('orangered', 'indianred4'))
  
  ggsave(paste0(path_resultat, "Process_b.png"), plot=gg, width = 15, height = 10, units = "cm")
  
  print(gg)
}

############################################################################

## Données : 


estimation_X <- function(solution, simulation){
  
  classes <- apply(solution$parametre_var$t_ik, 1, which.max)
  n <- solution$size$n
  d <- solution$size$d
  nsignaux <- solution$size$nsignaux
  
  X_hat <- sapply(seq(1,n), function(i,solution, classes) cbind(sapply(seq(1,nsignaux), function(t, i, classes) 
    simulation$data_facteur[[t]][,-(1:d)]%*%solution$coefficient[[classes[i]]][-(1:d)] + 
      solution$filtre_Kalman$mu_kt_B[[classes[i]]][(t+1),], i, classes, simplify = "array")), solution, classes)
  
  colnames(X_hat) <- colnames(simulation$data_X)
  
  
  return(list(X_hat=X_hat, classes=classes))
}


graphique_data_estimee <- function(solution, simulation, path_resultat='./'){
  
  X_hat <- estimation_X_melange(solution, simulation)$X_hat
  X <- simulation$data_X
  REFERENCE <- simulation$REFERENCE
  
  d <- solution$size$d
  
  nsignaux <- solution$size$nsignaux
  
  
  f <- function(j, solution, simulation, REFERENCE){
    
    table_data <- rbind( X_hat[REFERENCE$DIM==j,] %>%
                           as.data.frame %>%
                           mutate(temps=seq(1,nsignaux), type="estimation"), 
                         X[REFERENCE$DIM ==j, ]%>%
                           as.data.frame %>%
                           mutate(temps=seq(1,nsignaux), type="simulation"))
    
    g2 <- table_data %>%
      gather(key="individu", value='donnees', -temps, -type)%>%
      mutate(individu = as.numeric(substr(individu, 2, nchar(individu)))) %>%
      ggplot(aes(x=temps, y=donnees, color=type)) + geom_line(size=0.8) + facet_wrap(.~individu, ncol=7, scale='free') +
      theme_bw() 
    ggsave(paste0(path_resultat, "data_x_wrap_",j,".png"), plot=g2, width = 60, height = 30, units = "cm")
    
    print(g2)
    
    g1 <- table_data %>%
      gather(key="individu", value='donnees', -temps, -type)%>%
      mutate(individu = as.numeric(substr(individu, 2, nchar(individu)))) %>%
      ggplot(aes(x=temps, y=donnees, color=type, linetype=type, size=type)) +  geom_line(aes(group=interaction(individu,type)))+
      theme_bw() + scale_color_manual(values=c('orangered', 'indianred4')) +scale_size_manual(values=c(1, 0.7))
    
    ggsave(paste0(path_resultat, "data_x_",j,".png"), plot=g1, width = 30, height = 20, units = "cm")
    print(g1)
  }
  
  
  sapply(seq(1,d), f, solution, simulation, REFERENCE)
}



#### GRAPHIQUE Uta + b_kt 


graphique_facteur_processus_V3 <- function(solution, simulation, path_resultat="./"){
  
  K <- solution$size$nbClasse
  vec_classes_estimees <- apply(solution$parametre_var$t_ik, 1, which.max)
  vec_classes_simulees <- simulation$classes
  
  table_corresp_classes <- create_table_corresp_classes(solution, simulation)
  
  d <- solution$size$d
  nsignaux <- solution$size$nsignaux
  K <- solution$size$nbClasse
  
  
  
  b_simulation <- simulation$data_processus
  b_estimation <- solution$filtre_Kalman$mu_kt_B
  
  
  table_processus <- data.frame(simulation=NA, estimation=NA, classe=NA, temps=NA, dim=NA)
  nsignaux <- solution$size$nsignaux
  
  
  d <-solution$size$d 

  
  
  for(j in 1:d){
    
    for(k in 1:K){
      
      table_facteur_simu <- sapply(seq(1,nsignaux), function(t,j,k) (simulation$data_facteur[[t]] %*% simulation$coefficient[[k]])[j,], j,k)
      
      table_facteur_estimation <- sapply(seq(1,nsignaux), function(t,j,k) (simulation$data_facteur[[t]]%*%solution$coefficient[[table_corresp_classes$estimation[table_corresp_classes$simulation==k]]]), j,k)
      
      table_processus <- rbind(table_processus, data.frame(simulation = b_simulation[[k]][-1,j] + table_facteur_simu, 
                                                           estimation = b_estimation[[table_corresp_classes$estimation[table_corresp_classes$simulation==k]]][-1,j] + table_facteur_estimation,
                                                           classe=k, temps = seq(1,nsignaux), dim=j ))  
     
      
      }
    
  }
  table_processus <- table_processus[-1,]
  
  gg <- table_processus %>%
    gather(key='type', value='valeurs', -temps, -classe, -dim) %>%
    ggplot(aes(x=temps, y=valeurs, color=type)) + geom_line(size=0.8)+ 
    facet_grid(classe~dim, scale='free', labeller = label_both)+
    labs(title = "Processus b et effets exogènes estimés et simulés, dimension ") +
    theme_bw() + scale_color_manual(values=c('orangered', 'indianred4'))
  
  ggsave(paste0(path_resultat, "facteurs_exo_process_b.png"),plot=gg, width = 15, height = 10, units = "cm")
  
  print(gg)
}




graphique_effet_facteur_V3 <- function(solution, simulation, path_resultat='./', subtitle = ""){
  
  K <- solution$size$nbClasse
  vec_classes_estimees <- apply(solution$parametre_var$t_ik, 1, which.max)
  vec_classes_simulees <- simulation$classes
  
  table_corresp_classes <- create_table_corresp_classes(solution, simulation)
  
  
  table_effet <- data.frame(simulation=NA, estimation=NA,dim= NA , temps=NA, classe=NA)
  nsignaux <- solution$size$nsignaux
  d <-solution$size$d 
  for(k in 1:K){
    
    
    for(j in 1:d){
      table_effet <- rbind(table_effet, data.frame(simulation = sapply(seq(1,nsignaux), function(t,j) (simulation$data_facteur[[t]]%*%simulation$coefficient[[k]])[j,], j), 
                                                   estimation =sapply(seq(1,nsignaux), function(t,j) (simulation$data_facteur[[t]][,]%*%solution$coefficient[[table_corresp_classes$estimation[table_corresp_classes$simulation==k]]])[j,], j) ,
                                                   dim=j, temps = seq(1,nsignaux), classe=k ))   
    }
  }
  
  table_effet <- table_effet[-1,]
  
  gg <- table_effet %>%
    gather(key='type', value='valeurs', -temps, -dim, -classe) %>%
    ggplot(aes(x=temps, y=valeurs, color=type)) + geom_line(size=0.8)+ 
    facet_grid(classe~dim, labeller = label_both, scales = 'free')+
    labs(title = "Effets exogènes estimés et simulés", subtitle = subtitle) +
    theme_bw() + scale_color_manual(values=c('orangered', 'indianred4'))
  
  ggsave(paste0(path_resultat, "effet_facteurs_.png"),plot=gg, width = 15, height = 10, units = "cm")
  
  print(gg)
}


graph_centre_B_F <- function(solution){
  table_backward <- data.frame(centres=NA, classe=0, dimension = 0, type='backward', temps=NA)
  table_forward <- data.frame(centres=NA, classe=0, dimension = 0, type='forward', temps=NA)
  for(k in 1:K){
    for(j in 1:d){
      table_backward <- rbind(table_backward,   solution$filtre_Kalman$mu_kt_B[[k]][,j] %>%
                                as.data.frame %>%
                                rename('centres'= '.')%>%
                                mutate(classe = k, dimension=j, type='backward', temps = seq(0, nsignaux)))
      table_forward <- rbind(table_forward,   solution$filtre_Kalman$mu_kt_F[[k]][,j] %>%
                               as.data.frame %>%
                               rename('centres'= '.')%>%
                               mutate(classe = k, dimension=j, type='forward', temps = seq(0, nsignaux)))
      
      
    }
  }
  
  gg <- rbind(table_backward[-1,], table_forward[-1,]) %>%
    ggplot(aes(x=temps, y=centres, color=type, linetype=type)) + geom_point()+geom_line() + 
    facet_grid(classe~dimension, labeller = label_both) + theme_bw() +
    labs(title='Centres estimés via le filtre de Kalman')
  print(gg)
}

graphique_processus_simulation <- function(processus, path_res='./',sub=""){
  
  K = length(processus)
  
  if(is.null(dim(processus[[1]]))){
    d <- 1
    nsignaux <- length(processus[[1]]) -1
  } else {
    d <- dim(processus[[1]])[2]
    nsignaux <- dim(processus[[1]])[1] -1 
  }
  
  table_processus <- data.frame(simulation=NA, classe=NA, temps=NA, dim=NA)
  
  for(j in 1:d){
    for(k in 1:K){
      table_processus <- rbind(table_processus, data.frame(simulation = as.matrix(processus[[k]], nrow=nsignaux+1, ncol=d)[-1,j], 
                                                           classe=k, temps = seq(1,nsignaux), dim=j ))   
    }
  }
  
  gg <- table_processus[-1,] %>%
    mutate(classe = as.factor(classe)) %>%
    rename('Classes' = classe) %>%
    ggplot(aes(x=temps, y=simulation, color=Classes)) + geom_line(size=0.8) + facet_grid(dim~.) +
    theme_bw() + labs(title="Processus simulés", subtitle= sub) + scale_color_manual(values=cols)
  
  
  
  ggsave(paste0(path_res, 'processus_simules.png'), plot=gg,width = 15, height = 10, units = "cm")
  
  print(gg)
}



graphique_data_estimee_reelle <- function(solution, simulation, path_resultat="./"){
  
  X_hat <- estimation_X(solution, simulation)$X_hat 
  classe_hat <- estimation_X(solution, simulation)$classes
  X <- simulation$data_X
  
  REFERENCE <- simulation$REFERENCE
  
  d <- solution$size$d
  
  nsignaux <- solution$size$nsignaux
  
  classe_estim <- classe_hat%>%
    as.data.frame %>%
    rename('classe'='.') %>%
    mutate(individu = colnames(X))
  
  f <- function(j, X_hat, X, REFERENCE, classe_estim){
    
    lab_y = "Température (°C) "
    lab_x='Temps'
    
    table_data <- rbind( X_hat[REFERENCE$DIM==j,] %>%
                           as.data.frame %>%
                           mutate(temps=seq(1,nsignaux), type="estimations"), 
                         X[REFERENCE$DIM ==j, ]%>%
                           as.data.frame %>%
                           mutate(temps=seq(1,nsignaux), type="reelles"))
    
    
    
    g2 <- merge(table_data %>%
                  gather(key="individu", value='donnees', -temps, -type), classe_estim ) %>%
      mutate(classe = as.factor(classe)) %>%
      ggplot(aes(x=temps, y=donnees, color=classe, linetype=type)) + geom_line(size=0.8) + facet_wrap(.~individu, ncol=3) +
      theme_bw()+ labs(x=lab_x, y=lab_y, title=paste("Données réelles et estimées, dimension:", j)) +   scale_color_manual(values=c('darkslategrey', 'chartreuse3')) 
    ggsave(paste0(path_resultat, "data_x_wrap_",j,".png"), plot=gg,width = 30, height = 15, units = "cm")
    
    print(g2)
    
  }
  
  
  sapply(seq(1,d), f, X_hat, X, REFERENCE, classe_estim)
}




estimation_X_melange <- function(solution, simulation){
  
  n <- solution$size$n
  nsignaux <- solution$size$nsignaux
  d <- solution$size$d
  K <- solution$size$nbClasse
  proba <- solution$parametre_var$t_ik
  
  matrice_melange <- matrix(NA, nrow=(nsignaux+1), ncol=n)
  
  for(i in 1:n){
    effet_k <- 0 
    for(k in 1:K){
      ss <-   proba[i,k]*solution$filtre_Kalman$mu_kt_B[[k]]
      
      effet_k <- ss + effet_k
    }
    matrice_melange[,i] <- effet_k
    
  }
  
  classes <- apply(solution$parametre_var$t_ik, 1, which.max)
  n <- solution$size$n
  d <- solution$size$d
  nsignaux <- solution$size$nsignaux

  X_hat <- sapply(seq(1,n), function(i,solution, classes) cbind(sapply(seq(1,nsignaux), function(t, i, classes) 
    simulation$data_facteur[[t]]%*%solution$coefficient[[classes[i]]] + 
      + matrice_melange[(t+1),i], i, classes, simplify = "array")), solution, classes)
  
  colnames(X_hat) <- colnames(simulation$data_X)
  
  
  return(list(X_hat=X_hat, classes=classes))
}


