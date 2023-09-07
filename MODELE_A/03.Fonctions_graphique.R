graphique_simulation <- function(BDD, d, path_essais = "./", sub =" "){
  
  X <- BDD$data_X
  REFERENCE <- BDD$REFERENCE
  
  
  XX <- data.frame(temps = NA, dim=NA, individu=NA, x=NA)

  for(j in 1:d){
    XX <- rbind(XX, X[which(REFERENCE$DIM==j), ] %>%
                  as.data.frame %>%
                  mutate(temps = seq(1, (dim(X)[1]/d)), dim=j) %>%
                  gather(key="individu", value="x", -temps, -dim)
    )
  }
  
  if(d==1){
    ggg <- XX[-1,] %>%
      ggplot(aes(x=temps, y=x, group=individu)) + geom_line(size=0.6, color='grey8') + theme_bw() +
      theme(legend.position = "none") + 
      labs(x="Time", y="Simulated data", title = "Simulated Dataset", subtitle= sub)
    
  }else{
    ggg <- XX[-1,] %>%
      ggplot(aes(x=temps, y=x, group=individu)) + geom_line(size=0.6, color='grey8') + theme_bw() +
      facet_grid(dim~., labeller = label_both)+
      theme(legend.position = "none") + labs(x="Time", y="Simulated data", title = "Simulated Dataset", subtitle= sub)
    
  }

  print(ggg)
  ggsave(paste0(path_essais,"simulation_data_X.png"), plot=ggg, width = 15, height = 8, units = "cm")
  
  ##### 
  # PROCESSUS
  
  K = length(BDD$data_processus)
  
  if(is.null(dim(BDD$data_processus[[1]]))){
    d <- 1
    nsignaux <- length(BDD$data_processus[[1]]) -1
  } else {
    d <- dim(BDD$data_processus[[1]])[2]
    nsignaux <- dim(BDD$data_processus[[1]])[1] -1
  }
  
  table_processus <- data.frame(simulation=NA, classe=NA, temps=NA, dim=NA)
  
  for(j in 1:d){
    for(k in 1:K){
      table_processus <- rbind(table_processus, data.frame(simulation = as.matrix(BDD$data_processus[[k]], 
                                                                                  nrow=nsignaux+1, ncol=d)[-1,j],
                                                           classe=k, temps = seq(1,nsignaux), dim=j ))
    }
  }
  
  gg <- table_processus[-1,] %>%
    mutate(classe = as.factor(classe)) %>%
    rename('Classes' = classe) %>%
    ggplot(aes(x=temps, y=simulation, color=Classes)) + geom_line(size=0.8) + facet_grid(dim~.) +
    theme_bw() + labs(title="Class center dynamics", subtitle= sub, y="Class centers") +
    scale_color_manual(values=c('deeppink', 'lightseagreen', '#CCCC00', 'black'))
  
  
  
 ggsave(paste0(path_essais, 'processus_simules.png'), plot=gg, width = 15, height = 8, units = "cm")
  
  print(gg)
  
}


resultat_algorithme_VEM <- function(simulation_VEM, solution, chemin_sauvegarde='./'){
  
  graphique_matrice_classes(solution, simulation_VEM, chemin_sauvegarde)
  
  graphique_proccesus(solution, simulation_VEM,chemin_sauvegarde)

  
  graphique_effet_facteur(solution, simulation_VEM, chemin_sauvegarde)
  
}

#### FONCTION AFFICHAGE DE RESULTATS  

## 1 : Heatmap des classes 

graphique_matrice_classes <- function(solution, simulation, path_resultat = "./"){
  
  vec_classes_estimees <- apply(solution$parametre_var$t_ik, 1, which.max)
  vec_classes_simulees <- simulation$classes
  
  
  K <- solution$size$nbClasse
  
  
  gg <- table(vec_classes_estimees, vec_classes_simulees) %>%
    as.data.frame %>%
    rename(effectifs = Freq)%>%
    arrange(effectifs) %>%
    ggplot(aes(x=vec_classes_estimees, y=vec_classes_simulees, fill=effectifs)) + 
    geom_tile(colour = "grey50") + scale_fill_gradient2(low = "white", high = "tomato3") +
    geom_text(aes(label=effectifs)) +
    labs(x="Estimated classes", y="Simulated classes", 
         title = "The confusion matrix of estimated 
and simulated classes") + 
    theme_bw()  +  theme(legend.position = "none") 
  
  ggsave(paste0(path_resultat, "matrice_classes.png"), plot=gg, width = 8, height = 8, units = "cm") 
  
  
  print(gg)
}


## 3: Tracer les processus : 

graphique_proccesus <- function(solution, simulation , path_resultat = "./"){
  
  K <- solution$size$nbClasse
  vec_classes_estimees <- apply(solution$parametre_var$t_ik, 1, which.max)
  vec_classes_simulees <- simulation$classes
  
  if(length(unique(apply(solution$parametre_var$t_ik, 1, function(x) which.max(x)))) == 1){
    table_corresp_classes <- create_table_corresp_classes_proba(solution$parametre_var$t_ik, simulation$classes)
  }else{
    table_corresp_classes <- create_table_corresp_classes(solution, simulation)
  }
  
  d <- solution$size$d
  nsignaux <- solution$size$nsignaux
  K <- solution$size$nbClasse
  
  b_simulation <- simulation$data_processus
  b_estimation <- solution$filtre_Kalman$mu_kt_B
  
  table_processus <- data.frame(simulation=NA, estimation=NA, classe=NA, temps=NA, dim=NA)
  for(j in 1:d){
    for(k in 1:K){
      indice <- as.numeric(as.character(table_corresp_classes$estimation[table_corresp_classes$simulation==k]))
      table_processus <- rbind(table_processus, 
                               data.frame(simulation = b_simulation[[k]][-1,j], 
                                          estimation = b_estimation[[indice]][-1,j] , 
                                          classe=k, temps = seq(1,nsignaux), dim=j ))   
    }
    
  }
  
  if(d ==1){
    gg <-  table_processus[-1,] %>%
      gather(key='type', value='centres', -temps, -classe, -dim) %>%
      ggplot(aes(x=temps, y=centres, color=type)) + geom_line(size=0.8)+ 
      facet_grid(classe~., scale='free', labeller = label_both)+
      labs(title = "Estimated and simulated class centers", y="Class centers", 
           x="Time", color="Type")+
      theme_bw()+ scale_color_manual(values=c('orangered', 'indianred4'))
    
  }else{
    gg <-  table_processus[-1,] %>%
      gather(key='type', value='centres', -temps, -classe, -dim) %>%
      ggplot(aes(x=temps, y=centres, color=type)) + geom_line(size=0.8)+ 
      facet_grid(classe~dim, scale='free', labeller = label_both)+  
      labs(title = "Estimated and simulated class centers", y="Class centers", 
           x="Time", color="Type")+
      theme_bw()+ scale_color_manual(values=c('orangered', 'indianred4'))
    
  }
    
  ggsave(paste0(path_resultat, "Process_b.png"), plot=gg, width = 15, height = 8, unit='cm')
  
  print(gg)
}

## 4: Tracer l'effet des facteurs exogènes 
graphique_effet_facteur<- function(solution, simulation, path_resultat  = "./"){
  
  table_effet <- data.frame(simulation=NA, estimation=NA,dim= NA , temps=NA)
  nsignaux <- solution$size$nsignaux
  d <-solution$size$d 
 
  for(j in 1:d){
    table_effet <- rbind(table_effet, data.frame(simulation = sapply(seq(1,nsignaux), function(t,j) (simulation$data_facteur[[t]]%*%simulation$coefficient)[j,], j), 
                                                 estimation =sapply(seq(1,nsignaux), function(t,j) (simulation$data_facteur[[t]]%*%solution$coefficient)[j,], j) ,
                                                 dim=j, temps = seq(1,nsignaux) ))   
  }
  
  table_effet <- table_effet[-1,]
  
  gg <- table_effet %>%
    gather(key='type', value='valeurs', -temps, -dim) %>%
    ggplot(aes(x=temps, y=valeurs, color=type)) + geom_line(size=0.8)+
    labs(title = "Estimated and simulated effects of exogeneous factors", 
         x='Time', y='Exogeneous effects', color='Type')+ 
    theme_bw()+ scale_color_manual(values=c('orangered', 'indianred4'))
  
  ggsave(paste0(path_resultat, "Effets_facteurs_.png"), plot = gg, width = 15, height = 8, unit='cm')
  
  print(gg)
}


# Estimation 
estimation_X <- function(solution, simulation){
  
  classes <- apply(solution$parametre_var$t_ik, 1, which.max)
  n <- solution$size$n
  nsignaux <- solution$size$nsignaux
  d <- solution$size$d
  X_hat <- sapply(seq(1,n), function(i,solution, classes) cbind(sapply(seq(1,nsignaux), function(t, i, classes)
    simulation$data_facteur[[t]]%*%solution$coefficient +
      solution$filtre_Kalman$mu_kt_B[[classes[i]]][(t+1),], i, classes, simplify = "array")), solution, classes)
  colnames(X_hat) <- colnames(simulation$data_X)
  
  return(list(X_hat=X_hat, classes=classes))
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
  
  
  X_hat <- sapply(seq(1,n), function(i,solution) 
    cbind(sapply(seq(1,nsignaux),
                 function(t, i) simulation$data_facteur[[t]]%*%solution$coefficient + matrice_melange[(t+1),i], i, simplify = "array")), solution)

  colnames(X_hat) <- colnames(simulation$data_X)
  return(X_hat)
}







