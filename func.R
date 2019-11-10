fold_enrichment = function(graph, ref_graph){
  #estimated_number_of_nodes = 114.53 #external data
  estimated_number_of_nodes = 94.3157894736842 #external data
  
  correct = 0
  incorrect = 0
  
  N = length(TF_names) * length(TF_names) #all possible combination
  predict = nrow(graph) # all predicted edges
  biogrid = nrow(ref_graph) #all biogrid edges
  
  for (i in seq(1:nrow(graph)) ){
    link = graph[i,]
    mark = link[1]
    target = link[2]
    
    for(r in seq(1:nrow(ref_graph))){
      row = ref_graph[r,]
      if( row[1] == mark){
        if(row[2] == target){
          correct = correct + 1
          break
        }
      }
    }
    incorrect = incorrect + 1
  }
  incorrect = incorrect - correct
  
  score_fold = (correct*N)/( predict * biogrid)
  ratio = correct/(correct + incorrect)
  penalty_score = score_fold * ( estimated_number_of_nodes / (1 + (predict - estimated_number_of_nodes)^2) )
  return( c(score_fold, ratio, penalty_score) )
}


plot_score_list = function( liste_scores ){
  
  dd  <-  as.data.frame(matrix(unlist(liste_scores), nrow=length(unlist(liste_scores[1]))))
  
  fold_enr = as.numeric(as.vector(dd[1,]))
  prec = as.numeric(as.vector(dd[2,]))
  num_iter = seq(1,length(prec), by=1)
  
  ddf = data.frame(num_iter , fold_enr, prec)
  
  data.m <- melt(ddf, id.vars='num_iter')
  
  ggplot(data.m, aes(num_iter, value)) + 
    geom_bar(aes(fill = variable), position = "dodge", stat="identity") + 
    theme(axis.text.x = element_text(angle = 90, hjust = 1))
  
}
