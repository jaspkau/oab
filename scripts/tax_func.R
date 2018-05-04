
tax_assign = function(x, conf, rank){
  # x = data frame containig taxonomies and confidences
  # conf = column with confidences
  # rank = corresponding taxonomic ranks
  
  #make a list of confidences as some level
  list = x[[conf]]
  ##make a list of ranks at some level (correspoding to confidences)
  list2 = as.character(x[[rank]])
  ###see if the confidence is <95%, if so, change the rank to unknown
  list3 = ifelse(list < 0.90, "unidentified", list2) ###if list is <95 than assign unknown otherwise go with list2
  #replace the original rank with the new filtered rank based on threshold
  x[rank] = list3
  return(x)
}