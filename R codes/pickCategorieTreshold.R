  #@author Matej Mihelcic
  #@institution Rudjer Boskovic Institute, Zagreb, Croatia
  #@mail matmih1@gmail.com
  #@description function to select a category generality treshold

pickCategoryT<-function(x,p){
  ##1->small, 2->medium, 3->large
  categoryList=list();
  for(i in 1:length(x))
      categoryList[[i]]<-x[[i]][x[[i]][4]<=p,1:4];
  
  return(categoryList);
}