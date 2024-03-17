##
## CONTINUITY CORRECTION
##

cc <- function(
    data, 
    value = 0.5,
    type = c("single", "all")[1]
){
  
  type <- match.arg(type)
  
  if(type == "all"){
    
    if(any(c(data$y0,data$y1,data$c0,data$c1) == 0)){
      
      data$y0 <- data$y0 + value
      data$y1 <- data$y1 + value
      data$c0 <- data$c0 + value
      data$c1 <- data$c1 + value
    }
  }
  
  if(type == "single"){
    
    correction = ((((data$y0 == 0)|(data$y1 == 0))|(data$c0 == 0))| (data$c1 == 0)) * value
    
    data$y0 <- correction + data$y0
    data$y1 <- correction + data$y1
    data$c0 <- correction + data$c0
    data$c1 <- correction + data$c1
    
  }
  
  return(data)
  
}
