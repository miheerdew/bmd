stdize <-
  function(Mat){
    
    # Subtract mean of each row
    Mat = Mat - rowMeans(Mat)
    
    # Find sd of each row
    sds = sqrt(rowSums(Mat*Mat))
    
    # Return stdized
    return(Mat/sds)
  }

my_stdize <-
  function(Mat){
    
    # Subtract mean of each row
    Mat = Mat - rowMeans(Mat)
    
    # Find sd of each row
    sds = sqrt(rowSums(Mat*Mat))
    
    # Return stdized
    return(Mat/sds * sqrt(ncol(Mat) - 1))
  }