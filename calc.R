setwd("/Users/abhinav/Downloads/essentials/projects")
library(DescTools)
############################################################
AdditionSubtraction1 <- function() {
  result <- 0
  cat("Enter numbers by pressing 'ENTER' once, finish by pressing twice\n")
  raw_data <- scan(what = character(), quiet = TRUE)
  
  if (length(raw_data) <= 0) {
    cat("No numbers entered\n")
    return()
  }
  
  vect <- suppressWarnings(as.numeric(raw_data))
  if (any(is.na(vect))) {
    cat("Removed invalid non-numeric inputs\n")
    vect <- vect[!is.na(vect)]
  }
  
  if (length(vect) > 0) {
    for (i in vect){
      result <- result + i
    }
    cat("Result =", result, "\n")
  } else {
    cat("No valid numbers to calculate\n")
  }
}

Multiplication2 <- function() {
  result <- 1
  cat("Enter numbers by pressing 'ENTER' once, finish by pressing twice\n")
  raw_data <- scan(what = character(), quiet = TRUE)
  
  if (length(raw_data) <= 0) {
    cat("No numbers entered\n")
    return()
  }
  
  vect <- suppressWarnings(as.numeric(raw_data))
  if (any(is.na(vect))) {
    cat("Removed invalid non-numeric inputs\n")
    vect <- vect[!is.na(vect)]
  }
  
  if (length(vect) > 0) {
    for (i in vect){
      result <- result * i
    }
    cat("Result =", result, "\n")
  } else {
    cat("No valid numbers to calculate\n")
  }
}

Divide2Nos3 <- function() {
  cat("Enter 2 numbers by pressing 'ENTER' once, finish by pressing twice\n")
  raw_data <- scan(what = character(), n = 2, quiet = TRUE)
  
  if (length(raw_data) != 2) {
    cat("Please enter exactly 2 numbers\n")
    return()
  }
  
  vect <- suppressWarnings(as.numeric(raw_data))
  
  if (any(is.na(vect))) {
    cat("Invalid non-numeric inputs detected\n")
    return()
  }
  
  if (vect[2] == 0) {
    cat("Cannot divide by 0\n")
    return()
  }
  
  result <- vect[1] / vect[2]
  cat("Result =", result, "\n")
}

Exponent4 <- function() {
  cat("Enter a base and its power")
  base <- scan(what = character(), n = 1, quiet = TRUE)
  power <- scan(what = character(), n = 1, quiet = TRUE)
  
  base <- suppressWarnings(as.numeric(base))
  power <- suppressWarnings(as.numeric(power))
  
  if (any(is.na(c(base, power)))) {
    cat("Invalid non-numeric inputs detected\n")
    return()
  }
  
  if (base == 0 & power == 0) {
    cat("Not Defined.\n")
    return()
  }
  
  if (base == 0 & power < 0) {
    cat("Cannot divide by 0.\n")
    return()
  }
  
  result <- base ^ power
  cat("Result =", result, "\n")
}

Factorial5 <- function() {
  result <- 1
  cat("Enter a positive number")
  num <- scan(what = character(), n = 1, quiet = TRUE)
  
  num <- suppressWarnings(as.numeric(num))
  
  if (any(is.na(num))) {
    cat("Invalid non-numeric input detected\n")
    return()
  }
  
  if (num < 0) {
    cat("Cannot calculate factorial for -ve numbers\n")
  } else if (num == 0) {
    result <- 1
    cat("Result =", result, "\n")
  } else {
    for (i in 1:num) {
      result <- result * i
    }
    cat("Result =", result, "\n")
  }
}

fact <- function(x) {
  result <- 1
  for(i in 1:x){
    result <- result * i
  }
  return(result)
}

Permutation6 <- function() {
  cat("Enter a positive N and R")
  N <- scan(what = character(), n = 1, quiet = TRUE)
  R <- scan(what = character(), n = 1, quiet = TRUE)
  
  N <- suppressWarnings(as.integer(N))
  R <- suppressWarnings(as.integer(R))
  
  if (any(is.na(c(N, R)))) {
    cat("Invalid non-numeric input detected\n")
    return()
  }
  
  if ((N < 0) || (R < 0)) {
    cat("Enter valid N and R\n")
    return()
  }
  
  
  if (N == 0 && R == 0) {
    result <- 1
    cat("Result =", result, "\n")
  } else if (N < R) {
    result <- 0
    cat("Result =", result, "\n")
  } else {
    result <- fact(N) / fact(N-R)
    cat("Result =", result, "\n")
  }
  
}

Combination7 <- function() {
  cat("Enter a positive N and R")
  N <- scan(what = character(), n = 1, quiet = TRUE)
  R <- scan(what = character(), n = 1, quiet = TRUE)
  
  N <- suppressWarnings(as.integer(N))
  R <- suppressWarnings(as.integer(R))
  
  if (any(is.na(c(N, R)))) {
    cat("Invalid non-numeric input detected\n")
    return()
  }
  
  if ((N < 0) || (R < 0)) {
    cat("Enter valid N and R\n")
    return()
  }
  
  
  if (R == 0) {
    result <- 1
    cat("Result =", result, "\n")
  } else if (N < R) {
    result <- 0
    cat("Result =", result, "\n")
  } else if (N == R) {
    result <- 1
    cat("Result =", result, "\n")
  } else {
    result <- fact(N) / (fact(N-R) * fact(R))
    cat("Result =", result, "\n")
  }
  
}

AP8 <- function() {
  cat("Enter first term, common diff, and no. of terms\n")
  a <- scan(what = character(), n = 1, quiet = TRUE)
  d <- scan(what = character(), n = 1, quiet = TRUE)
  n <- scan(what = character(), n = 1, quiet = TRUE)
  
  a <- suppressWarnings(as.numeric(a))
  d <- suppressWarnings(as.numeric(d))
  n <- suppressWarnings(as.integer(n))
  
  if (any(is.na(c(a, d, n)))) {
    cat("Invalid non-numeric inputs detected\n")
    return()
  }
  
  if (n <= 0) {
    cat("Number of terms must be a +ve\n")
    return()
  }
  
  ap <- a + (0:(n-1)) * d
  cat("AP sequence is =\n", ap)
  
}

GP9 <- function() {
  cat("Enter first term, common ratio, and no. of terms\n")
  a <- scan(what = character(), n = 1, quiet = TRUE)
  r <- scan(what = character(), n = 1, quiet = TRUE)
  n <- scan(what = character(), n = 1, quiet = TRUE)
  
  a <- suppressWarnings(as.numeric(a))
  r <- suppressWarnings(as.numeric(r))
  n <- suppressWarnings(as.integer(n))
  
  if (any(is.na(c(a, r, n)))) {
    cat("Invalid non-numeric inputs detected.\n")
    return()
  }
  
  if (n <= 0) {
    cat("Number of terms must be a +ve.\n")
    return()
  }
  
  gp <- a * r^(0:(n-1))
  cat("GP sequence is =\n", gp)
  
}

HP10 <- function () {
  cat("Enter first term, common diff, and no. of terms\n")
  a <- scan(what = character(), n = 1, quiet = TRUE)
  d <- scan(what = character(), n = 1, quiet = TRUE)
  n <- scan(what = character(), n = 1, quiet = TRUE)
  
  a <- suppressWarnings(as.numeric(a))
  d <- suppressWarnings(as.numeric(d))
  n <- suppressWarnings(as.integer(n))
  
  if (any(is.na(c(a, d, n)))) {
    cat("Invalid non-numeric inputs detected\n")
    return()
  }
  
  if (n <= 0) {
    cat("Number of terms must be a +ve\n")
    return()
  }
  
  ap <- a + (0:(n-1)) * d
  if (any(ap == 0)) {
    cat("HP couldn't be formed (division by zero)\n")
  } else {
    hp <- 1 / ap
    cat("HP sequence is =\n", hp)
  }
}

SimpleInterest11 <- function () {
  options(scipen = 999)
  cat("Enter Principal Amount")
  Pri <- scan(what = character(), n = 1, quiet = TRUE)
  cat("Enter ROI (% per annum)")
  Rat <- scan(what = character(), n = 1, quiet = TRUE)
  cat("Enter Time (years)")
  Tim <- scan(what = character(), n = 1, quiet = TRUE)
  
  Pri <- suppressWarnings(as.numeric(Pri))
  Rat <- suppressWarnings(as.numeric(Rat))
  Tim <- suppressWarnings(as.numeric(Tim))
  
  if (any(is.na(c(Pri, Rat, Tim)))) {
    cat("Invalid non-numeric inputs detected.\n")
    return()
  }
  
  if (Pri < 0 || Rat < 0 || Tim < 0) {
    cat("Principal, Rate, and Time cannot be negative\n")
    return()
  }
  
  SI <- (Pri * Rat * Tim) / 100
  
  cat("\n--- Simple Interest Calculation ---\n")
  cat("Principal (P)   :", round(Pri, 2), "INR\n")
  cat("Rate (R)        :", round(Rat, 2), "%\n")
  cat("Time (T)        :", Tim, "years\n")
  cat("Simple Interest :", round(SI, 2), "INR\n")
  cat("Total Amount    :", round(Pri + SI, 2), "INR\n")
}

CompoundInterest12 <- function() {
  options(scipen = 999)
  cat("Enter Principal Amount")
  Pri <- scan(what = character(), n = 1, quiet = TRUE)
  cat("Enter ROI (% per annum)")
  Rat <- scan(what = character(), n = 1, quiet = TRUE)
  cat("Enter Time (years)")
  Tim <- scan(what = character(), n = 1, quiet = TRUE)
  cat("Enter compounding freq")
  Com <- scan(what = character(), n = 1, quiet = TRUE)
  
  Pri <- suppressWarnings(as.numeric(Pri))
  Rat <- suppressWarnings(as.numeric(Rat))
  Tim <- suppressWarnings(as.numeric(Tim))
  Com <- suppressWarnings(as.numeric(Com))
  
  if (any(is.na(c(Pri, Rat, Tim, Com)))) {
    cat("Invalid non-numeric inputs detected\n")
    return()
  }
  
  
  if (Pri < 0 || Rat < 0 || Tim < 0 || Com <= 0) {
    cat("Principal, Rate, Time and Compounding freq. cannot be negative\n")
    return()
  }
  
  r <- Rat / 100
  
  CI <- (Pri * ((1 + (r / Com))^(Com * Tim))) - Pri
  
  cat("\n--- Compond Interest Calculation ---\n")
  cat("Principal (P)   :", round(Pri, 2), "INR\n")
  cat("Rate (R)        :", round(Rat, 2), "%\n")
  cat("Time (T)        :", Tim, "years\n")
  cat("Comp. freq.     :", Com, "times\n")
  cat("Simple Interest :", round(CI, 2), "INR\n")
  cat("Total Amount    :", round(Pri + CI, 2), "INR\n")
}

QuadraticRoots13 <- function() {
  cat("Enter quadratic eq. coeffs.")
  a <- scan(what = character(), n = 1, quiet = TRUE)
  b <- scan(what = character(), n = 1, quiet = TRUE)
  c <- scan(what = character(), n = 1, quiet = TRUE)
  
  a <- suppressWarnings(as.numeric(a))
  b <- suppressWarnings(as.numeric(b))
  c <- suppressWarnings(as.numeric(c))
  
  if (any(is.na(c(a, b, c)))) {
    cat("Invalid non-numeric input detected\n")
    return()
  }
  
  if (a == 0) {
    cat("'a' cannot be 0\n")
    return()
  }
  
  d <- (b ^ 2) - (4*a*c)
  if (d >= 0) {
    root1 <- (-b + sqrt(d)) / (2*a)
    root2 <- (-b - sqrt(d)) / (2*a)
    result <- round(c(root1, root2), 4)
    cat("Result =\n")
    cat("\n",result[1],"\n", result[2])
  } else {
    real <- (-b) / (2*a)
    imag <- sqrt(-d) / (2*a)
    root1 <- paste(round(real,4),"+",round(imag, 4),"i")
    root2 <- paste(round(real,4),"-",round(imag, 4),"i")
    result <- c(root1, root2)
    cat("Result =\n")
    cat("\n",result[1],"\n", result[2])
  }
}

HCF14 <- function() {
  cat("Enter 2 positive natural nos.")
  
  num1 <- scan(what = character(), n = 1, quiet = TRUE)
  num2 <- scan(what = character(), n = 1, quiet = TRUE)
  
  num1 <- suppressWarnings(as.integer(num1))
  num2 <- suppressWarnings(as.integer(num2))
  
  if (any(is.na(c(num1, num2)))) {
    cat("Invalid non-numeric input detected\n")
    return()
  }
  
  if ((num1 < 0) || (num2 < 0)) {
    cat("Enter valid numbers\n")
    return()
  }
  
  if (num1 == 0 && num2 == 0){
    result <- 0
    cat("Result =", result, "\n")
  } else if (num1 == 0 || num2 == 0) {
    result <- max(c(num1, num2))
    cat("Result =", result, "\n")
  } else if (num1 == num2) {
    result <- num1
    cat("Result =", result, "\n")
  } else {
    while (num2 != 0) {
      temp <- num2
      num2 <- num1 %% num2
      num1 <- temp
    }
    result <- num1
    cat("Result =", result, "\n")
  }
  
}

gr_cm_di <- function(a, b) {
  while (b != 0) {
    temp <- b
    b <- a %% b
    a <- temp
  }
  result <- a
  return(result)
}

LCM15 <- function() {
  cat("Enter 2 positive natural nos.")
  
  num1 <- scan(what = character(), n = 1, quiet = TRUE)
  num2 <- scan(what = character(), n = 1, quiet = TRUE)
  
  num1 <- suppressWarnings(as.integer(num1))
  num2 <- suppressWarnings(as.integer(num2))
  
  if (any(is.na(c(num1, num2)))) {
    cat("Invalid non-numeric input detected\n")
    return()
  }
  
  if ((num1 < 0) || (num2 < 0)) {
    cat("Enter valid numbers\n")
    return()
  }
  
  if (num1 == 0 || num2 == 0){
    result <- 0
    cat("Result =", result, "\n")
  } else if (num1 == num2) {
    result <- num1
    cat("Result =", result, "\n")
  } else {
    result <- abs(num1 * num2) / (gr_cm_di(num1, num2))
    cat("Result =", result, "\n")
  }
  
}

MatrixAdd16 <- function() {
  cat("Enter dimension of 2 matrices")
  
  rows <- scan(what = character(), n = 1, quiet = TRUE)
  cols <- scan(what = character(), n = 1, quiet = TRUE)
  
  rows <- suppressWarnings(as.integer(rows))
  cols <- suppressWarnings(as.integer(cols))
  
  if (any(is.na(c(rows, cols)))) {
    cat("Invalid non-numeric input detected\n")
    return()
  }
  
  if ((rows < 0) || (cols < 0)) {
    cat("Enter valid dimensions\n")
    return()
  }
  
  cat("Enter 1st matrix elements row-wise")
  vect1 <- scan(what = character(), n = (rows*cols), quiet = TRUE)
  cat("Enter 2nd matrix elements row-wise")
  vect2 <- scan(what = character(), n = (rows*cols), quiet = TRUE)
  vect1 <- suppressWarnings(as.numeric(vect1))
  vect2 <- suppressWarnings(as.numeric(vect2))
  if (any(is.na(c(vect1, vect2)))) {
    cat("Invalid non-numeric input detected\n")
    return()
  }
  
  
  mat1 <- matrix(data = vect1, nrow = rows, ncol = cols, byrow = TRUE)
  mat2 <- matrix(data = vect2, nrow = rows, ncol = cols, byrow = TRUE)
  
  result <- mat1 + mat2
  cat("Result =\n") 
  print(result)
}

MatrixSubtract17 <- function() {
  cat("Enter dimension of 2 matrices")
  
  rows <- scan(what = character(), n = 1, quiet = TRUE)
  cols <- scan(what = character(), n = 1, quiet = TRUE)
  
  rows <- suppressWarnings(as.integer(rows))
  cols <- suppressWarnings(as.integer(cols))
  
  if (any(is.na(c(rows, cols)))) {
    cat("Invalid non-numeric input detected\n")
    return()
  }
  
  if ((rows < 0) || (cols < 0)) {
    cat("Enter valid dimensions\n")
    return()
  }
  
  cat("Enter 1st matrix elements row-wise")
  vect1 <- scan(what = character(), n = (rows*cols), quiet = TRUE)
  cat("Enter 2nd matrix elements row-wise")
  vect2 <- scan(what = character(), n = (rows*cols), quiet = TRUE)
  vect1 <- suppressWarnings(as.numeric(vect1))
  vect2 <- suppressWarnings(as.numeric(vect2))
  if (any(is.na(c(vect1, vect2)))) {
    cat("Invalid non-numeric input detected.\n")
    return()
  }
  
  
  mat1 <- matrix(data = vect1, nrow = rows, ncol = cols, byrow = TRUE)
  mat2 <- matrix(data = vect2, nrow = rows, ncol = cols, byrow = TRUE)
  
  result <- mat1 - mat2
  cat("Result =\n") 
  print(result)
}

MatrixTranspose18 <- function() {
  cat("Enter dimension of a matrix")
  
  rows <- scan(what = character(), n = 1, quiet = TRUE)
  cols <- scan(what = character(), n = 1, quiet = TRUE)
  
  rows <- suppressWarnings(as.integer(rows))
  cols <- suppressWarnings(as.integer(cols))
  
  if (any(is.na(c(rows, cols)))) {
    cat("Invalid non-numeric input detected\n")
    return()
  }
  
  if ((rows < 0) || (cols < 0)) {
    cat("Enter valid dimensions\n")
    return()
  }
  
  cat("Enter matrix elements row-wise")
  vect1 <- scan(what = character(), n = (rows*cols), quiet = TRUE)
  vect1 <- suppressWarnings(as.numeric(vect1))
  if (any(is.na(vect1))) {
    cat("Invalid non-numeric input detected\n")
    return()
  }
  mat1 <- matrix(data = vect1, nrow = rows, ncol = cols, byrow = TRUE)
  result <- t(mat1) 
  cat("Result =\n") 
  print(result)
}

MatrixMulti19 <- function() {
  cat("Enter dimension of 1st matrix")
  rows1 <- scan(what = character(), n = 1, quiet = TRUE)
  cols1 <- scan(what = character(), n = 1, quiet = TRUE)
  cat("Enter columns in 2nd matrix")
  cols2 <- scan(what = character(), n = 1, quiet = TRUE)
  
  rows1 <- suppressWarnings(as.integer(rows1))
  cols1 <- suppressWarnings(as.integer(cols1))
  cols2 <- suppressWarnings(as.integer(cols2))
  
  if (any(is.na(c(rows1, cols1, cols2)))) {
    cat("Invalid non-numeric input detected\n")
    return()
  }
  
  if ((rows1 < 0) || (cols1 < 0) || (cols2 < 0)) {
    cat("Enter valid dimensions\n")
    return()
  }
  
  cat("Enter 1st matrix elements row-wise")
  vect1 <- scan(what = character(), n = (rows1*cols1), quiet = TRUE)
  cat("Enter 2nd matrix elements row-wise")
  vect2 <- scan(what = character(), n = (cols1*cols2), quiet = TRUE)
  vect1 <- suppressWarnings(as.numeric(vect1))
  vect2 <- suppressWarnings(as.numeric(vect2))
  if (any(is.na(c(vect1, vect2)))) {
    cat("Invalid non-numeric input detected.\n")
    return()
  }
  
  
  mat1 <- matrix(data = vect1, nrow = rows1, ncol = cols1, byrow = TRUE)
  mat2 <- matrix(data = vect2, nrow = cols1, ncol = cols2, byrow = TRUE)
  
  result <- mat1 %*% mat2
  cat("Result =\n") 
  print(result)
  
}

MatrixInv20 <- function() {
  cat("Enter dimension of a sqaure matrix")
  
  dimn <- scan(what = character(), n = 1, quiet = TRUE)
  
  dimn <- suppressWarnings(as.integer(dimn))
  
  
  if (any(is.na(dimn))) {
    cat("Invalid non-numeric input detected\n")
    return()
  }
  
  if (dimn < 0) {
    cat("Enter valid dimensions\n")
    return()
  }
  
  cat("Enter matrix elements row-wise")
  vect1 <- scan(what = character(), n = (dimn*dimn), quiet = TRUE)
  vect1 <- suppressWarnings(as.numeric(vect1))
  if (any(is.na(vect1))) {
    cat("Invalid non-numeric input detected\n")
    return()
  }
  mat1 <- matrix(data = vect1, nrow = dimn, ncol = dimn, byrow = TRUE)
  
  if (det(mat1) == 0) {
    cat("Inverse does not exist")
    return()
  } else {
    result <- solve(mat1)
    cat("Result =\n") 
    print(result)
  }
}

Mathematics <- function() {
  repeat{
    cat("\n-----Mathematics menu-----\n")
    cat("1> Addition & Subtraction\n")
    cat("2> Multiplication\n")
    cat("3> Division\n")
    cat("4> Exponent\n")
    cat("5> Factorial\n")
    cat("6> Permutation\n")
    cat("7> Combination\n")
    cat("8> Create an AP\n")
    cat("9> Create a GP\n")
    cat("10> Create an HP\n")
    cat("11> Simple Interest\n")
    cat("12> Compound Interest\n")
    cat("13> Quadratic Equation Roots\n")
    cat("14> HCF\n")
    cat("15> LCM\n")
    cat("16> Matrix Addition\n")
    cat("17> Matrix Subtraction\n")
    cat("18> Matrix Transpose\n")
    cat("19> Matrix Multiplication\n")
    cat("20> Matrix Inverse\n")
    cat("21> Exit\n")
    input <- readline("Enter an integer: ")
    choice <- suppressWarnings(as.integer(input))
    
    if (is.na(choice)) {
      cat("\nInvalid choice, please enter 1-20.\n")
      next
    }
    
    switch (as.character(choice),
            "1" = AdditionSubtraction1(),
            "2" = Multiplication2(),
            "3" = Divide2Nos3(),
            "4" = Exponent4(),
            "5" = Factorial5(),
            "6" = Permutation6(),
            "7" = Combination7(),
            "8" = AP8(),
            "9" = GP9(),
            "10" = HP10(),
            "11" = SimpleInterest11(),
            "12" = CompoundInterest12(),
            "13" = QuadraticRoots13(),
            "14" = HCF14(),
            "15" = LCM15(),
            "16" = MatrixAdd16(),
            "17" = MatrixSubtract17(),
            "18" = MatrixTranspose18(),
            "19" = MatrixMulti19(),
            "20" = MatrixInv20(),
            "21" = { cat("\nExited Mathematics Section\n"); break },
            { cat("\nInvalid choice, please enter 1-20.\n") }
    )
  }
   
}
############################################################
Raw1 <- function() {
  cat("Enter each observation by pressing 'ENTER' once, finish by pressing twice\n")
  raw_data <- scan(what = character(), quiet = TRUE)
  
  if (length(raw_data) <= 0) {
    cat("No numbers entered\n")
    return()
  }
  
  vect <- suppressWarnings(as.numeric(raw_data))
  if (any(is.na(vect))) {
    cat("Removed invalid non-numeric inputs\n")
    vect <- vect[!is.na(vect)]
  }
  
  if (length(vect) > 0) {
    sum_of_obs <- 0
    for (i in vect) {
      sum_of_obs <- sum_of_obs + i
    }
    num_of_obs <- length(vect)
    AM <- sum_of_obs / num_of_obs
    
    prod_of_obs <- 1
    if (any(vect <= 0)) {
      GM <- "Can't compute GM"
    } else {
      for (i in vect){
        prod_of_obs <- prod_of_obs * i
      }
      GM <- prod_of_obs ^ (1 / num_of_obs)
    }
    
    sum_of_reci <- 0
    if (any(vect == 0)) {
      HM <- "Can't compute HM"
    } else {
      for (i in vect){
        sum_of_reci <- sum_of_reci + (1/i)
      }
      HM <- (sum_of_reci / num_of_obs) ^ (-1)
    }
    
    sorted_vect <- sort(vect)
    if (num_of_obs %% 2 == 1) { 
      Median <- sorted_vect[(num_of_obs + 1) / 2]
    } else { 
      mid1 <- sorted_vect[num_of_obs / 2]
      mid2 <- sorted_vect[(num_of_obs / 2) + 1]
      Median <- (mid1 + mid2) / 2
    }
    
    if (num_of_obs %% 2 == 0) {
      lower_half <- sorted_vect[1:(num_of_obs / 2)]
      upper_half <- sorted_vect[((num_of_obs / 2) + 1):num_of_obs]
    } else {
      lower_half <- sorted_vect[1:((num_of_obs-1) / 2)]
      upper_half <- sorted_vect[((num_of_obs + 3) / 2):num_of_obs]
    }
    
    len_lower <- length(lower_half)
    if (len_lower %% 2 == 1) {
      Q1 <- lower_half[(len_lower + 1) / 2]
    } else {
      Q1 <- (lower_half[len_lower/2] + lower_half[(len_lower/2) + 1]) / 2
    }
    
    len_upper <- length(upper_half)
    if (len_upper %% 2 == 1) {
      Q3 <- upper_half[(len_upper + 1) / 2]
    } else {
      Q3 <- (upper_half[len_upper/2] + upper_half[(len_upper/2) + 1]) / 2
    }
    

    freq_table <- table(vect)
    max_freq <- max(freq_table)
    mode_values <- as.numeric(names(freq_table)[freq_table == max_freq])
    if (length(mode_values) == 1) {
      Mode <- mode_values
    } else {
      Mode <- mode_values
    }
    
    
    Range <- max(vect) - min(vect)
    IQR <- Q3 - Q1
    Mean_dev_Mean <- sum(abs(vect - AM)) / num_of_obs
    Mean_dev_Median <- sum(abs(vect - Median)) / num_of_obs
    Variance <- sum((vect - AM)^2) / num_of_obs
    Stan_Dev <- Variance ^ (.5)
    CV <- AM / Stan_Dev 
    mu3 <- sum((vect - AM)^3) / num_of_obs
    mu4 <- sum((vect - AM)^4) / num_of_obs
    if (mu3 < 0) {
      BETA1 <- -1 * (mu3^2) / (Variance^3)
    } else if (mu3 == 0){
      BETA1 <- 0
    } else {
      BETA1 <- (mu3^2) / (Variance^3)
    }
    BETA2 <- mu4 / (Variance^2)
    if (mu3 < 0) {
      GAMMA1 <- -sqrt(-1*BETA1)
    } else if (mu3 == 0) {
      GAMMA1 <- 0
    } else {
      GAMMA1 <- sqrt(BETA1)
    }
    GAMMA2 <- BETA2 - 3
    
    cat(
      "\nDescriptive Statistics\n",
      "----------------------\n",
      "Count (n): ", num_of_obs, "\n",
      "Min: ", min(vect), "\n",
      "Max: ", max(vect), "\n",
      "Arithmetic Mean (AM): ", AM, "\n",
      "Geometric Mean (GM): ", if (is.character(GM)) GM else GM, "\n",
      "Harmonic Mean (HM): ", if (is.character(HM)) HM else HM, "\n",
      "Q1: ", Q1, "\n",
      "Median: ", Median, "\n",
      "Q3: ", Q3, "\n",
      "Range: ", Range, "\n",
      "IQR: ", IQR, "\n",
      "Mode(s): ", Mode, "\n",
      "Mean Deviation about Mean: ", Mean_dev_Mean, "\n",
      "Mean Deviation about Median: ", Mean_dev_Median, "\n",
      "Variance: ", Variance, "\n",
      "Standard Deviation: ", Stan_Dev, "\n",
      "Coefficient of Variation: ", CV, "\n",
      "Beta1 (skewness): ", BETA1, "\n",
      "Gamma1 (skewness): ", GAMMA1, "\n",
      "Beta2 (kurtosis): ", BETA2, "\n",
      "Gamma2 (kurtosis): ", GAMMA2, "\n"
    )
  } else {
    cat("No valid numbers to calculate\n")
    return()
  }
}

Frequency2 <- function() {
  cat("Enter xi (distinct values) separated by spaces, then press ENTER\n")
  xi_input <- scan(what = character(), quiet = TRUE)
  if (length(xi_input) <= 0) {
    cat("No xi entered\n")
    return()
  }
  xi <- suppressWarnings(as.numeric(xi_input))
  if (any(is.na(xi))) {
    cat("Removed invalid non-numeric xi\n")
    xi <- xi[!is.na(xi)]
  }
  if (length(xi) <= 0) {
    cat("No valid xi to proceed\n")
    return()
  }
  
  cat("Enter frequencies fi for each xi, then press ENTER\n")
  fi_input <- scan(what = character(), quiet = TRUE)
  if (length(fi_input) <= 0) {
    cat("No frequencies entered\n")
    return()
  }
  fi <- suppressWarnings(as.numeric(fi_input))
  if (any(is.na(fi))) {
    cat("Removed invalid non-numeric fi\n")
    fi <- fi[!is.na(fi)]
  }
  
  
  if (length(fi) != length(xi)) {
    cat("Length mismatch, length(fi) must be equal to length(xi)\n")
    return()
  }
  if (any(fi < 0)) {
    cat("Invalid frequencies: fi must be non-negative\n")
    return()
  }
  if (sum(fi) == 0) {
    cat("All frequencies are zero\n")
    return()
  }
  
  keep <- fi > 0
  xi <- xi[keep]
  fi <- fi[keep]
  
  n <- sum(fi)
  wsum <- sum(fi * xi)
  AM <- wsum / n
  
  if (any(xi <= 0)) {
    GM <- "Can't compute GM"
  } else {
    GM <- 10^((sum(log10(xi) * fi)) / sum(fi))
  }
  
  
  if (any(xi == 0)) {
    HM <- "Can't compute HM"
  } else if (any(xi < 0)) {
    HM <- "Can't compute HM"
  } else {
    HM <- (sum((1/xi) * fi) / sum(fi)) ^ (-1)
  }
  
  vect <- rep(xi, fi)
  vect <- sort(vect)
  num_of_obs <- length(vect)
  
  if (num_of_obs %% 2 == 1) {
    Median <- vect[(num_of_obs + 1) / 2]
  } else {
    Median <- mean(vect[(num_of_obs/2):(num_of_obs/2 + 1)])
  }
  
  if (num_of_obs %% 2 == 0) {
    lower_half <- vect[1:(num_of_obs / 2)]
    upper_half <- vect[(num_of_obs / 2 + 1):num_of_obs]
  } else {
    lower_half <- vect[1:((num_of_obs - 1) / 2)]
    upper_half <- vect[((num_of_obs + 3) / 2):num_of_obs]
  }
  
  len_lower <- length(lower_half)
  if (len_lower %% 2 == 1) {
    Q1 <- lower_half[(len_lower + 1) / 2]
  } else {
    Q1 <- mean(lower_half[(len_lower/2):(len_lower/2 + 1)])
  }
  
  len_upper <- length(upper_half)
  if (len_upper %% 2 == 1) {
    Q3 <- upper_half[(len_upper + 1) / 2]
  } else {
    Q3 <- mean(upper_half[(len_upper/2):(len_upper/2 + 1)])
  }
  
  max_freq <- max(fi)
  Mode <- xi[fi == max_freq]
  
  Range <- max(xi) - min(xi)
  IQR <- Q3 - Q1
  Mean_dev_Mean <- sum(fi * abs(xi - AM)) / n
  Mean_dev_Median <- mean(abs(vect - Median))
  Variance <- sum(fi * (xi - AM)^2) / n
  Stan_Dev <- sqrt(Variance)
  CV <- AM / Stan_Dev
  mu3 <- sum(fi * (xi - AM)^3) / n
  mu4 <- sum(fi * (xi - AM)^4) / n
  
  if (mu3 < 0) {
    BETA1 <- -1 * (mu3^2) / (Variance^3)
  } else if (mu3 == 0) {
    BETA1 <- 0
  } else {
    BETA1 <- (mu3^2) / (Variance^3)
  }
  BETA2 <- mu4 / (Variance^2)
  if (mu3 < 0) {
    GAMMA1 <- -sqrt(-1 * BETA1)
  } else if (mu3 == 0) {
    GAMMA1 <- 0
  } else {
    GAMMA1 <- sqrt(BETA1)
  }
  GAMMA2 <- BETA2 - 3
  
  cat(
    "\nDescriptive Statistics\n",
    "----------------------\n",
    "Count (n): ", num_of_obs, "\n",
    "Min: ", min(vect), "\n",
    "Max: ", max(vect), "\n",
    "Arithmetic Mean (AM): ", AM, "\n",
    "Geometric Mean (GM): ", if (is.character(GM)) GM else GM, "\n",
    "Harmonic Mean (HM): ", if (is.character(HM)) HM else HM, "\n",
    "Q1: ", Q1, "\n",
    "Median: ", Median, "\n",
    "Q3: ", Q3, "\n",
    "Range: ", Range, "\n",
    "IQR: ", IQR, "\n",
    "Mode(s): ", Mode, "\n",
    "Mean Deviation about Mean: ", Mean_dev_Mean, "\n",
    "Mean Deviation about Median: ", Mean_dev_Median, "\n",
    "Variance: ", Variance, "\n",
    "Standard Deviation: ", Stan_Dev, "\n",
    "Coefficient of Variation: ", CV, "\n",
    "Beta1 (skewness): ", BETA1, "\n",
    "Gamma1 (skewness): ", GAMMA1, "\n",
    "Beta2 (kurtosis): ", BETA2, "\n",
    "Gamma2 (kurtosis): ", GAMMA2, "\n"
  )
}

Bivariate3 <- function() {
  cat("Enter each Xi (Independent) by pressing 'ENTER' once, finish by pressing twice\n")
  Xi <- scan(what = character(), quiet = TRUE)
  cat("Enter each Yi (Dependent) by pressing 'ENTER' once, finish by pressing twice\n")
  Yi <- scan(what = character(), quiet = TRUE)
  
  if (length(Xi) <= 0 || length(Yi) <= 0) {
    cat("No numbers entered\n")
    return()
  }
  
  Xi <- suppressWarnings(as.numeric(Xi))
  if (any(is.na(Xi))) {
    cat("Removed invalid non-numeric inputs\n")
    Xi <- Xi[!is.na(Xi)]
  }
  
  Yi <- suppressWarnings(as.numeric(Yi))
  if (any(is.na(Yi))) {
    cat("Removed invalid non-numeric inputs\n")
    Yi <- Yi[!is.na(Yi)]
  }
  if (length(Xi) != length(Yi)) {
    cat("Length mismatch: length(Xi) must equal length(Yi)\n")
    return()
  }
  
  num_of_obs <- length(Xi)
  sum_of_Xi <- 0
  for (i in Xi) {
    sum_of_Xi <- sum_of_Xi + i
  }
  sum_of_Xisq <- 0 
  for (i in Xi) {
    sum_of_Xisq <- sum_of_Xisq + i^2
  }
  sum_of_Yi <- 0
  for (i in Yi) {
    sum_of_Yi <- sum_of_Yi + i
  }
  sum_of_Yisq <- 0 
  for (i in Yi) {
    sum_of_Xisq <- sum_of_Yisq + i^2
  }
  sum_of_XiYi <- sum(Xi * Yi)
  
  X_bar <- sum_of_Xi / num_of_obs
  Y_bar <- sum_of_Yi / num_of_obs
  VarX <- sum((Xi - X_bar)^2) / num_of_obs
  SDX <- sqrt(VarX)
  VarY <- sum((Yi - Y_bar)^2) / num_of_obs
  SDY <- sqrt(VarY)
  COV <- (sum_of_XiYi / num_of_obs) - (X_bar * Y_bar)
  CORR <- COV / sqrt(VarX * VarY)
  RegCoeff <- CORR*(SDY/SDX)
  RegEqn <- paste("Y =",round(Y_bar,3) - round(RegCoeff*(X_bar),3), "+", round(RegCoeff,3),"X")
  Intercept <- Y_bar - RegCoeff*(X_bar) + RegCoeff*0
  PE <- round(((1 - CORR^2)/sqrt(num_of_obs)) * 0.6745, 4)
  if (abs(CORR) > 6*PE) {
    Significant <- paste("YES")
  } else {
    Significant <- paste("NO")
  }
  cat(
    "\nBivariate Summary\n",
    "-----------------\n",
    "n: ", num_of_obs, "\n",
    "Mean of X: ", X_bar, "\n",
    "Mean of Y: ", Y_bar, "\n",
    "Var(X): ", VarX, "\n",
    "SD(X): ", SDX, "\n",
    "Var(Y): ", VarY, "\n",
    "SD(Y): ", SDY, "\n",
    "Cov(X,Y): ", COV, "\n",
    "Corr(X,Y): ", CORR, "\n",
    "Regression (Y on X): ", RegEqn, "\n",
    "Slope (b): ", RegCoeff, "\n",
    "Intercept (a): ", Intercept, "\n",
    "Probable Error (PE) of r: ", PE, "\n",
    "Significant (|r| > 6*PE): ", Significant, "\n"
  )
  
  
}

DescStats <- function() {
  repeat {
    cat("\n-----Descriptive Statistics menu-----\n")
    cat("\n1] Raw Data")
    cat("\n2] Frequency Data")
    cat("\n3] Bivariate Data")
    cat("\n4] Exit\n")
    input <- readline("Enter your choice: ")
    choice <- suppressWarnings(as.integer(input))
    
    if (is.na(choice)) {
      cat("\nInvalid choice, please choose a valid option!!\n")
      next
    }
    
    switch(as.character(choice),
           "1" = Raw1(),
           "2" = Frequency2(),
           "3" = Bivariate3(),
           "4" = { cat("\nExited Desc.Stats. section\n"); break },
           { cat("\nInvalid choice, please choose a valid option!!\n") }
    )
  }
   
}
############################################################
BisectionMethod1 <- function() {
  cat("Choose a function to find root\n")
  cat("1) f(x) = x^3 - x - 2\n")
  cat("2) f(x) = cos(x) - x^3\n")
  cat("3) f(x) = x^2 - 2\n")
  cat("Enter 1, 2, or 3")
  
  choice <- scan(what = integer(), n = 1, quiet = TRUE)
  
  if (length(choice) == 0 || !(choice %in% c(1, 2, 3))) {
    cat("Invalid choice for functions\n")
    return()
  }
  
  if (choice == 1) {
    f <- function(x) x^3 - x - 2
    e <- 1e-8
  } else if (choice == 2) {
    f <- function(x) cos(x) - x^3
    e <- 1e-8
  } else if (choice == 3) {
    f <- function(x) x^2 - 2
    e <- 1e-8
  }
  
  repeat {
    cat("Enter 2 initial approximations")
    a <- scan(what = character(), n = 1, quiet = TRUE)
    b <- scan(what = character(), n = 1, quiet = TRUE)
    
    a <- suppressWarnings(as.numeric(a))
    b <- suppressWarnings(as.numeric(b))
    
    if (any(is.na(c(a, b)))) {
      cat("Invalid non-numeric inputs detected\n")
      return( )
    }
    
    if (f(a) * f(b) < 0.0) {
      break
    } else {
      cat("f(a) * f(b) >= 0, Please enter a valid interval\n")
    }
  }
  
  k <- 0
  repeat {
    m <- (a + b) / 2
    k <- k + 1
    
    cat(sprintf("%d | m = %.4f | a = %.4f | b = %.4f | f(m) = %.4e\n", k, m, a, b, f(m)))
    
    fm <- f(m)
    if (abs(fm) < e) {
      cat(sprintf("\nHence, root = %.4f\n", m))
      cat("Total iterations = ", k, "\n")
      curve(f, from = m - 2, to = m + 2, n = 1000, col = "blue", lwd = 2)
      abline(h = 0, col = "black")
      points(m, 0, col = "red", pch = 19)
      text(m, 0, label = round(m,4),  col = "red", pos = 3)
      break
    }
    
    if (f(a) * fm < 0.0) {
      b <- m
    } else {
      a <- m
    }
  }
}

SecantMethod2 <- function() {
  cat("Choose a function to find root\n")
  cat("1) f(x) = x^3 - x - 2\n")
  cat("2) f(x) = cos(x) - x^3\n")
  cat("3) f(x) = x^2 - 2\n")
  cat("Enter 1, 2, or 3")
  
  choice <- scan(what = integer(), n = 1, quiet = TRUE)
  
  if (length(choice) == 0 || !(choice %in% c(1, 2, 3))) {
    cat("Invalid choice. Exiting.\n")
    return( )
  }
  
  if (choice == 1) {
    f <- function(x) x^3 - x - 2
    e <- 1e-8
  } else if (choice == 2) {
    f <- function(x) cos(x) - x^3
    e <- 1e-8
  } else if (choice == 3) {
    f <- function(x) x^2 - 2
    e <- 1e-8
  }
  
  repeat {
    cat("Enter 2 initial approximations")
    x0 <- scan(what = character(), n = 1, quiet = TRUE)
    x1 <- scan(what = character(), n = 1, quiet = TRUE)
    
    x0 <- suppressWarnings(as.numeric(x0))
    x1 <- suppressWarnings(as.numeric(x1))
    
    if (any(is.na(c(x0, x1)))) {
      cat("Invalid non-numeric inputs detected\n")
      return()
    }
    
    fx0 <- f(x0)
    fx1 <- f(x1)
    if (abs(fx1 - fx0) > 1e-12) {
      break
    } else {
      cat("f(x1) - f(x0) is too small, choose different approximations\n")
    }
  }
  
  k <- 0
  
  repeat {
    fx0 <- f(x0)
    fx1 <- f(x1)
    
    x2 <- x1 - (fx1 * (x1 - x0) / (fx1 - fx0))
    k <- k + 1
    
    cat(sprintf("%d | x2 = %.4f | x0 = %.4f | x1 = %.4f\n", k, x2, x0, x1))
    fx2 <- f(x2)
    if (abs(fx2) < e) {
      cat(sprintf("\nHence, root = %.4f\n", x2))
      cat("Total iterations = ", k, "\n")
      curve(f, from = x2-2, to = x2+2, n = 1000, col = "blue", lwd = 2)
      abline(h = 0, col = "black")
      points(x2, 0, col = "red", pch = 19)
      text(x2, 0, label = round(x2,4),  col = "red", pos = 3)
      break
    }
    
    x0 <- x1
    x1 <- x2
  }
}

RegulaFalsi3 <- function() {
  cat("Choose a function to find root\n")
  cat("1) f(x) = x^3 - x - 2\n")
  cat("2) f(x) = cos(x) - x^3\n")
  cat("3) f(x) = x^2 - 2\n")
  cat("Enter 1, 2, or 3")
  
  choice <- scan(what = integer(), n = 1, quiet = TRUE)
  
  if (length(choice) == 0 || !(choice %in% c(1, 2, 3))) {
    cat("Invalid choice. Exiting.\n")
    return( )
  }
  
  if (choice == 1) {
    f <- function(x) x^3 - x - 2
    e <- 1e-8
  } else if (choice == 2) {
    f <- function(x) cos(x) - x^3
    e <- 1e-8
  } else if (choice == 3) {
    f <- function(x) x^2 - 2
    e <- 1e-8
  }
  
  repeat {
    cat("Enter 2 initial approximations: ")
    a <- scan(what = character(), n = 1, quiet = TRUE)
    b <- scan(what = character(), n = 1, quiet = TRUE)
    
    a <- suppressWarnings(as.numeric(a))
    b <- suppressWarnings(as.numeric(b))
    
    if (any(is.na(c(a, b)))) {
      cat("Invalid non-numeric inputs detected.\n")
      return()
    }
    
    if (f(a) * f(b) < 0.0) {
      break
    } else {
      cat("f(a) and f(b) must have opposite signs. Try again.\n")
    }
  }
  
  k <- 0
  
  repeat {
    fa <- f(a)
    fb <- f(b)
    c <- a - fa * (b - a) / (fb - fa)
    k <- k + 1
    
    cat(sprintf("%d | c = %.4f | a = %.4f | b = %.4f\n", k, c, a, b))
    
    fc <- f(c)
    if (abs(fc) < e) {
      cat(sprintf("\nHence, root = %.4f\n", c))
      cat("Total iterations = ", k, "\n")
      curve(f, from = c-2, to = c+2, n = 1000, col = "blue", lwd = 2)
      abline(h = 0, col = "black")
      points(c, 0, col = "red", pch = 19)
      text(c, 0, label = round(c,4),  col = "red", pos = 3)
      break
    }
    
    if (fa * fc < 0.0) {
      b <- c
    } else {
      a <- c
    }
  }
}

NewtonRaph4 <- function() {
  cat("Choose a function to find root\n")
  cat("1) f(x) = x^3 - x - 2\n")
  cat("2) f(x) = cos(x) - x^3\n")
  cat("3) f(x) = x^2 - 2\n")
  cat("Enter 1, 2, or 3")
  
  choice <- scan(what = integer(), n = 1, quiet = TRUE)
  
  if (length(choice) == 0 || !(choice %in% c(1, 2, 3))) {
    cat("Invalid choice. Exiting.\n")
    return( )
  }
  
  if (choice == 1) {
    f <- function(x) x^3 - x - 2
    fdash <- function(x) 3*x^2 - 1
    e <- 1e-8
  } else if (choice == 2) {
    f <- function(x) cos(x) - x^3
    fdash <- function(x) -sin(x) - 3*x^2
    e <- 1e-8
  } else if (choice == 3) {
    f <- function(x) x^2 - 2
    fdash <- function(x) 2*x
    e <- 1e-8
  }
  
  repeat {
    cat("Enter initial approximation: ")
    x0 <- scan(what = character(), n = 1, quiet = TRUE)
    x0 <- suppressWarnings(as.numeric(x0))
    
    if (is.na(x0)) {
      cat("Invalid non-numeric input detected.\n")
      return()
    }
    
    d0 <- fdash(x0)
    if (abs(d0) > 1e-12) {
      break
    } else {
      cat("Derivative is too small at this point. Choose a different value.\n")
    }
  }
  
  k <- 0
  x <- x0
  
  repeat {
    fx <- f(x)
    dfx <- fdash(x)
    
    x_next <- x - fx / dfx
    k <- k + 1
    
    cat(sprintf("%d | %.6f | %.6f | %.6f\n", k, x_next, fx, dfx))
    
    fx_next <- f(x_next)
    if (abs(fx_next) < e) {
      cat(sprintf("\nHence, root = %.6f\n", x_next))
      cat("Total iterations = ", k, "\n")
      curve(f, from = x_next-2, to = x_next+2, n = 1000, col = "blue", lwd = 2)
      abline(h = 0, col = "black")
      points(x_next, 0, col = "red", pch = 19)
      text(x_next, 0, label = round(x_next,4),  col = "red", pos = 3)
      break
    }
    
    if (abs(fdash(x_next)) < 1e-12) {
      cat("Derivative became too small. Method may not converge.\n")
      break
    }
    
    x <- x_next
  }
}

ForwardSub5 <- function() {
  cat("Enter order of square matrix\n")
  dimn_raw <- scan(what = character(), n = 1, quiet = TRUE)
  if (length(dimn_raw) != 1) {
    cat("Please enter exactly 1 number\n") 
    return()
  }
  dimn <- suppressWarnings(as.integer(dimn_raw))
  if (is.na(dimn) || dimn <= 0) {
    cat("Invalid order, must be a positive integer\n")
    return()
  }
  
  total <- dimn * dimn
  cat("Enter the coefficient matrix (row-wise)\n")
  vals_raw <- scan(what = character(), n = total, quiet = TRUE)
  if (length(vals_raw) != total) {
    cat("Invalid length")
    return()
  }
  vals <- suppressWarnings(as.numeric(vals_raw))
  if (any(is.na(vals))) {
    cat("Invalid numeric inputs\n")
    return()
  }
  L <- matrix(vals, nrow = dimn, ncol = dimn, byrow = TRUE)
  
  cat("Enter the elements of RHS\n")
  b_raw <- scan(what = character(), n = dimn, quiet = TRUE)
  if (length(b_raw) != dimn) {
    cat("Invalid length for b\n")
    return()
  }
  b <- suppressWarnings(as.numeric(b_raw))
  if (any(is.na(b))) {
    cat("Invalid numeric inputs in b\n")
    return()
  }
  
  if (any(diag(L) == 0)) {
    cat("Zero found on diagonal, forward substitution not possible\n")
    return()
  }
  
  cat("\nMatrix L\n")
  print(L)
  cat("Vector b\n") 
  print(as.matrix(b))
  result <- backsolve(L, b, upper.tri = FALSE)
  print(as.matrix(result))
}

BackwardSub6 <- function() { 
  cat("Enter order of square matrix\n")
  dimn_raw <- scan(what = character(), n = 1, quiet = TRUE)
  if (length(dimn_raw) != 1) {
    cat("Please enter exactly 1 number\n") 
    return()
  }
  dimn <- suppressWarnings(as.integer(dimn_raw))
  if (is.na(dimn) || dimn <= 0) {
    cat("Invalid order, must be a positive integer\n")
    return()
  }
  
  total <- dimn * dimn
  cat("Enter the coefficient matrix (row-wise)\n")
  vals_raw <- scan(what = character(), n = total, quiet = TRUE)
  if (length(vals_raw) != total) {
    cat("Invalid length")
    return()
  }
  vals <- suppressWarnings(as.numeric(vals_raw))
  if (any(is.na(vals))) {
    cat("Invalid numeric inputs\n")
    return()
  }
  L <- matrix(vals, nrow = dimn, ncol = dimn, byrow = TRUE)
  
  cat("Enter the elements of RHS\n")
  b_raw <- scan(what = character(), n = dimn, quiet = TRUE)
  if (length(b_raw) != dimn) {
    cat("Invalid length for b\n")
    return()
  }
  b <- suppressWarnings(as.numeric(b_raw))
  if (any(is.na(b))) {
    cat("Invalid numeric inputs in b\n")
    return()
  }
  
  if (any(diag(L) == 0)) {
    cat("Zero found on diagonal, backward substitution not possible\n")
    return()
  }
  
  cat("\nMatrix L\n")
  print(L)
  cat("Vector b\n") 
  print(as.matrix(b))
  result <- backsolve(L, b, upper.tri = TRUE)
  print(as.matrix(result))
}

Lagrange7 <- function() {
  cat("Enter each Xi by pressing 'ENTER' once, finish by pressing twice\n")
  Xi_raw <- scan(what = character(), quiet = TRUE)
  cat("Enter each FXi by pressing 'ENTER' once, finish by pressing twice\n")
  FXi_raw <- scan(what = character(), quiet = TRUE)
  
  if (length(Xi_raw) == 0 || length(FXi_raw) == 0) {
    cat("No numbers entered\n")
    return()
  }
  
  Xi <- suppressWarnings(as.numeric(Xi_raw))
  if (any(is.na(Xi))) {
    cat("Removed invalid non-numeric inputs in Xi\n")
    Xi <- Xi[!is.na(Xi)]
  }
  
  FXi <- suppressWarnings(as.numeric(FXi_raw))
  if (any(is.na(FXi))) {
    cat("Removed invalid non-numeric inputs in Yi\n")
    FXi <- FXi[!is.na(FXi)]
  }
  
  if (length(Xi) != length(FXi)) {
    cat("Length mismatch: length(Xi) must equal length(Yi)\n")
    return( )
  }
  
  n <- length(Xi)
  if (n <= 0) {
    cat("No valid data points after cleaning inputs.\n")
    return()
  }
  
  if (any(duplicated(Xi))) {
    dup_idx <- which(duplicated(Xi) | duplicated(Xi, fromLast = TRUE))
    cat("Error: duplicate Xi values detected, causes division by zero\n")
    return()
  }
  
  cat("\nEnter the value of x to interpolate at\n")
  x_raw <- scan(what = character(), nmax = 1, quiet = TRUE)
  if (length(x_raw) == 0) {
    cat("No x provided for interpolation\n")
    return()
  }
  x <- suppressWarnings(as.numeric(x_raw))
  if (is.na(x)) {
    cat("Invalid x\n")
    return()
  }
  
  result <- 0
  for (i in seq_len(n)) {
    term <- FXi[i]
    for (j in seq_len(n)) {
      if (j != i) {
        term <- term * ((x - Xi[j]) / (Xi[i] - Xi[j]))
      }
    }
    result <- result + term
  }
  
  cat("Value at x =", x, "is FXi =", result)
}

NewtonForward8 <- function() {
  cat("Under Construction")
  return()
}

NewtonBackward9 <- function() {
  cat("Under Construction")
  return()
}

# Trapezoidal10 <- function() {
#   cat("Choose a function to find area\n")
#   cat("1) f(x) = 1 / (1 + x)\n")
#   cat("2) f(x) = | cos(x) |\n")
#   cat("3) f(x) = x^2 + (x^4)/10\n")
#   cat("Enter 1, 2, or 3")
# 
#   choice <- scan(what = integer(), n = 1, quiet = TRUE)
# 
#   if (length(choice) == 0 || !(choice %in% c(1, 2, 3))) {
#     cat("Invalid choice. Exiting.\n")
#     return( )
#   }
# 
#   if (choice == 1) {
#     g <- function(x) 1 / (1 + x)
#   } else if (choice == 2) {
#     g <- function(x) abs(cos(x))
#   } else if (choice == 3) {
#     g <- function(x) x^2 + (x^4)/10
#   }
# 
#   #g <- function(x) 3*x - 4
#   cat("Enter the lower limit and upper limit \n")
#   ab <- scan(what = numeric(), nmax = 2, quiet = TRUE)
#   if (length(ab) < 2 || any(is.na(ab))) {
#     cat("Invalid input for a and b\n")
#     return()
#   }
#   a <- ab[1]
#   b <- ab[2]
# 
#   result <- ((g(a) + g(b)) / 2) * (b - a)
#   curve(g, from = min(c(a,b)) - 1, to = max(c(a,b)) +1, n = 1000, col = "blue", lwd = 2)
#   points(c(a, b), c(g(a), g(b)), pch = 19, col = "red", cex = 1.2)
#   abline(v = c(a, b), lty = 3, col = "red")
#   segments(a, g(a), b, g(b), lty = 3, col = "red")
#   polygon(x = c(a, a, b, b), y = c(0, g(a), g(b), 0),
#     col = adjustcolor("orange", alpha.f = 0.4), border = NA)
#   abline(h = 0, col = "black")
#   cat("Result = ", round(result, 4))
# }

# Simpson13rd11 <- function() {
#   cat("Choose a function to find area\n")
#   cat("1) f(x) = 1 / (1 + x)\n")
#   cat("2) f(x) = | cos(x) |\n")
#   cat("3) f(x) = x^2 + (x^4)/10\n")
#   cat("Enter 1, 2, or 3")
# 
#   choice <- scan(what = integer(), n = 1, quiet = TRUE)
# 
#   if (length(choice) == 0 || !(choice %in% c(1, 2, 3))) {
#     cat("Invalid choice. Exiting.\n")
#     return( )
#   }
# 
#   if (choice == 1) {
#     g <- function(x) 1 / (1 + x)
#   } else if (choice == 2) {
#     g <- function(x) abs(cos(x))
#   } else if (choice == 3) {
#     g <- function(x) x^2 + (x^4)/10
#   }
# 
#   #g <- function(x) 3*x^2 + 11
#   #g <- function(x) x^3-2*x^2+2
#   cat("Enter the lower limit and upper limit\n")
#   ab <- scan(what = numeric(), nmax = 2, quiet = TRUE)
#   if (length(ab) < 2 || any(is.na(ab))) {
#     cat("Invalid input for a and b\n")
#     return()
#   }
#   a <- ab[1]
#   b <- ab[2]
# 
#   h <- (b - a) / 2
#   x0 <- a
#   x1 <- a + h
#   x2 <- b
# 
#   result <- (h / 3) * (g(x0) + 4 * g(x1) + g(x2))
#   curve(g, from = min(c(a,b)) - 1, to = max(c(a,b)) +1, n = 1000, col = "blue", lwd = 2)
#   points(c(x0, x1, x2), c(g(x0), g(x1), g(x2)), pch = 19, col = "red", cex = 1.2)
#   abline(v = c(x0, x1, x2), lty = 3, col = "red")
#   segments(x0, g(x0), x1, g(x1), lty = 3, col = "red")
#   segments(x1, g(x1), x2, g(x2), lty = 3, col = "red")
#   polygon(x = c(x0, x0, x1, x2, x2), y = c(0, g(x0), g(x1), g(x2), 0),
#           col = adjustcolor("orange", alpha.f = 0.4), border = NA)
#   abline(h = 0, col = "black")
#   cat("Result = ", round(result, 4))
# }

# Simpson38th12 <- function() {
#   cat("Choose a function to find area\n")
#   cat("1) f(x) = 1 / (1 + x)\n")
#   cat("2) f(x) = | cos(x) |\n")
#   cat("3) f(x) = x^2 + (x^4)/10\n")
#   cat("Enter 1, 2, or 3")
# 
#   choice <- scan(what = integer(), n = 1, quiet = TRUE)
# 
#   if (length(choice) == 0 || !(choice %in% c(1, 2, 3))) {
#     cat("Invalid choice. Exiting.\n")
#     return( )
#   }
# 
#   if (choice == 1) {
#     g <- function(x) 1 / (1 + x)
#   } else if (choice == 2) {
#     g <- function(x) abs(cos(x))
#   } else if (choice == 3) {
#     g <- function(x) x^2 + (x^4)/10
#   }
# 
#   #g <- function(x) 2*x^3-17*x^2+23
#   cat("Enter the lower limit and upper limit\n")
#   ab <- scan(what = numeric(), nmax = 2, quiet = TRUE)
#   if (length(ab) < 2 || any(is.na(ab))) {
#     cat("Invalid input for a and b\n")
#     return()
#   }
#   a <- ab[1]
#   b <- ab[2]
# 
#   h <- (b - a) / 3
#   x0 <- a
#   x1 <- a + h
#   x2 <- a + 2 * h
#   x3 <- b
# 
#   result <- (3 * h / 8) * (g(x0) + 3 * g(x1) + 3 * g(x2) + g(x3))
#   curve(g, from = min(c(a,b)) - 1, to = max(c(a,b)) +1, n = 1000, col = "blue", lwd = 2)
#   points(c(x0, x1, x2, x3), c(g(x0), g(x1), g(x2), g(x3)), pch = 19, col = "red", cex = 1.2)
#   abline(v = c(x0, x1, x2, x3), lty = 3, col = "red")
#   segments(x0, g(x0), x1, g(x1), lty = 3, col = "red")
#   segments(x1, g(x1), x2, g(x2), lty = 3, col = "red")
#   segments(x2, g(x2), x3, g(x3), lty = 3, col = "red")
#   polygon(x = c(x0, x0, x1, x2, x3, x3), y = c(0, g(x0), g(x1), g(x2), g(x3), 0),
#           col = adjustcolor("orange", alpha.f = 0.4), border = NA)
#   abline(h = 0, col = "black")
#   cat("Result = ", round(result, 4))
# }

Trapezoidal10 <- function() {
  cat("Choose a function to find area\n")
  cat("1) f(x) = 1 / (1 + x)\n")
  cat("2) f(x) = | cos(x) |\n")
  cat("3) f(x) = x^2 + (x^4)/10\n")
  cat("Enter 1, 2, or 3")
  
  choice <- scan(what = integer(), n = 1, quiet = TRUE)
  
  if (length(choice) == 0 || !(choice %in% c(1, 2, 3))) {
    cat("Invalid choice for function\n")
    return( )
  }
  
  if (choice == 1) {
    g <- function(x) 1 / (1 + x)
  } else if (choice == 2) {
    g <- function(x) abs(cos(x))
  } else if (choice == 3) {
    g <- function(x) x^2 + (x^4)/10
  }
  
  #g <- function(x) 3*x - 4
  cat("Enter the lower limit and upper limit \n")
  ab <- scan(what = numeric(), nmax = 2, quiet = TRUE)
  if (length(ab) < 2 || any(is.na(ab))) {
    cat("Invalid input for a and b\n")
    return()
  }
  a <- ab[1]
  b <- ab[2]
  
  cat("Enter number of intervals\n")
  n <- scan(what = character(), nmax = 1, quiet = TRUE)
  if (length(n) <= 0 || any(is.na(n))) {
    cat("Invalid input for n detected\n")
    return()
  }
  
  n <- suppressWarnings(as.integer(n))
  if (any(is.na(n)) || n < 1) {
    cat("Invalid input for n detected\n")
    return()
  }
  
  h <- (b-a) / n
  xi <- seq(from = a, to = b, by = h)
  yi <- g(xi)
  
  result <- (h/2) * (yi[1] + yi[n+1] + 2*sum(yi[2:n]))
  curve(g, from = min(c(a,b))-1, to = max(c(a,b))+1, n = 1000, col = "blue", lwd = 2)
  abline(h = 0, col = "black")
  for (i in 1:n) {
    polygon(x = c(xi[i], xi[i], xi[i+1], xi[i+1]), y = c(0, yi[i], yi[i+1], 0),
      col = adjustcolor("orange", alpha.f = 0.35), border = NA)
    segments(xi[i], yi[i], xi[i+1], yi[i+1], col = "red", lty = 3)
    abline(v = xi[i], col = adjustcolor("red", 0.4), lty = 3)
  }
  abline(v = xi[n+1], col = adjustcolor("red", 0.4), lty = 3)
  points(xi, yi, pch = 19, col = "red", cex = 0.8)
  title(main = round(result, 6), col.main = "red")
  cat("Result = ", round(result, 6))
}

Simpson13rd11 <- function() {
  cat("Choose a function to find area\n")
  cat("1) f(x) = 1 / (1 + x)\n")
  cat("2) f(x) = | cos(x) |\n")
  cat("3) f(x) = x^2 + (x^4)/10\n")
  cat("Enter 1, 2, or 3")
  
  choice <- scan(what = integer(), n = 1, quiet = TRUE)
  
  if (length(choice) == 0 || !(choice %in% c(1, 2, 3))) {
    cat("Invalid choice. Exiting.\n")
    return( )
  }
  
  if (choice == 1) {
    g <- function(x) 1 / (1 + x)
  } else if (choice == 2) {
    g <- function(x) abs(cos(x))
  } else if (choice == 3) {
    g <- function(x) x^2 + (x^4)/10
  }
  
  #g <- function(x) 3*x^2 + 11
  #g <- function(x) x^3-2*x^2+2
  cat("Enter the lower limit and upper limit\n")
  ab <- scan(what = numeric(), nmax = 2, quiet = TRUE)
  if (length(ab) < 2 || any(is.na(ab))) {
    cat("Invalid input for a and b\n")
    return()
  }
  a <- ab[1]
  b <- ab[2]
  
  cat("Enter number of intervals\n")
  n <- scan(what = character(), nmax = 1, quiet = TRUE)
  if (length(n) <= 0 || any(is.na(n))) {
    cat("Invalid input for n detected\n")
    return()
  }
  
  n <- suppressWarnings(as.integer(n))
  if (any(is.na(n)) || n < 2 || (n %% 2 != 0)) {
    cat("Invalid input for n detected")
    return()
  }
  
  h <- (b-a) / n
  xi <- seq(from = a, to = b, by = h)
  yi <- g(xi)
  
  oddsum <- sum(yi[seq(from = 2, to = n, by = 2)]) 
  evensum <- sum(yi) - oddsum - yi[1] - yi[n+1]
  result <- (h/3)*(yi[1] + yi[n+1] + 4*oddsum + 2*evensum)
  
  curve(g, from = min(c(a,b))-1, to = max(c(a,b))+1, n = 1000, col = "blue", lwd = 2)
  abline(h = 0, col = "black")
  for (i in 1:n) {
    polygon(x = c(xi[i], xi[i], xi[i+1], xi[i+1]), y = c(0, yi[i], yi[i+1], 0),
            col = adjustcolor("orange", alpha.f = 0.35), border = NA)
    segments(xi[i], yi[i], xi[i+1], yi[i+1], col = "red", lty = 3)
    abline(v = xi[i], col = adjustcolor("red", 0.4), lty = 3)
  }
  abline(v = xi[n+1], col = adjustcolor("red", 0.4), lty = 3)
  points(xi, yi, pch = 19, col = "red", cex = 0.8)
  title(main = round(result, 6), col.main = "red")
  cat("Result = ", round(result, 6))
  
}

Simpson38th12 <- function() {
  cat("Choose a function to find area\n")
  cat("1) f(x) = 1 / (1 + x)\n")
  cat("2) f(x) = | cos(x) |\n")
  cat("3) f(x) = x^2 + (x^4)/10\n")
  cat("Enter 1, 2, or 3")
  
  choice <- scan(what = integer(), n = 1, quiet = TRUE)
  
  if (length(choice) == 0 || !(choice %in% c(1, 2, 3))) {
    cat("Invalid choice. Exiting.\n")
    return( )
  }
  
  if (choice == 1) {
    g <- function(x) 1 / (1 + x)
  } else if (choice == 2) {
    g <- function(x) abs(cos(x))
  } else if (choice == 3) {
    g <- function(x) x^2 + (x^4)/10
  }
  
  #g <- function(x) x^3-2*x^2+2
  cat("Enter the lower limit and upper limit\n")
  ab <- scan(what = numeric(), nmax = 2, quiet = TRUE)
  if (length(ab) < 2 || any(is.na(ab))) {
    cat("Invalid input for a and b\n")
    return()
  }
  a <- ab[1]
  b <- ab[2]
  
  cat("Enter number of intervals\n")
  n <- scan(what = character(), nmax = 1, quiet = TRUE)
  if (length(n) <= 0 || any(is.na(n))) {
    cat("Invalid input for n detected\n")
    return()
  }
  
  n <- suppressWarnings(as.integer(n))
  if (any(is.na(n)) || n < 3 || (n %% 3 != 0)) {
    cat("Invalid input for n detected")
    return()
  }
  
  h <- (b-a) / n
  xi <- seq(from = a, to = b, by = h)
  yi <- g(xi)
  
  multi3 <- sum(yi[seq(from = 3, to = n, by = 3)]) 
  others <- sum(yi) - multi3 - yi[1] - yi[n+1]
  result <- (3*h/8)*(yi[1] + yi[n+1] + 3*others + 2*multi3)
  
  curve(g, from = min(c(a,b))-1, to = max(c(a,b))+1, n = 1000, col = "blue", lwd = 2)
  abline(h = 0, col = "black")
  for (i in 1:n) {
    polygon(x = c(xi[i], xi[i], xi[i+1], xi[i+1]), y = c(0, yi[i], yi[i+1], 0),
            col = adjustcolor("orange", alpha.f = 0.35), border = NA)
    segments(xi[i], yi[i], xi[i+1], yi[i+1], col = "red", lty = 3)
    abline(v = xi[i], col = adjustcolor("red", 0.4), lty = 3)
  }
  abline(v = xi[n+1], col = adjustcolor("red", 0.4), lty = 3)
  points(xi, yi, pch = 19, col = "red", cex = 0.8)
  title(main = round(result, 6), col.main = "red")
  cat("Result = ", round(result, 6))
  
}

Numerical <- function() {
  repeat {
    cat("\n-----Numerical menu-----\n")
    cat("\n1} Bisection method")
    cat("\n2} Secant method")
    cat("\n3} Regula Falsi method")
    cat("\n4} Newton Raphson method")
    cat("\n5} Forward Subs. method")
    cat("\n6} Backwad Subs. method")
    cat("\n7} Lagrange Interp. method")
    cat("\n8} Newton Forward Interp. method")
    cat("\n9} Newton Backward Interp. method")
    cat("\n10} Trapezoidal int. method")
    cat("\n11} Simpson's 1/3rd int. method")
    cat("\n12} Simpson's 3/8th int. method")
    cat("\n13} Exit\n")
    input <- readline("Enter your choice: ")
    choice <- suppressWarnings(as.integer(input))
    
    if (is.na(choice)) {
      cat("\nInvalid choice, please choose a valid option!!\n")
      next
    }
    
    switch(as.character(choice),
           "1"  = BisectionMethod1(),
           "2"  = SecantMethod2(),
           "3"  = RegulaFalsi3(),
           "4"  = NewtonRaph4(),
           "5"  = ForwardSub5(),
           "6"  = BackwardSub6(),
           "7"  = Lagrange7(),
           "8"  = NewtonForward8(),
           "9"  = NewtonBackward9(),
           "10" = Trapezoidal10(),
           "11" = Simpson13rd11(),
           "12" = Simpson38th12(),
           "13" = { cat("\nExited Numerical section\n"); break },
           { cat("\nInvalid choice, please choose a valid option!!\n") }
    )
  }
   
}
############################################################
BuffonNeedle1 <- function() {
  cat("Enter number of iterations\n")
  iterate <- scan(what = character(), nmax = 1, quiet = TRUE)
  if (length(iterate) != 1 || any(is.na(iterate))) {
    cat("Invalid input for n detected\n")
    return()
  }
  
  iterate <- suppressWarnings(as.integer(iterate))
  if (any(is.na(iterate)) || iterate <= 0) {
    cat("Invalid input for n detected")
    return()
  }
  
  count <- 0
  
  cat("Enter 'l'\n")
  l <- scan(what = character(), nmax = 1, quiet = TRUE)
  if (length(l) != 1 || any(is.na(l))) {
    cat("Invalid input for l\n")
    return()
  }
  
  l <- suppressWarnings(as.integer(l))
  if (any(is.na(l)) || l <= 1) {
    cat("Invalid input for l\n")
    return()
  }
  
  cat("Enter 'd'\n")
  d <- scan(what = character(), nmax = 1, quiet = TRUE)
  if (length(d) != 1 || any(is.na(d))) {
    cat("Invalid input for d\n")
    return()
  }
  
  d <- suppressWarnings(as.integer(d))
  if (any(is.na(d)) || d <= 1) {
    cat("Invalid input for d\n")
    return()
  }
  
  y <- runif(n = iterate, min = 0, max = d/2)
  t <- runif(n = iterate, min = 0, max = pi)
  
  for (i in 1:iterate) {
    if (y[i] <= (l/2)*sin(t[i])) {
      count <- count + 1
    }
  }
  
  p = count/iterate
  pie <- (2*l)/(p*d)
  cat("Approximate value of pi =", round(pie, 5))
}

MidSqRNG2 <- function() {
  cat("Enter number of RNG\n")
  RNG <- scan(what = character(), nmax = 1, quiet = TRUE)
  if (length(RNG) != 1 || any(is.na(RNG))) {
    cat("Invalid input for n detected\n")
    return()
  }
  
  RNG <- suppressWarnings(as.integer(RNG))
  if (any(is.na(RNG)) || RNG <= 0) {
    cat("Invalid input for n detected\n")
    return()
  }
  
  cat("Enter a 4 digit seed\n")
  SEED <- scan(what = character(), nmax = 1, quiet = TRUE)
  if (length(SEED) != 1 || any(is.na(SEED))) {
    cat("Invalid input for SEED\n")
    return()
  }
  
  SEED <- suppressWarnings(as.integer(SEED))
  if (any(is.na(SEED)) || nchar(SEED) != 4) {
    cat("Invalid input for 4 digit SEED detected\n")
    return()
  }
  
  result <- numeric(RNG)
  current <- SEED
  
  for (i in 1:RNG) {
    a <- as.integer(current ^ 2)
    #cat("\n",a)
    a <- sprintf("%08d", a)
    current <- as.integer(substr(a, 3, 6))
    result[i] <- current
  }
  
  cat(result)
}

CongRNG3 <- function() {
  cat("Enter number of RNG\n")
  RNG <- scan(what = character(), nmax = 1, quiet = TRUE)
  if (length(RNG) != 1 || any(is.na(RNG))) {
    cat("Invalid input for n detected\n")
    return()
  }
  
  RNG <- suppressWarnings(as.integer(RNG))
  if (any(is.na(RNG)) || RNG <= 0) {
    cat("Invalid input for n detected\n")
    return()
  }
  
  result <- numeric(RNG)
  
  cat("Enter positive 'a' the multiplier\n")
  a <- scan(what = character(), nmax = 1, quiet = TRUE)
  if (length(a) != 1 || any(is.na(a))) {
    cat("Invalid input for multiplier\n")
    return()
  }
  
  a <- suppressWarnings(as.integer(a))
  if (any(is.na(a)) || a < 0) {
    cat("Invalid input for multiplier\n")
    return()
  }
  
  cat("Enter positive 'c' the increment\n")
  c <- scan(what = character(), nmax = 1, quiet = TRUE)
  if (length(c) != 1 || any(is.na(c))) {
    cat("Invalid input for increment\n")
    return()
  }
  
  c <- suppressWarnings(as.integer(c))
  if (any(is.na(c)) || c < 0) {
    cat("Invalid input for increment\n")
    return()
  }
  
  cat("Enter positive 'm' the modulo divisor\n")
  m <- scan(what = character(), nmax = 1, quiet = TRUE)
  if (length(m) != 1 || any(is.na(m))) {
    cat("Invalid input for divisor\n")
    return()
  }
  
  m <- suppressWarnings(as.integer(m))
  if (any(is.na(m)) || m < 0) {
    cat("Invalid input for divisor\n")
    return()
  }
  
  cat("Enter a seed (x)\n")
  SEED <- scan(what = character(), nmax = 1, quiet = TRUE)
  if (length(SEED) != 1 || any(is.na(SEED))) {
    cat("Invalid input for SEED\n")
    return()
  }
  
  SEED <- suppressWarnings(as.integer(SEED))
  if (any(is.na(SEED))) {
    cat("Invalid input for SEED detected\n")
    return()
  }
  
  prime_vect <- Primes(m)
  prime_vect_m <- c()
  for (i in prime_vect) {
    if ((m %% i) == 0) {
      prime_vect_m <- c(prime_vect_m, i)
    }
  }
  
  cond1 <- (GCD(c, m) == 1)
  cond2 <- all(((a-1) %% prime_vect_m) == 0)
  cond3 <- (m %% 4 != 0) || ((a - 1) %% 4 == 0)
  
  if (cond1 && cond2 && cond3) {
    for (i in 1:RNG) {
      SEED <- (a * SEED + c) %% m
      result[i] <- SEED
    }
    cat("Knuth conditions met\n")
    print(result)
  } else {
    for (i in 1:RNG) {
      SEED <- (a * SEED + c) %% m
      result[i] <- SEED
    }
    cat("Knuth conditions not met\n")
    print(result)
  }
}

MultiRNG4 <- function() {
  cat("Enter number of RNG\n")
  RNG <- scan(what = character(), nmax = 1, quiet = TRUE)
  if (length(RNG) != 1 || any(is.na(RNG))) {
    cat("Invalid input for n detected\n")
    return()
  }
  
  RNG <- suppressWarnings(as.integer(RNG))
  if (any(is.na(RNG)) || RNG <= 0) {
    cat("Invalid input for n detected")
    return()
  }
  
  result <- numeric(RNG)
  
  cat("Enter positive 'a' the multiplier\n")
  a <- scan(what = character(), nmax = 1, quiet = TRUE)
  if (length(a) != 1 || any(is.na(a))) {
    cat("Invalid input for multiplier\n")
    return()
  }
  
  a <- suppressWarnings(as.integer(a))
  if (any(is.na(a)) || a < 0) {
    cat("Invalid input for multiplier\n")
    return()
  }
  
  
  cat("Enter positive 'm' the modulo divisor\n")
  m <- scan(what = character(), nmax = 1, quiet = TRUE)
  if (length(m) != 1 || any(is.na(m))) {
    cat("Invalid input for divisor\n")
    return()
  }
  
  m <- suppressWarnings(as.integer(m))
  if (any(is.na(m)) || m < 0) {
    cat("Invalid input for divisor\n")
    return()
  }
  
  cat("Enter a seed (x)\n")
  SEED <- scan(what = character(), nmax = 1, quiet = TRUE)
  if (length(SEED) != 1 || any(is.na(SEED))) {
    cat("Invalid input for SEED\n")
    return()
  }
  
  SEED <- suppressWarnings(as.integer(SEED))
  if (any(is.na(SEED))) {
    cat("Invalid input for SEED detected\n")
    return()
  }
  
  for (i in 1:RNG) {
    SEED <- (a * SEED) %% m
    result[i] <- SEED
  }
  print(result)
  
}

Simulation <- function() {
  repeat {
    cat("\n-----Simulation menu-----\n")
    cat("\n1] Buffon Needle Pi approximation")
    cat("\n2] Mid-Square RNG")
    cat("\n3] Congruential RNG")
    cat("\n4] Multiplicative RNG")
    cat("\n5] Exit\n")
    input <- readline("Enter your choice: ")
    choice <- suppressWarnings(as.integer(input))
    
    if (is.na(choice)) {
      cat("\nInvalid choice, please choose a valid option!!\n")
      next
    }
    
    switch (as.character(choice),
      "1" = BuffonNeedle1(),
      "2" = MidSqRNG2(),
      "3" = CongRNG3(),
      "4" = MultiRNG4(),
      "5" = { cat("\nExited Calculator\n"); break },
      { cat("\nInvalid choice, please enter 1-5.\n") }
    )
  }
}
############################################################
MAIN <- function() {
  repeat {
    cat("\n-----Calculator menu-----\n")
    cat("1) Mathematics\n")
    cat("2) Descriptive Stats.\n")
    cat("3) Numerical\n")
    cat("4) Simulation\n")
    cat("5) Exit\n")
    input <- readline("Enter an integer: ")
    choice <- suppressWarnings(as.integer(input))
    
    if (is.na(choice)) {
      cat("\nInvalid choice, please enter 1-5.\n")
      next
    }
    
    switch (as.character(choice),
            "1" = Mathematics(),
            "2" = DescStats(),
            "3" = Numerical(),
            "4" = Simulation(),
            "5" = { cat("\nExited Calculator\n"); break },
            { cat("\nInvalid choice, please enter 1-5.\n") }
    )
  }
   
}

MAIN()