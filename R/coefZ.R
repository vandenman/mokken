# Aangepast op 6 maart 2020
#"coefZ.old" <-
# Functie voor toevoegen van lowerbound
#function(X){
#   X <- check.data(X)
#    N <- nrow(X)
#    S <- var(X)
#    Sij <- outer(apply(X,2,var),apply(X,2,var),"*")
#    Zij <- (S * sqrt(N-1))/sqrt(Sij)
#    diag(S) <- diag(Sij) <- diag(Zij) <- 0
#    Zi <- (apply(S,1,sum) * sqrt(N-1))/ sqrt(apply(Sij,1,sum))       
#    Z  <- (sum(S)/2 * sqrt(N-1))/ sqrt(sum(Sij)/2)
#    return(list(Zij=Zij,Zi=Zi,Z=Z))
#}


"coefZ" <-
function(X, lowerbound = 0){
     #X <- check.data(X)
     N <- nrow(X)
     S <- var(X)
     Smax <- var(apply(X,2,sort))
     Sij <- outer(apply(X,2,var),apply(X,2,var),"*")
     Zij <- sqrt(N-1) * (S/Smax - lowerbound)/(sqrt(Sij)/Smax)
     diag(S) <- diag(Sij) <- diag(Zij) <- 0
     Zi <- sqrt(N-1) * (apply(S,1,sum) / apply(Smax,1,sum) - lowerbound)/ (sqrt(apply(Sij,1,sum))  /  (apply(Smax,1,sum) )    )
     Z  <- sqrt(N-1) * (sum(S) / (sum(Smax)) - lowerbound)/ (sqrt(sum(Sij)/2) / (sum(Smax) / 2) )
     return(list(Zij=Zij,Zi=Zi,Z=Z))
}

