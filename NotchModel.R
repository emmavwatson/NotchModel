
#requires dplyr
#requires reshape2
#require leftout_template40x40.txt:
leftout_template40x40 <- read.delim("leftout_template40x40.txt")
#define size of field: 
cellnum <- 40
#adjust poising factor of each group with the code call at the bottom:
#fwa_multi_wt_vs_wt <- iteration_notch_sim(2.2, 1)

notch_model1_noise32_FINALnoplots <- function(iter, Kd, Kh, Kn, uDm, uHm, uFm, uDp, uHp, uFp, uN, v, mul, mul2, mul3, mul4, leftoutY) {
    
    Heatmap_Matrix_N <- data.frame(matrix(runif(cellnum*cellnum, min=0, max=1), ncol=cellnum))
    Heatmap_Matrix_Dm <- data.frame(matrix(runif(cellnum*cellnum, min=0, max=1), ncol=cellnum))
    Heatmap_Matrix_Hm <- data.frame(matrix(runif(cellnum*cellnum, min=0, max=0.1), ncol=cellnum))
    Heatmap_Matrix_Fm <- data.frame(matrix(runif(cellnum*cellnum, min=0, max=1), ncol=cellnum))
    Heatmap_Matrix_Dp <- data.frame(matrix(runif(cellnum*cellnum, min=0, max=1), ncol=cellnum))
    Heatmap_Matrix_Hp <- data.frame(matrix(runif(cellnum*cellnum, min=0, max=0.1), ncol=cellnum))
    Heatmap_Matrix_Fp <- data.frame(matrix(runif(cellnum*cellnum, min=0, max=1), ncol=cellnum))
    
    Heatmap_Matrix_N_new <- data.frame(matrix(nrow=cellnum, ncol=cellnum))
    Heatmap_Matrix_Dm_new <- data.frame(matrix(nrow=cellnum, ncol=cellnum))
    Heatmap_Matrix_Hm_new <- data.frame(matrix(nrow=cellnum, ncol=cellnum))
    Heatmap_Matrix_Fm_new <- data.frame(matrix(nrow=cellnum, ncol=cellnum))
    Heatmap_Matrix_Dp_new <- data.frame(matrix(nrow=cellnum, ncol=cellnum))
    Heatmap_Matrix_Hp_new <- data.frame(matrix(nrow=cellnum, ncol=cellnum))
    Heatmap_Matrix_Fp_new <- data.frame(matrix(nrow=cellnum, ncol=cellnum))
    
    #iter <- 1000
    #Kd <- 0.1
    #Kh <- 0.01
    #Kn <- 1
    #uDm <- 0.01
    #uHm <- 0.1
    #uFm <- 0.1
    #uDp <- 0.1
    #uHp <- 0.1
    #uFp <- 0.1
    #uN <- 0.01
    #v <- 30
    #mul <- 3
    #mul2 <- 0.1
    #mul3 <- 0.3
    #mul4 <- 1
    Heatmap_Matrix_N1 <- as.matrix(Heatmap_Matrix_N)
    
    for (j in seq(from=2, to=(cellnum-2), by=2)) {
        for (i in 2:(cellnum-1)) {
            N0 <- Heatmap_Matrix_N[j,i]
            Dm <- Heatmap_Matrix_Dm[j,i]
            Dp <- Heatmap_Matrix_Dp[j,i]
            Hm <- Heatmap_Matrix_Hm[j,i]
            Hp <- Heatmap_Matrix_Hp[j,i]
            Fm <- Heatmap_Matrix_Fm[j,i]
            Fp <- Heatmap_Matrix_Fp[j,i]
            D1 <- Heatmap_Matrix_Dp[j-1,i] 
            D2 <- Heatmap_Matrix_Dp[j-1,i-1] 
            D3 <- Heatmap_Matrix_Dp[j,i-1] 
            D4 <- Heatmap_Matrix_Dp[j,i+1] 
            D5 <- Heatmap_Matrix_Dp[j+1,i] 
            D6 <- Heatmap_Matrix_Dp[j+1,i-1]
            N1 <- Heatmap_Matrix_N[j-1,i] 
            N2 <- Heatmap_Matrix_N[j-1,i-1] 
            N3 <- Heatmap_Matrix_N[j,i-1] 
            N4 <- Heatmap_Matrix_N[j,i+1] 
            N5 <- Heatmap_Matrix_N[j+1,i] 
            N6 <- Heatmap_Matrix_N[j+1,i-1] 
            
            Navg <- (N1+N2+N3+N4+N5+N6)/Kd
            Davg <- (D1+D2+D3+D4+D5+D6)/Kd
            Dm_new <- Dm + uDm*(-1*Dm + 1/(1+((Hp/Kh)*(Hp/Kh))))
            Dp_new <- Dp + uDp*((-1*Dp + Dm) )
            Hm_new <- Hm + uHm*(-1*Hm + ((N0/Kn)*(N0/Kn))/(1+((N0/Kn))*(N0/Kn)))
            Hp_new <- Hp + uHp*(-1*Hp + Hm)
            Fm_new <- Fm + uFm*(-1*Fm + 1/(1 + ((Hp/Kh)*(Hp/Kh)) + ((Hp/Kh)*(Hp/Kh))))
            Fp_new <- Fp + uFp*(-1*Fp + Fm)
            N0_new <- N0 + uN*(-1*(1+ v*Fp)*N0 + mul3*((N0/Kn)*(N0/Kn))/(1+((N0/Kn))*(N0/Kn)) + mul4*(Davg*Davg)/(1+(Davg*Davg)) - mul2*((Navg*Navg)/(1+(Navg*Navg))))
            
            N0_new <- N0_new + (runif(1, min=-0.05, max=0.05))*N0_new
            N0_new <- ifelse(N0_new > 1, 1, N0_new)
            
            Heatmap_Matrix_N_new[j,i] <- N0_new
            Heatmap_Matrix_Dm_new[j,i] <- Dm_new
            Heatmap_Matrix_Hm_new[j,i] <- Hm_new
            Heatmap_Matrix_Fm_new[j,i] <- Fm_new
            Heatmap_Matrix_Dp_new[j,i] <- Dp_new
            Heatmap_Matrix_Hp_new[j,i] <- Hp_new
            Heatmap_Matrix_Fp_new[j,i] <- Fp_new
        }
    }
    for (j in seq(from=3, to=(cellnum-1), by=2)) {
        for (i in 2:(cellnum-1)) {
            N0 <- Heatmap_Matrix_N[j,i]
            Dm <- Heatmap_Matrix_Dm[j,i]
            Dp <- Heatmap_Matrix_Dp[j,i]
            Hm <- Heatmap_Matrix_Hm[j,i]
            Hp <- Heatmap_Matrix_Hp[j,i]
            Fm <- Heatmap_Matrix_Fm[j,i]
            Fp <- Heatmap_Matrix_Fp[j,i]
            D1 <- Heatmap_Matrix_Dp[j-1,i] 
            D2 <- Heatmap_Matrix_Dp[j-1,i+1] 
            D3 <- Heatmap_Matrix_Dp[j,i-1] 
            D4 <- Heatmap_Matrix_Dp[j,i+1] 
            D5 <- Heatmap_Matrix_Dp[j+1,i] 
            D6 <- Heatmap_Matrix_Dp[j+1,i+1] 
            N1 <- Heatmap_Matrix_N[j-1,i] 
            N2 <- Heatmap_Matrix_N[j-1,i+1] 
            N3 <- Heatmap_Matrix_N[j,i-1] 
            N4 <- Heatmap_Matrix_N[j,i+1] 
            N5 <- Heatmap_Matrix_N[j+1,i] 
            N6 <- Heatmap_Matrix_N[j+1,i+1] 
            
            Navg <- (N1+N2+N3+N4+N5+N6)/Kd
            Davg <- (D1+D2+D3+D4+D5+D6)/Kd
            Dm_new <- Dm + uDm*(-1*Dm + 1/(1+((Hp/Kh)*(Hp/Kh))))
            Dp_new <- Dp + uDp*((-1*Dp + Dm) )
            Hm_new <- Hm + uHm*(-1*Hm + ((N0/Kn)*(N0/Kn))/(1+((N0/Kn))*(N0/Kn)))
            Hp_new <- Hp + uHp*(-1*Hp + Hm)
            Fm_new <- Fm + uFm*(-1*Fm + 1/(1 + ((Hp/Kh)*(Hp/Kh)) + ((Hp/Kh)*(Hp/Kh))))
            Fp_new <- Fp + uFp*(-1*Fp + Fm)
            N0_new <- N0 + uN*(-1*(1+ v*Fp)*N0 + mul3*((N0/Kn)*(N0/Kn))/(1+((N0/Kn))*(N0/Kn)) + mul4*(Davg*Davg)/(1+(Davg*Davg)) - mul2*((Navg*Navg)/(1+(Navg*Navg))))
           
            N0_new <- N0_new + (runif(1, min=-0.05, max=0.05))*N0_new
            N0_new <- ifelse(N0_new > 1, 1, N0_new)
            
            Heatmap_Matrix_N_new[j,i] <- N0_new
            Heatmap_Matrix_Dm_new[j,i] <- Dm_new
            Heatmap_Matrix_Hm_new[j,i] <- Hm_new
            Heatmap_Matrix_Fm_new[j,i] <- Fm_new
            Heatmap_Matrix_Dp_new[j,i] <- Dp_new
            Heatmap_Matrix_Hp_new[j,i] <- Hp_new
            Heatmap_Matrix_Fp_new[j,i] <- Fp_new
        }
    }
    Heatmap_Matrix_N_new[1,] <- Heatmap_Matrix_N[1,]
    Heatmap_Matrix_Dp_new[1,] <- Heatmap_Matrix_Dp[1,]
    Heatmap_Matrix_N_new[,1] <- Heatmap_Matrix_N[,1]
    Heatmap_Matrix_Dp_new[,1] <- Heatmap_Matrix_Dp[,1]
    Heatmap_Matrix_N_new[cellnum,] <- Heatmap_Matrix_N[cellnum,]
    Heatmap_Matrix_Dp_new[cellnum,] <- Heatmap_Matrix_Dp[cellnum,]
    Heatmap_Matrix_N_new[,cellnum] <- Heatmap_Matrix_N[,cellnum]
    Heatmap_Matrix_Dp_new[,cellnum] <- Heatmap_Matrix_Dp[,cellnum]
    Heatmap_Matrix_Dm_new[1,] <- Heatmap_Matrix_Dm[1,]
    Heatmap_Matrix_Dm_new[,1] <- Heatmap_Matrix_Dm[,1]
    Heatmap_Matrix_Dm_new[cellnum,] <- Heatmap_Matrix_Dm[cellnum,]
    Heatmap_Matrix_Dm_new[,cellnum] <- Heatmap_Matrix_Dm[,cellnum]
    Heatmap_Matrix_Hm_new[1,] <- Heatmap_Matrix_Hm[1,]
    Heatmap_Matrix_Hm_new[,1] <- Heatmap_Matrix_Hm[,1]
    Heatmap_Matrix_Hm_new[cellnum,] <- Heatmap_Matrix_Hm[cellnum,]
    Heatmap_Matrix_Hm_new[,cellnum] <- Heatmap_Matrix_Hm[,cellnum]
    Heatmap_Matrix_Hp_new[1,] <- Heatmap_Matrix_Hp[1,]
    Heatmap_Matrix_Hp_new[,1] <- Heatmap_Matrix_Hp[,1]
    Heatmap_Matrix_Hp_new[cellnum,] <- Heatmap_Matrix_Hp[cellnum,]
    Heatmap_Matrix_Hp_new[,cellnum] <- Heatmap_Matrix_Hp[,cellnum]
    Heatmap_Matrix_Fm_new[1,] <- Heatmap_Matrix_Fm[1,]
    Heatmap_Matrix_Fm_new[,1] <- Heatmap_Matrix_Fm[,1]
    Heatmap_Matrix_Fm_new[cellnum,] <- Heatmap_Matrix_Fm[cellnum,]
    Heatmap_Matrix_Fm_new[,cellnum] <- Heatmap_Matrix_Fm[,cellnum]
    Heatmap_Matrix_Fp_new[1,] <- Heatmap_Matrix_Fp[1,]
    Heatmap_Matrix_Fp_new[,1] <- Heatmap_Matrix_Fp[,1]
    Heatmap_Matrix_Fp_new[cellnum,] <- Heatmap_Matrix_Fp[cellnum,]
    Heatmap_Matrix_Fp_new[,cellnum] <- Heatmap_Matrix_Fp[,cellnum]
    Heatmap_Matrix_N_new <- as.matrix(Heatmap_Matrix_N_new)
    Heatmap_Matrix_Dp_new <- as.matrix(Heatmap_Matrix_Dp_new)
    Heatmap_Matrix_Dm_new <- as.matrix(Heatmap_Matrix_Dm_new)
    Heatmap_Matrix_Hm_new <- as.matrix(Heatmap_Matrix_Hm_new)
    Heatmap_Matrix_Hp_new <- as.matrix(Heatmap_Matrix_Hp_new)
    Heatmap_Matrix_Fm_new <- as.matrix(Heatmap_Matrix_Fm_new)
    Heatmap_Matrix_Fp_new <- as.matrix(Heatmap_Matrix_Fp_new)
    
    
    for(m in 1:iter) {
        Heatmap_Matrix_N <- Heatmap_Matrix_N_new
        Heatmap_Matrix_Dm <- Heatmap_Matrix_Dm_new
        Heatmap_Matrix_Hm <- Heatmap_Matrix_Hm_new
        Heatmap_Matrix_Fm <- Heatmap_Matrix_Fm_new
        Heatmap_Matrix_Dp <- Heatmap_Matrix_Dp_new
        Heatmap_Matrix_Hp <- Heatmap_Matrix_Hp_new
        Heatmap_Matrix_Fp <- Heatmap_Matrix_Fp_new
        
        Heatmap_Matrix_N_new <- data.frame(matrix(nrow=cellnum, ncol=cellnum))
        Heatmap_Matrix_Dm_new <- data.frame(matrix(nrow=cellnum, ncol=cellnum))
        Heatmap_Matrix_Hm_new <- data.frame(matrix(nrow=cellnum, ncol=cellnum))
        Heatmap_Matrix_Fm_new <- data.frame(matrix(nrow=cellnum, ncol=cellnum))
        Heatmap_Matrix_Dp_new <- data.frame(matrix(nrow=cellnum, ncol=cellnum))
        Heatmap_Matrix_Hp_new <- data.frame(matrix(nrow=cellnum, ncol=cellnum))
        Heatmap_Matrix_Fp_new <- data.frame(matrix(nrow=cellnum, ncol=cellnum))
        
        
        for (j in seq(from=2, to=(cellnum-2), by=2)) {
            for (i in 2:(cellnum-1)) {
                
                N0 <- Heatmap_Matrix_N[j,i]
                Dm <- Heatmap_Matrix_Dm[j,i]
                Dp <- Heatmap_Matrix_Dp[j,i]
                Hm <- Heatmap_Matrix_Hm[j,i]
                Hp <- Heatmap_Matrix_Hp[j,i]
                Fm <- Heatmap_Matrix_Fm[j,i]
                Fp <- Heatmap_Matrix_Fp[j,i]
                D1 <- Heatmap_Matrix_Dp[j-1,i] 
                D2 <- Heatmap_Matrix_Dp[j-1,i-1] 
                D3 <- Heatmap_Matrix_Dp[j,i-1] 
                D4 <- Heatmap_Matrix_Dp[j,i+1] 
                D5 <- Heatmap_Matrix_Dp[j+1,i] 
                D6 <- Heatmap_Matrix_Dp[j+1,i-1] 
                N1 <- Heatmap_Matrix_N[j-1,i] 
                N2 <- Heatmap_Matrix_N[j-1,i-1] 
                N3 <- Heatmap_Matrix_N[j,i-1] 
                N4 <- Heatmap_Matrix_N[j,i+1] 
                N5 <- Heatmap_Matrix_N[j+1,i] 
                N6 <- Heatmap_Matrix_N[j+1,i-1] 
                
                Navg <- (N1+N2+N3+N4+N5+N6)/Kd
                Davg <- (D1+D2+D3+D4+D5+D6)/Kd
                Dm_new <- Dm + uDm*(-1*Dm + 1/(1+((Hp/Kh)*(Hp/Kh))))
                Dp_new <- Dp + uDp*((-1*Dp + Dm) )
                Hm_new <- Hm + uHm*(-1*Hm + ((N0/Kn)*(N0/Kn))/(1+((N0/Kn))*(N0/Kn)))
                Hp_new <- Hp + uHp*(-1*Hp + Hm)
                Fm_new <- Fm + uFm*(-1*Fm + 1/(1 + ((Hp/Kh)*(Hp/Kh)) + ((Hp/Kh)*(Hp/Kh))))
                Fp_new <- Fp + uFp*(-1*Fp + Fm)
                N0_new <- N0 + uN*(-1*(1+ v*Fp)*N0 + mul3*((N0/Kn)*(N0/Kn))/(1+((N0/Kn))*(N0/Kn)) + mul4*(Davg*Davg)/(1+(Davg*Davg)) - mul2*((Navg*Navg)/(1+(Navg*Navg))))
                
                N0_new <- N0_new + (runif(1, min=-0.05, max=0.05))*N0_new
                N0_new <- ifelse(N0_new > 1, 1, N0_new)
                
                Heatmap_Matrix_N_new[j,i] <- N0_new
                Heatmap_Matrix_Dm_new[j,i] <- Dm_new
                Heatmap_Matrix_Hm_new[j,i] <- Hm_new
                Heatmap_Matrix_Fm_new[j,i] <- Fm_new
                Heatmap_Matrix_Dp_new[j,i] <- Dp_new
                Heatmap_Matrix_Hp_new[j,i] <- Hp_new
                Heatmap_Matrix_Fp_new[j,i] <- Fp_new
            }
        }
        for (j in seq(from=3, to=(cellnum-1), by=2)) {
            for (i in 2:(cellnum-1)) {
                
                N0 <- Heatmap_Matrix_N[j,i]
                Dm <- Heatmap_Matrix_Dm[j,i]
                Dp <- Heatmap_Matrix_Dp[j,i]
                Hm <- Heatmap_Matrix_Hm[j,i]
                Hp <- Heatmap_Matrix_Hp[j,i]
                Fm <- Heatmap_Matrix_Fm[j,i]
                Fp <- Heatmap_Matrix_Fp[j,i]
                D1 <- Heatmap_Matrix_Dp[j-1,i] 
                D2 <- Heatmap_Matrix_Dp[j-1,i+1] 
                D3 <- Heatmap_Matrix_Dp[j,i-1] 
                D4 <- Heatmap_Matrix_Dp[j,i+1] 
                D5 <- Heatmap_Matrix_Dp[j+1,i] 
                D6 <- Heatmap_Matrix_Dp[j+1,i+1] 
                N1 <- Heatmap_Matrix_N[j-1,i] 
                N2 <- Heatmap_Matrix_N[j-1,i+1] 
                N3 <- Heatmap_Matrix_N[j,i-1] 
                N4 <- Heatmap_Matrix_N[j,i+1] 
                N5 <- Heatmap_Matrix_N[j+1,i] 
                N6 <- Heatmap_Matrix_N[j+1,i+1] 
                
                Navg <- (N1+N2+N3+N4+N5+N6)/Kd
                Davg <- (D1+D2+D3+D4+D5+D6)/Kd
                Dm_new <- Dm + uDm*(-1*Dm + 1/(1+((Hp/Kh)*(Hp/Kh))))
                Dp_new <- Dp + uDp*((-1*Dp + Dm) )
                Hm_new <- Hm + uHm*(-1*Hm + ((N0/Kn)*(N0/Kn))/(1+((N0/Kn))*(N0/Kn)))
                Hp_new <- Hp + uHp*(-1*Hp + Hm)
                Fm_new <- Fm + uFm*(-1*Fm + 1/(1 + ((Hp/Kh)*(Hp/Kh)) + ((Hp/Kh)*(Hp/Kh))))
                Fp_new <- Fp + uFp*(-1*Fp + Fm)
                N0_new <- N0 + uN*(-1*(1+ v*Fp)*N0 + mul3*((N0/Kn)*(N0/Kn))/(1+((N0/Kn))*(N0/Kn)) + mul4*(Davg*Davg)/(1+(Davg*Davg)) - mul2*((Navg*Navg)/(1+(Navg*Navg))))
                
                N0_new <- N0_new + (runif(1, min=-0.05, max=0.05))*N0_new
                N0_new <- ifelse(N0_new > 1, 1, N0_new)
                
                Heatmap_Matrix_N_new[j,i] <- N0_new
                Heatmap_Matrix_Dm_new[j,i] <- Dm_new
                Heatmap_Matrix_Hm_new[j,i] <- Hm_new
                Heatmap_Matrix_Fm_new[j,i] <- Fm_new
                Heatmap_Matrix_Dp_new[j,i] <- Dp_new
                Heatmap_Matrix_Hp_new[j,i] <- Hp_new
                Heatmap_Matrix_Fp_new[j,i] <- Fp_new
            }
        }
        #leftout <- data.frame(matrix( 
        #   c(sample(2:(cellnum-1), b, replace=T), sample(2:(cellnum-1), b, replace=T)), 
        #  nrow=b, 
        # ncol=2))
        leftout2 <- leftoutY[leftoutY$row %% 2 == 0, ]
        leftout2 <- subset(leftout2, row != 1 & row != 40 & column != 1 & column != 40)
        for(k in 1:nrow(leftout2)) {
            j <- leftout2[k,1]
            i <- leftout2[k,2]
            N0 <- Heatmap_Matrix_N[j,i]
            Dm <- Heatmap_Matrix_Dm[j,i]
            Dp <- Heatmap_Matrix_Dp[j,i]
            Hm <- Heatmap_Matrix_Hm[j,i]
            Hp <- Heatmap_Matrix_Hp[j,i]
            Fm <- Heatmap_Matrix_Fm[j,i]
            Fp <- Heatmap_Matrix_Fp[j,i]
            D1 <- Heatmap_Matrix_Dp[j-1,i] 
            D2 <- Heatmap_Matrix_Dp[j-1,i-1] 
            D3 <- Heatmap_Matrix_Dp[j,i-1] 
            D4 <- Heatmap_Matrix_Dp[j,i+1] 
            D5 <- Heatmap_Matrix_Dp[j+1,i] 
            D6 <- Heatmap_Matrix_Dp[j+1,i-1] 
            N1 <- Heatmap_Matrix_N[j-1,i] 
            N2 <- Heatmap_Matrix_N[j-1,i-1] 
            N3 <- Heatmap_Matrix_N[j,i-1] 
            N4 <- Heatmap_Matrix_N[j,i+1] 
            N5 <- Heatmap_Matrix_N[j+1,i] 
            N6 <- Heatmap_Matrix_N[j+1,i-1] 
            
            Navg <- (N1+N2+N3+N4+N5+N6)/Kd
            Davg <- (D1+D2+D3+D4+D5+D6)/Kd
            Dm_new <- Dm + uDm*(-1*Dm + 1/(1+((Hp/Kh)*(Hp/Kh))))
            Dp_new <- Dp + uDp*((-1*Dp + Dm) )
            Hm_new <- Hm + uHm*(-1*Hm + ((N0/Kn)*(N0/Kn))/(1+((N0/Kn))*(N0/Kn)))
            Hp_new <- Hp + uHp*(-1*Hp + Hm)
            Fm_new <- Fm + uFm*(-1*Fm + 1/(1 + ((Hp/Kh)*(Hp/Kh)) + ((Hp/Kh)*(Hp/Kh))))
            Fp_new <- Fp + uFp*(-1*Fp + Fm)
            N0_new <- N0 + uN*(-1*(1+ v*Fp)*N0 + mul3*((N0/Kn)*(N0/Kn))/(1+((N0/Kn))*(N0/Kn)) + mul*((Davg*Davg)/(1+(Davg*Davg))) - mul2*((Navg*Navg)/(1+(Navg*Navg))))
         
            N0_new <- N0_new + (runif(1, min=-0.05, max=0.05))*N0_new
            N0_new <- ifelse(N0_new > 1, 1, N0_new)
            
            Heatmap_Matrix_N_new[j,i] <- N0_new
            Heatmap_Matrix_Dm_new[j,i] <- Dm_new
            Heatmap_Matrix_Hm_new[j,i] <- Hm_new
            Heatmap_Matrix_Fm_new[j,i] <- Fm_new
            Heatmap_Matrix_Dp_new[j,i] <- Dp_new
            Heatmap_Matrix_Hp_new[j,i] <- Hp_new
            Heatmap_Matrix_Fp_new[j,i] <- Fp_new
        }
        leftout3 <- leftoutY[leftoutY$row %% 2 == 1, ]
        leftout3 <- subset(leftout3, row != 1 & row != 40 & column != 1 & column != 40)
        for(k in 1:nrow(leftout3)) {
            j <- leftout3[k,1]
            i <- leftout3[k,2]
            N0 <- Heatmap_Matrix_N[j,i]
            Dm <- Heatmap_Matrix_Dm[j,i]
            Dp <- Heatmap_Matrix_Dp[j,i]
            Hm <- Heatmap_Matrix_Hm[j,i]
            Hp <- Heatmap_Matrix_Hp[j,i]
            Fm <- Heatmap_Matrix_Fm[j,i]
            Fp <- Heatmap_Matrix_Fp[j,i]
            D1 <- Heatmap_Matrix_Dp[j-1,i] 
            D2 <- Heatmap_Matrix_Dp[j-1,i+1] 
            D3 <- Heatmap_Matrix_Dp[j,i-1] 
            D4 <- Heatmap_Matrix_Dp[j,i+1] 
            D5 <- Heatmap_Matrix_Dp[j+1,i] 
            D6 <- Heatmap_Matrix_Dp[j+1,i+1] 
            N1 <- Heatmap_Matrix_N[j-1,i] 
            N2 <- Heatmap_Matrix_N[j-1,i+1] 
            N3 <- Heatmap_Matrix_N[j,i-1] 
            N4 <- Heatmap_Matrix_N[j,i+1] 
            N5 <- Heatmap_Matrix_N[j+1,i] 
            N6 <- Heatmap_Matrix_N[j+1,i+1] 
            
            Navg <- (N1+N2+N3+N4+N5+N6)/Kd
            Davg <- (D1+D2+D3+D4+D5+D6)/Kd
            Dm_new <- Dm + uDm*(-1*Dm + 1/(1+((Hp/Kh)*(Hp/Kh))))
            Dp_new <- Dp + uDp*((-1*Dp + Dm) )
            Hm_new <- Hm + uHm*(-1*Hm + ((N0/Kn)*(N0/Kn))/(1+((N0/Kn))*(N0/Kn)))
            Hp_new <- Hp + uHp*(-1*Hp + Hm)
            Fm_new <- Fm + uFm*(-1*Fm + 1/(1 + ((Hp/Kh)*(Hp/Kh)) + ((Hp/Kh)*(Hp/Kh))))
            Fp_new <- Fp + uFp*(-1*Fp + Fm)
            N0_new <- N0 + uN*(-1*(1+ v*Fp)*N0 + mul3*((N0/Kn)*(N0/Kn))/(1+((N0/Kn))*(N0/Kn)) + mul*((Davg*Davg)/(1+(Davg*Davg))) - mul2*((Navg*Navg)/(1+(Navg*Navg))))
            
            N0_new <- N0_new + (runif(1, min=-0.05, max=0.05))*N0_new
            N0_new <- ifelse(N0_new > 1, 1, N0_new)
            
            Heatmap_Matrix_N_new[j,i] <- N0_new
            Heatmap_Matrix_Dm_new[j,i] <- Dm_new
            Heatmap_Matrix_Hm_new[j,i] <- Hm_new
            Heatmap_Matrix_Fm_new[j,i] <- Fm_new
            Heatmap_Matrix_Dp_new[j,i] <- Dp_new
            Heatmap_Matrix_Hp_new[j,i] <- Hp_new
            Heatmap_Matrix_Fp_new[j,i] <- Fp_new
        }
        
        Heatmap_Matrix_N_new[1,] <- Heatmap_Matrix_N[1,]
        Heatmap_Matrix_Dp_new[1,] <- Heatmap_Matrix_Dp[1,]
        Heatmap_Matrix_N_new[,1] <- Heatmap_Matrix_N[,1]
        Heatmap_Matrix_Dp_new[,1] <- Heatmap_Matrix_Dp[,1]
        Heatmap_Matrix_N_new[cellnum,] <- Heatmap_Matrix_N[cellnum,]
        Heatmap_Matrix_Dp_new[cellnum,] <- Heatmap_Matrix_Dp[cellnum,]
        Heatmap_Matrix_N_new[,cellnum] <- Heatmap_Matrix_N[,cellnum]
        Heatmap_Matrix_Dp_new[,cellnum] <- Heatmap_Matrix_Dp[,cellnum]
        Heatmap_Matrix_Dm_new[1,] <- Heatmap_Matrix_Dm[1,]
        Heatmap_Matrix_Dm_new[,1] <- Heatmap_Matrix_Dm[,1]
        Heatmap_Matrix_Dm_new[cellnum,] <- Heatmap_Matrix_Dm[cellnum,]
        Heatmap_Matrix_Dm_new[,cellnum] <- Heatmap_Matrix_Dm[,cellnum]
        Heatmap_Matrix_Hm_new[1,] <- Heatmap_Matrix_Hm[1,]
        Heatmap_Matrix_Hm_new[,1] <- Heatmap_Matrix_Hm[,1]
        Heatmap_Matrix_Hm_new[cellnum,] <- Heatmap_Matrix_Hm[cellnum,]
        Heatmap_Matrix_Hm_new[,cellnum] <- Heatmap_Matrix_Hm[,cellnum]
        Heatmap_Matrix_Hp_new[1,] <- Heatmap_Matrix_Hp[1,]
        Heatmap_Matrix_Hp_new[,1] <- Heatmap_Matrix_Hp[,1]
        Heatmap_Matrix_Hp_new[cellnum,] <- Heatmap_Matrix_Hp[cellnum,]
        Heatmap_Matrix_Hp_new[,cellnum] <- Heatmap_Matrix_Hp[,cellnum]
        Heatmap_Matrix_Fm_new[1,] <- Heatmap_Matrix_Fm[1,]
        Heatmap_Matrix_Fm_new[,1] <- Heatmap_Matrix_Fm[,1]
        Heatmap_Matrix_Fm_new[cellnum,] <- Heatmap_Matrix_Fm[cellnum,]
        Heatmap_Matrix_Fm_new[,cellnum] <- Heatmap_Matrix_Fm[,cellnum]
        Heatmap_Matrix_Fp_new[1,] <- Heatmap_Matrix_Fp[1,]
        Heatmap_Matrix_Fp_new[,1] <- Heatmap_Matrix_Fp[,1]
        Heatmap_Matrix_Fp_new[cellnum,] <- Heatmap_Matrix_Fp[cellnum,]
        Heatmap_Matrix_Fp_new[,cellnum] <- Heatmap_Matrix_Fp[,cellnum]
        Heatmap_Matrix_N_new <- as.matrix(Heatmap_Matrix_N_new)
        Heatmap_Matrix_Dp_new <- as.matrix(Heatmap_Matrix_Dp_new)
        Heatmap_Matrix_Dm_new <- as.matrix(Heatmap_Matrix_Dm_new)
        Heatmap_Matrix_Hm_new <- as.matrix(Heatmap_Matrix_Hm_new)
        Heatmap_Matrix_Hp_new <- as.matrix(Heatmap_Matrix_Hp_new)
        Heatmap_Matrix_Fm_new <- as.matrix(Heatmap_Matrix_Fm_new)
        Heatmap_Matrix_Fp_new <- as.matrix(Heatmap_Matrix_Fp_new)
        
    }
    return(Heatmap_Matrix_N_new)
}

iteration_notch_sim <- function(grpA, grpB) {
fwa_multi <- matrix(nrow=1, ncol=4)
fwa_multi <- data.frame(fwa_multi)
colnames(fwa_multi) <- c("fwa_A_ON","fwa_A_OFF", "fwa_B_ON", "fwa_B_OFF")
leftoutY <- sample_n(leftout_template40x40, 800)
fwa6 <- notch_model1_noise32_FINALnoplots(1000, 0.1, 0.01, 1, 0.01, 0.1, 0.1, 0.1, 0.1, 0.1, 0.01, 30, grpA, 0, 0, grpB, leftoutY)
fwa_B <- fwa6
for(k in 1:nrow(leftoutY)) {
    j <- leftoutY[k,1]
    i <- leftoutY[k,2]
    fwa_B[j,i] <- NA }
fwa_B[1,] <- NA
fwa_B[40,] <- NA
fwa_B[,1] <- NA
fwa_B[,40] <- NA
fwa_C <- melt(fwa_B)
fwa_C <- fwa_C[complete.cases(fwa_C), ]
fwa_B_ON <- nrow(subset(fwa_C, value > 0.5))
fwa_B_OFF <- nrow(subset(fwa_C, value < 0.5))
fwa_A <- matrix(nrow=nrow(leftoutY), ncol=3)
for(k in 1:nrow(leftoutY)) {
    j <- leftoutY[k,1]
    i <- leftoutY[k,2]
    fwa_A[k,1] <- j
    fwa_A[k,2] <- i
    fwa_A[k,3] <- fwa6[j,i] }
fwa_A <- data.frame(fwa_A)
fwa_A <- subset(fwa_A, X1 != 1 & X1 != 40 & X2 != 1 & X2 != 40)
fwa_A_ON <- nrow(subset(fwa_A, X3 > 0.5))
fwa_A_OFF <- nrow(subset(fwa_A, X3 < 0.5))
fwa_multi[1,1] <- fwa_A_ON
fwa_multi[1,2] <- fwa_A_OFF
fwa_multi[1,3] <- fwa_B_ON
fwa_multi[1,4] <- fwa_B_OFF
return(fwa_multi)
}

fwa_multi_wt_vs_wt <- iteration_notch_sim(2.2,1)
write.table(fwa_multi_wt_vs_wt, file = paste(Sys.time(), "fwa_1q_vs_wt.txt", sep = "_"), sep="\t", row.names=FALSE, col.names=TRUE, quote = FALSE)