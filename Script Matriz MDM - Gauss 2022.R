##### setup  ######
{
  #setwd("C:/Users/capel/Dropbox/TCC/Genetic Algorithms") #Setar diretório de trabalho
  library("readxl")
  library("dplyr")
  library("GA")
  library("rgenoud")
}

{
  mt <- read_excel("Matriz MDM - Gauss 2022.xlsx")
  mt <- mt %>% select(-1)
  mt <- as.matrix(mt)
  dimnames(mt) <- (list(c(1:ncol(mt)),c(1:ncol(mt))))
  nm <- 7 #definir número de módulos previamente 
  #mt
  #ncol(mt)
}

#x <- c(0,0,0,0,0,0,0,0,0,0,16,15,14,9,5,3,0) 
##### Cromossomos binário para decimal #####

#Vetor a ser gerado no GA
#o projeto não aceita ser calculado com apenas 1 módulo
vga <- as.integer(c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,
                    1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,
                    1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,0,
                    1,1,1,1,1,1,1,0,0,1,1,0,0,0,0,0,
                    1,1,1,1,0,1,0,0,0,0,0,0,0,0,0,0,
                    0,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,
                    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0))

vga

#Cromossomo x final 
x <- as.integer(c(sum(vga[1:16]), sum(vga[17:32]), sum(vga[33:48]), sum(vga[49:64]), sum(vga[65:80]), 
                  sum(vga[81:96]), sum(vga[97:112])))
x
length(x)


############### Cálculo de MI ###############
MI <- function(vga){
  ##### Parametros iniciais ####
  offset <- 0.12 #offset para penalizações
  
  x <- as.integer(c(sum(vga[1:16]), sum(vga[17:32]), sum(vga[33:48]), sum(vga[49:64]), sum(vga[65:80]), 
                    sum(vga[81:96]), sum(vga[97:112])))
  nm <- 7
  nc <- ncol(mt)
  BC <- 0
  WC <- 0
  WCi1 <- 0
  SWC <- WC #soma do WC
  WCj1 <- 0
  BS <- 0
  WS <- 0 
  WSi1 <- 0
  WSj1 <- 0 
  ST <- 0
  q <- 0 
  p <- 0 
  qp <- 0 
  SR1 <- 0
  SR2 <- 0
  
  ##### MI1 #####
  for ( i in 1:(nm-1)){
    j <- i+1
    if (i == (nm-1)){
      WCi1 <- {
        sum((mt[(x[nm-(i-1)]+1):(x[nm-(j-1)]) , (x[nm-(i-1)]+1): (x[nm-(j-1)])])!=0)
      }
      BC <- {
        sum((mt[(x[nm-(j-1)]+1):ncol(mt),(x[nm-(i-1)]+1):x[nm-(j-1)]])!=0) +  
          sum((mt[(x[nm-(i-1)]+1):x[nm-(j-1)],(x[nm-(j-1)]+1):ncol(mt)])!=0)  
      }
      WCj1 <- {
        sum((mt[(x[nm-(i)]+1):(ncol(mt)) , (x[nm-(i)]+1):(ncol(mt))])!=0)
      }
      WSi1 <- {
        sum(mt[(x[nm-(i-1)]+1):(x[nm-(j-1)]) , (x[nm-(i-1)]+1): (x[nm-(j-1)])])
      }
      BS <- {
        sum(mt[(x[nm-(j-1)]+1):ncol(mt),(x[nm-(i-1)]+1):x[nm-(j-1)]]) +  
          sum(mt[(x[nm-(i-1)]+1):x[nm-(j-1)],(x[nm-(j-1)]+1):ncol(mt)])  
      }
      WSj1 <- {
        sum(mt[(x[nm-(i)]+1):(ncol(mt)) , (x[nm-(i)]+1):(ncol(mt))])
      }
      
      ST <- ST + ((BC/WCi1)+(BC/WCj1)+(BS/WSi1)+(BS/WSj1))
      
      break
      
    } else {
      
      WCi1 <- {
        sum((mt[(x[nm-(i-1)]+1):(x[nm-(j-1)]) , (x[nm-(i-1)]+1): (x[nm-(j-1)])])!=0)
      }
      BC <- {
        sum((mt[(x[nm-(j-1)]+1):ncol(mt),(x[nm-(i-1)]+1):x[nm-(j-1)]])!=0) +  
          sum((mt[(x[nm-(i-1)]+1):x[nm-(j-1)],(x[nm-(j-1)]+1):ncol(mt)])!=0)  
      }
      WCj1 <- {
        sum((mt[(x[nm-(i)]+1):(x[nm-(j)]) , (x[nm-(i)]+1):(x[nm-(j)])])!=0)
      }
      WSi1 <- {
        sum(mt[(x[nm-(i-1)]+1):(x[nm-(j-1)]) , (x[nm-(i-1)]+1): (x[nm-(j-1)])])
      }
      BS <- {
        sum(mt[(x[nm-(j-1)]+1):ncol(mt),(x[nm-(i-1)]+1):x[nm-(j-1)]]) +  
          sum(mt[(x[nm-(i-1)]+1):x[nm-(j-1)],(x[nm-(j-1)]+1):ncol(mt)])  
      }
      WSj1 <- {
        sum(mt[(x[nm-(i)]+1):(x[nm-(j)]) , (x[nm-(i)]+1):(x[nm-(j)])])
      }
      
      ST <- ST + ((BC/WCi1)+(BC/WCj1)+(BS/WSi1)+(BS/WSj1))
      
    }
  }
  
  MI1 <- 1-(ST/(2*nm*(nm-1)))
  MI1
  ##### MI2 #####
  for (i in 1:nm){
    j <- i+1
    
    if (i == (nm)){
      WC <- sum((mt[(x[nm-(i-1)]+1):nc , (x[nm-(i-1)]+1):nc])!=0)
      SWC <- SWC + WC
      q <- nc
      p <- (x[nm-(i-1)]+1)
      qp <- qp + (q-p+1)^2
      
      break
    } else{  
      WC <- {
        sum((mt[(x[nm-(i-1)]+1):(x[nm-(j-1)]) , (x[nm-(i-1)]+1): (x[nm-(j-1)])])!=0)
      }
      q <- x[nm-(j-1)]
      p <- (x[nm-(i-1)]+1)
      qp <- qp + (q-p+1)^2
      SWC <- SWC + WC
      
    }
  }
  
  MI2 <- (SWC + nc)/qp
  
  
  ##### MI3 #####
  for (i in 1:(ncol(mt)-1)){
    j <- i+1
    
    SR1 <- SR1 + ((1-(j-i)/(ncol(mt)-1))*(mt[i,j] +mt[j,i])/max(mt))
    SR2 <- SR2 + ((mt[i,j] +mt[j,i])/max(mt))
  }
  MI3 <- SR1/SR2
  
  ##### Índice de modularidade MI #####
  
  MI <- (MI1 * (1/3)) + (MI2 * (1/3)) + (MI3 * (1/3))  
  
  ##### Penalizações #####
  
  #Ultima variável de X precisa ser 0
  
  P1 <- ifelse(x[nm] != 0, offset*1.5, 0)  
  
  #Resultados de MI1, MI2, MI3 e MI devem ser entre 0 e 1
  
  P2 <- ifelse(MI1 <= 0, offset*3, 0)
  P3 <- ifelse(MI2 <= 0, offset*3, 0)
  P4 <- ifelse(MI3 <= 0, offset*3, 0)
  P5 <- ifelse(MI1 > 1, offset*3, 0)
  P6 <- ifelse(MI2 > 1, offset*3, 0)
  P7 <- ifelse(MI3 > 1, offset*3, 0)
  
  #As veriaveis de x devem ser em ordem crescente da direita para a esquerda
  
  P81 <- ifelse(x[nm-1] <= x[nm], offset * (1/3), 0 )
  P82 <- ifelse(x[nm-2] <= x[nm-1], offset * (1/3),0 )
  P83 <- ifelse(x[nm-3] <= x[nm-2], offset * (1/3),0 )
  P84 <- ifelse(x[nm-4] <= x[nm-3], offset * (1/3),0 )
  P85 <- ifelse(x[nm-5] <= x[nm-4], offset * (1/3),0 )
  P86 <- ifelse(x[nm-6] <= x[nm-5], offset * (1/3),0 )
  
  P8 <- (P81+P82+P83+P84+P85+P86)
  
  
  ##### Retorno da função #####
  #print(P1)
  #print(P2)
  #print(P3)
  #print(P4)
  #print(P5)
  #print(P6)
  #print(P7)
  #print(P8)
  #print(MI1)
  #print(MI2)
  #print(MI3)
  return (MI - P1 - P2 - P3 - P4 - P5 - P6 - P7 - P8)
  #return (MI)
}

MI(vga)

##### GA #####


GA <- ga(type = "binary", fitness = MI, nBits = ((ncol(mt)-1)*nm), upper = max, maxiter = 200, 
         popSize = 300, pcrossover = 0.8, pmutation = 0.1, seed = 105)

plot(GA)
summary(GA)
GA
vga <- as.integer(GA@solution[1,])


#Cromossomo x final 
x <- as.integer(c(sum(vga[1:16]), sum(vga[17:32]), sum(vga[33:48]), sum(vga[49:64]), sum(vga[65:80]), 
                  sum(vga[81:96]), sum(vga[97:112])))
x
vga
