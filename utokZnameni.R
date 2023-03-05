setwd("./")
library(epitools)
library(ggplot2)
library(plotly)
library(zoo)
library(reshape2)

## zabijime proklety do sun.k.int dnu po znameni
kill.em <- function(x, s.int=14, s.d=100) {
  min(x + as.integer(runif(1,0,s.int-1)),s.d)
}

get_cohorts <- function(scen=0, sun.p=1e6, sun.d=100, sun.d.pr=0.01, sun.s.pr=0.65,
                        sun.h.p=1, sun.h.int=14, sun.m.int=14, sun.k.p=1e-3, sun.k.int = 14) {

  sun.deaths <- sun.p * sun.d.pr
  sun.signs <- sun.p * sun.s.pr

  ## populacni dataframe
  # radek = 'ID' jedince
  # 1. sloupec - den umrti
  # 2. sloupec - den znameni slunce
  df <- data.frame(matrix(nrow = sun.p, ncol = 2))
  colnames(df) <- c("deaths","signs")

  ## dataframe z umrtim
  # radek = umrti jedince
  # 1. sloupec - den umrti
  # 2. sloupec - 'ID' jedince
  deaths <- data.frame(matrix(nrow = sun.deaths, ncol =2))
  colnames(deaths) <- c("day","person")

  ## dataframe se znamenim
  # radek = oznaceni jedince
  # 1. sloupec - den oznaceni
  # 2. sloupec - 'ID' jedince
  signs <- data.frame(matrix(nrow = sun.signs, ncol =2))
  colnames(signs) <- c("day","person")

  ## losujem umrlce [den]
  deaths$day <- as.integer(runif(sun.deaths,1,sun.d))
  ## losujem umrlce [ID]
  deaths$person <- sample(1:sun.p, sun.deaths)

  ## losujem slunickare [den]
  signs$day <- as.integer(runif(sun.signs,1,sun.d))
  ## losujem slunickare [ID]
  signs$person <- sample(1:sun.p,sun.signs)

  ## fouknem to do dataframe-u
  df$deaths[deaths$person] <- deaths$day
  df$signs[signs$person] <- signs$day

  ## zabijem par s prokletym znakem
  if (scen == 2 || scen == 3 || scen == 4) {
    # indexy prokletejch
    cursed <- sample(signs$person,sun.k.p*sun.signs)
    # umrou do 14 dni
    df$deaths[cursed] <- sapply(df[cursed,2],kill.em,s.int=sun.k.int,s.d=sun.d)
  }

  #dropnem zivy & nepoznamenany
  df.drp <- subset(df, deaths | signs>0)
  rownames(df.drp) <- seq(1,length(df.drp$deaths),1)

  ## dataframe pro statistiku
  # radek = den sledovani
  # 1. sloupec - pocet zemrelych bez znaku
  # 2. sloupec - pocet zemrelych se znakem
  # 3. sloupec - velikost kohorty bez znaku
  # 4. sloupec - velikost kohorty se znakem
  # 5. sloupec - pocet novych se znakem
  # 6. sloupec - riziko umrti bez znaku
  # 7. sloupec - riziko umrti se znakem
  df.st <- data.frame(matrix(0,nrow = sun.d, ncol = 7))
  colnames(df.st) <- c("nsdead","sdead","nonsig","sunsig","newsig","Rns","Rs")

  # pocet vyrazenych v scenari 1
  sc1.out <- 0
  # pocet slunci po smrti
  s.after <- 0

  for (i in 1:length(df.drp$deaths)) {
    # den umrti
    d.day <- df.drp$deaths[i]
    # den znameni
    s.day <- df.drp$signs[i]

    # ve scenari 3 posunem slunickare 14 dni dopredu
    if (scen == 3)
      if(!is.na(s.day))
        s.day <- min(s.day + sun.m.int, sun.d)

    # umrli bez znaku - prictem k NS umrlejm v danym dni
    if (!is.na(d.day) && is.na(s.day)) {
      df.st$nsdead[d.day] <- df.st$nsdead[d.day]+1
      next
    }

    # zivy se znakem
    if (is.na(d.day) && !is.na(s.day)) {
      df.st$newsig[s.day] <- df.st$newsig[s.day]+1
      next
    }

    # umrly se znakem
    if (!is.na(d.day) && !is.na(s.day)) {
      if (s.day > d.day) {
        df.st$nsdead[d.day] <- df.st$nsdead[d.day]+1
        s.after <- s.after + 1
        next
      }
      # poctiva klasifikace
      if (scen == 0 || scen == 2 || scen == 3) {
        df.st$sdead[d.day] <- df.st$sdead[d.day]+1
        df.st$newsig[s.day] <- df.st$newsig[s.day]+1
        next
      }
      # healthy vaccinee effect
      else if (scen == 1) {
        if ( (d.day-s.day) <= sun.h.int && runif(1) <= sun.h.p ) {
          df.st$nsdead[d.day] <- df.st$nsdead[d.day]+1
          sc1.out <- sc1.out + 1
          next
        }
        else {
          df.st$sdead[d.day] <- df.st$sdead[d.day]+1
          df.st$newsig[s.day] <- df.st$newsig[s.day]+1
          next
        }
      }
      # spatna klasifikace slunickaru -> slunickar mrtvej do 14 dnu je 'pro nas' bez slunicka
      else if (scen == 4) {
        if ( (d.day-s.day) <= sun.m.int ) {
          df.st$nsdead[d.day] <- df.st$nsdead[d.day] + 1
          df.st$newsig[s.day] <- df.st$newsig[s.day] + 1
          next
        } else {
          df.st$sdead[d.day] <- df.st$sdead[d.day] + 1
          df.st$newsig[s.day] <- df.st$newsig[s.day] + 1
          next
        }
      }
    }

  }

  #print(c("vyrazenych.healty_vac:",sc1.out,"vyrazenych.sun_after_death:",s.after))

  # nastavime prvni den + hned spoctem riziko
  df.st$nonsig[1] = sun.p - df.st$newsig[1] - df.st$nsdead[1]
  df.st$sunsig[1] = df.st$newsig[1] - df.st$sdead[1]

  df.st$Rns[1] <- df.st$nsdead[1]/df.st$nonsig[1]
  if (scen != 3) {
    df.st$Rs[1] <- df.st$sdead[1]/df.st$sunsig[1]
  }

  # korekce kohort bez/se znamenim + riziko
  for (i in 2:sun.d) {
    df.st$nonsig[i] = df.st$nonsig[i-1] - df.st$newsig[i] - df.st$nsdead[i]
    df.st$sunsig[i] = df.st$sunsig[i-1] + df.st$newsig[i] - df.st$sdead[i]
    df.st$Rns[i] <- df.st$nsdead[i]/df.st$nonsig[i]
    if (scen == 3 && i <= sun.m.int) next #prvnich par dni jsou tu 0 kvuli posunu
    df.st$Rs[i] <- df.st$sdead[i]/df.st$sunsig[i]
  }

  # zprumerujem klouzavym 7d prumerem at to neni tak strapaty
  df.st$Rns <- rollmean(df.st$Rns, k=7,fill=0, align="center")
  df.st$Rs <- rollmean(df.st$Rs, k=7,fill=0, align="center")

  return(df.st)
}

plot.risk.ens <- function(lens.s, lens.ns, save=FALSE, fn="") {

  ## Risks plot with ensembles
  mens.s <- melt(as.matrix(lens.s[7:(length(lens.s[,1])-7),]))
  mens.ns <- melt(as.matrix(lens.ns[7:(length(lens.ns[,1])-7),]))

  p <- ggplot() +
    ylab("Riziko úmrtí") + #ylim(c(0,1.25e-4)) +
    xlab("Den")  + #ylim(c(0,0.6e-4)) +
    ylim(c(0,1.75e-4)) +
    #scale_x_continuous(limits=c(7,365), breaks =c(7,seq(30, 350, by = 30),358)) +
    scale_x_continuous(limits=c(7,93), breaks =c(7,seq(10, 90, by = 10),93)) +
    #labs(title = "Riziko úmrtí") +
    geom_line(data = mens.s, aes(x = Var1, y = value, group = Var2, colour = "sun"), alpha = 0.3) +
    geom_line(data = mens.ns, aes(x = Var1, y = value, group = Var2, colour = "nosun"), alpha = 0.3) +
    scale_color_manual(name = "Skupina", limits=c("sun","nosun"), labels = c("se znamením","bez znamení"),
                       values = c("sun" = "red", "nosun" = "blue"))

  # osy
  p <- p + theme(axis.text=element_text(size=14),
                 axis.title=element_text(size=16,face="bold"))
  # legenda
  p <- p + theme(legend.key.size = unit(1, 'cm'), #change legend key size
                 legend.key.height = unit(1, 'cm'), #change legend key height
                 legend.key.width = unit(1, 'cm'), #change legend key width
                 legend.title = element_text(size=16), #change legend title font size
                 legend.text = element_text(size=14), #change legend text font size
                 legend.position = "top") #change legend position
  #legend.position = c(.9, .1)
  #p <- plotly_build(p)
  #p$x$data[[1]]$name <- "bez znamení"
  #p$x$data[[2]]$name <- "bez znamení"
  #p$x$data[[3]]$name <- "se znamením"
  #p$x$data[[4]]$name <- "se znamením"
  if(save) {
    ggsave(
      paste("Risk","_ensemble",fn,".png",sep=""),
      plot = p,
      device = "png",
      scale = 1,
      width = 2*1024,
      height = 2*768,
      units = "px")
  }
  p
}


## pocet dnu sledovani
sun.days <- 100
## pocatecni kohorta
sun.pop <- 1e6
## pocet 'all-cause' umrti
sun.deaths.pr <- 0.01
## pocet znameni
sun.signs.pr <- .65
## prst pro healthy vaccinee effect
sun.healthy.p <- 1
## interval[dny] pro healthy vaccinee effect
sun.healthy.int <- 14
## interval pro misklasifikaci
sun.misc.int <- 14
## prst ze slunce zabiji
sun.kill.p <- 1e-3
## behem tolika dnu po udeleni
sun.kill.int <- 14
## pocet opakovani 'modelu' pro ansambl
ens.c <- 100

## ansambly lidi se sluncem bez slunce

ens.s <- data.frame(matrix(0,nrow = sun.days, ncol = ens.c))
ens.ns <- data.frame(matrix(0,nrow = sun.days, ncol = ens.c))

# scenar 0 - znameni nezabiji, poctiva klasifikace, chodi vsichni nahodne
for (i in 1:ens.c) {
  stat <- get_cohorts(scen = 0, sun.p = sun.pop, sun.d = sun.days, sun.d.pr = sun.deaths.pr,
                      sun.s.pr = sun.signs.pr, sun.h.p=sun.healthy.p, sun.h.int=sun.healthy.int,
                      sun.m.int=sun.misc.int, sun.k.p=sun.kill.p, sun.k.int = sun.kill.int)

  ens.s[,i] <- stat$Rs
  ens.ns[,i] <- stat$Rns
  cat('*')
}
ens.s0 <- ens.s
ens.ns0 <- ens.ns
plot.risk.ens(ens.s0,ens.ns0)

# scenar 1 - lidi na umreni nechodi pro znameni, poctiva klasifikace
for (i in 1:ens.c) {
  stat <- get_cohorts(scen = 1, sun.p = sun.pop, sun.d = sun.days, sun.d.pr = sun.deaths.pr,
                      sun.s.pr = sun.signs.pr, sun.h.p=sun.healthy.p, sun.h.int=sun.healthy.int,
                      sun.m.int=sun.misc.int, sun.k.p=sun.kill.p, sun.k.int = sun.kill.int)

  ens.s[,i] <- stat$Rs
  ens.ns[,i] <- stat$Rns
  cat('*')
}
ens.s1 <- ens.s
ens.ns1 <- ens.ns
plot.risk.ens(ens.s1,ens.ns1)

# scenar 1.5 - lidi na umreni chodi pro znameni 50:50, poctiva klasifikace
for (i in 1:ens.c) {
  stat <- get_cohorts(scen = 1, sun.p = sun.pop, sun.d = sun.days, sun.d.pr = sun.deaths.pr,
                      sun.s.pr = sun.signs.pr, sun.h.p=sun.healthy.p/2, sun.h.int=sun.healthy.int,
                      sun.m.int=sun.misc.int, sun.k.p=sun.kill.p, sun.k.int = sun.kill.int)

  ens.s[,i] <- stat$Rs
  ens.ns[,i] <- stat$Rns
  cat('*')
}
ens.s1.5 <- ens.s
ens.ns1.5 <- ens.ns
plot.risk.ens(ens.s1.5,ens.ns1.5)


# scenar 2 - 1:1000 znameni zabiji, poctiva klasifikace
for (i in 1:ens.c) {
  stat <- get_cohorts(scen = 2, sun.p = sun.pop, sun.d = sun.days, sun.d.pr = sun.deaths.pr,
                      sun.s.pr = sun.signs.pr, sun.h.p=sun.healthy.p, sun.h.int=sun.healthy.int,
                      sun.m.int=sun.misc.int, sun.k.p=sun.kill.p, sun.k.int = sun.kill.int)

  ens.s[,i] <- stat$Rs
  ens.ns[,i] <- stat$Rns
  cat('*')
}
ens.s2 <- ens.s
ens.ns2 <- ens.ns
plot.risk.ens(ens.s2,ens.ns2)

# scenar 3 - 1:1000 znameni zabiji, znameni pocitame od 14+ dne
for (i in 1:ens.c) {
  stat <- get_cohorts(scen = 3, sun.p = sun.pop, sun.d = sun.days, sun.d.pr = sun.deaths.pr,
                      sun.s.pr = sun.signs.pr, sun.h.p=sun.healthy.p, sun.h.int=sun.healthy.int,
                      sun.m.int=sun.misc.int, sun.k.p=sun.kill.p, sun.k.int = sun.kill.int)

  ens.s[,i] <- stat$Rs
  ens.ns[,i] <- stat$Rns
  cat('*')
}
ens.s3 <- ens.s
ens.ns3 <- ens.ns
plot.risk.ens(ens.s3,ens.ns3)

# scenar 4 - 1:1000 znameni zabiji, pokud nekdo umre do 14 dnu od udeleni znameni, nepocita se mezi umrti se znamenim
for (i in 1:ens.c) {
  stat <- get_cohorts(scen = 4, sun.p = sun.pop, sun.d = sun.days, sun.d.pr = sun.deaths.pr,
                      sun.s.pr = sun.signs.pr, sun.h.p=sun.healthy.p, sun.h.int=sun.healthy.int,
                      sun.m.int=sun.misc.int, sun.k.p=sun.kill.p, sun.k.int = sun.kill.int)

  ens.s[,i] <- stat$Rs
  ens.ns[,i] <- stat$Rns
  cat('*')
}
ens.s4.1 <- ens.s
ens.ns4.1 <- ens.ns

save.image("./sun.RData")

plot.risk.ens(ens.s0,ens.ns0)
plot.risk.ens(ens.s1,ens.ns1)
plot.risk.ens(ens.s1.5,ens.ns1.5)
plot.risk.ens(ens.s2,ens.ns2)
plot.risk.ens(ens.s3,ens.ns3)
plot.risk.ens(ens.s4,ens.ns4)
