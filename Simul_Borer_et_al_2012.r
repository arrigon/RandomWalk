### Random walk + diffusion models,
## Plant - insect interactions
## In short, this script simulates the colonization of a landscape by two lineages of plants and insects,
## starting from refugial cells located at each extremity of the landscape
## the idea is simply to show that co-dispersing organisms will share their biogeographic patterns,
## until the complete landscape is colonized. This pattern then vanishes unless co-dispersal is mandatory.
## Check out Borer et al. 2012 for the complete details.

possiblepos=function(dd, cell){
    if(dd == 0){
      lpos = c(cell$x, cell$y)
    } else {
      lgt = 2*dd+1
      
      xincr = rep((cell$x-dd):(cell$x+dd), lgt)
      yincr = rep((cell$y-dd):(cell$y+dd), each = lgt)

      lpos = cbind(xincr, yincr)
      dists = as.matrix(dist(lpos, "manhattan"))
      test = dists[1+(nrow(dists)-1)/2, ]
      lpos = lpos[test<=dd, ]
      rmv = which(lpos[, 1]==cell$x & lpos[, 2]==cell$y)
      lpos = lpos[-rmv, ]
      }
    lpos
    }



####### Plant dispersal
#######################
disperse = function(oldmap, tx.migr, abil, flyers, col, ...){
  maplims = rbind(range(oldmap$x), range(oldmap$y))
  newmap = oldmap
  newmap$remain = newmap$nmax - rowSums(newmap[, c(3,4)])

  #### Step 1: Population growth 
  newpop = newmap[, col] + round((newmap[,]$nsurv*newmap[, col])/2)
  newpop[newpop > newmap$remain] = newmap[newpop > newmap$remain, col] + newmap[newpop > newmap$remain, ]$remain
  newmap[, col] = newpop
  newmap$remain = newmap$nmax - rowSums(newmap[, c(3,4)])

  #### Step2: Migration 
  ## A. get cells that can disperse
  if(sum(newmap[, col], na.rm=T) > 0){
    allcells = newmap[newmap[, col] > 3, ] # at least 3 indivs for being able to disperse
  
    ## B. make them dispersing
    for(j in 1:nrow(allcells)){
      
      # cell specific stats
      cell = allcells[j, ]

      # Make migration
      migrevents = rbinom(1, cell[, col], tx.migr) # compute number of possible destinations
      radius = rpois(migrevents, abil) # compute dispersal distances
      nmigr = rpois(migrevents, flyers) # compute number of migrants
      while(sum(nmigr) > cell[, col]){ # check that not more than pop size will leave
	nmigr = nmigr - 1
	nmigr[nmigr < 1] = 1
	}
      outmigr = sum(nmigr) # compute nb of guy that effectively leave
      
      ## C. substract migrants from source cell
      cell[, col]=cell[, col] - outmigr # adjust pop size substract migrants to pop size
      newmap[newmap$x == cell$x & newmap$y == cell$y, col] = cell[,col]      
      newmap$remain = newmap$nmax - rowSums(newmap[, c(3,4)])

      ## D. Finalize migration to target cells
      if(migrevents>0){
	for(i in 1:length(radius)){
	  # list all reachable cells, from emigration origin
	  newloc=possiblepos(radius[i],cell) #select randomly "migrevens" cells

	  #  select "migrevents" cells, at random from this list
	  if(is.matrix(newloc)) newloc = newloc[sample(1:nrow(newloc),1), ]

	  # make sure migrants will not move out of map (kill them if out moving)
	  if(newloc[1] >= maplims[1, 1] & newloc[1] <= maplims[1, 2] & newloc[2] >= maplims[2, 1] & newloc[2]<=maplims[2, 2]){

	    # increment target cell with number of migrants
	    poptarget = nmigr[i] + newmap[newmap$x == newloc[1] & newmap$y == newloc[2], col]

	      # make sure that migrants will not overload cell capacity (kill them if overload)
	      if(poptarget < newmap[newmap$x == newloc[1] & newmap$y == newloc[2], ]$remain){
		newmap[newmap$x==newloc[1]&newmap$y==newloc[2],col] = poptarget
		} else {
		newmap[newmap$x == newloc[1] & newmap$y == newloc[2], col] = newmap[newmap$x == newloc[1] & newmap$y == newloc[2], col] + newmap[newmap$x == newloc[1] & newmap$y == newloc[2], ]$remain
	      }
	    }
	  }
	}
      }

      #### Step 3. Deaths (as a function of nsurv)
      # cells cleaning (remove NA or negative values, if any; might happen due to overkills)
      newmap[is.na(newmap[, col]), col] = 0
      newmap[newmap[, col] < 0, col] = 0 #just to make sure

      # Compute number of dying dudes, note that nsurv (0 or 1) kills or keeps the complete cell alive
      # This scheme shuts down the complete cell IF nsurv = 0, 
      # because the number of kills starts right away with the total number of dudes in the cell...
      kills = round(newmap[, col] * (1 - newmap$nsurv), 1) # BANG!
      kills = kills + rbinom(nrow(newmap), newmap[, col], mort)

      # update pop sizes
      newmap[, col] = newmap[, col] - kills
      newmap[newmap[, col] < 0, col] = 0
      }

      #### Step 4. Update map
      newmap$remain = newmap$nmax - rowSums(newmap[, c(3,4)])
      assign("oldmap", newmap, env=.GlobalEnv)
      }



###### Insect dispersal
###### Conditional on plant
###########################
disperseI = function(oldmap, oldmapI, tx.migrI, abilI, flyersI, col, ...){
  maplims = rbind(range(oldmap$x), range(oldmap$y))

  # newmap: based on insect map
  newmap = oldmapI

  # compute number of plants already available as hosts)
  nbplants = c(oldmap$A + oldmap$B)

  # convert background plant as 0/1, for determining insect survival chances
  # link insect survival to plant presence / absence
  nbplants[nbplants >=1 ] = 1
  nbplants[nbplants < 0] = 0
  newmap$nsurv = nbplants

  newmap$remain = newmap$nmax - rowSums(newmap[,c(3,4)])


  #### Step 1 Population growth 
  newpop = newmap[, col] + round((newmap[,]$nsurv * newmap[, col])/2)
  newpop[newpop > newmap$remain] = newmap[newpop > newmap$remain, col] + newmap[newpop > newmap$remain, ]$remain
  newmap[, col] = newpop
  newmap$remain = newmap$nmax - rowSums(newmap[, c(3,4)])

  #### Step 2 Migration
  if(sum(newmap[, col], na.rm=T) > 0){

    # get cells that can disperse
    allcells = newmap[newmap[,col] > 3,] # at least 3 indivs for being able to disperse
  
    ## make them dispersing
    for(j in 1:nrow(allcells)){
      
      # cell specific stats
      cell=allcells[j, ]

      # Make migration
      migrevents = rbinom(1, cell[, col], tx.migrI) # calcule nb destinations possibles
      radius = rpois(migrevents, abilI) # 
      nmigr = rpois(migrevents, flyersI) # calcule le nb de zozos prêts à partir
      while(sum(nmigr)>cell[, col]){ # check that not more than effectif pop will leave
	nmigr = nmigr - 1
	nmigr[nmigr < 1] = 1
	}
      outmigr = sum(nmigr) # compute nb of guy that effectively leave
      cell[, col] = cell[, col] - outmigr
      
      ## substract migrants from source cell
      newmap[newmap$x == cell$x & newmap$y == cell$y, col]=cell[, col]      
      newmap$remain = newmap$nmax - rowSums(newmap[, c(3,4)])

      ## Make migration to target cells
      # list reachable cells, select length(radius) at random from this list
      if(migrevents>0){
	for(i in 1:length(radius)){
	  newloc=possiblepos(radius[i],cell)
	  if(is.matrix(newloc)) newloc = newloc[sample(1:nrow(newloc), 1), ]
	  if(newloc[1] >= maplims[1, 1] & newloc[1] <= maplims[1, 2] & newloc[2] >= maplims[2, 1] & newloc[2] <= maplims[2, 2]){
	    poptarget = nmigr[i] + newmap[newmap$x == newloc[1] & newmap$y == newloc[2], col]
	      if(poptarget < newmap[newmap$x == newloc[1] & newmap$y == newloc[2], ]$remain){
		newmap[newmap$x == newloc[1] & newmap$y == newloc[2], col] = poptarget
		} else {
		newmap[newmap$x == newloc[1] & newmap$y == newloc[2], col] = newmap[newmap$x == newloc[1] & newmap$y == newloc[2], col] + newmap[newmap$x == newloc[1] & newmap$y == newloc[2], ]$remain
	      }
	    }
	  }
	}
      }

      ## Make kills (as a function of nsurv)
      newmap[is.na(newmap[, col]), col] = 0
      newmap[newmap[, col] < 0, col] = 0

      # This scheme shuts down the complete cell IF nsurv = 0, 
      # because the number of kills starts right away with the total number of dudes in the cell...
      kills = round(newmap[, col] * (1 - newmap$nsurv), 1) # BANG!
      kills = kills + rbinom(nrow(newmap), newmap[, col], mort)

      newmap[, col] = newmap[, col] - kills
      newmap[newmap[, col] < 0, col] = 0
      }
      ## Update map
      newmap$remain = newmap$nmax - rowSums(newmap[, c(3,4)])
      assign("oldmapI", newmap,env=.GlobalEnv)
      }


############### SIMULATION
# slow x rapid: tx.migr = 0.05 x 0.1 abil = 1 x 2 flyers = 2 x 2 mort = 0.05 x 0.05
# slow x slow: tx.migr = 0.05 x 0.05 abil = 1 x 1 flyers = 2 x 2 mort = 0.05 x 0.05
# rapid x rapid: tx.migr = 0.1 x 0.05 abil = 2 x 1 flyers = 2 x 2 mort = 0.05 x 0.05
# rapid x slow: tx.migr = 0.1 x 0.1 abil = 2 x 2 flyers = 2 x 2 mort = 0.05 x 0.05

## Init
## Nb destinations / migration, rbinom(ntrials=nb indivs / cell, proba succès=tx.migr)
tx.migr=.1
tx.migrI=.1

## radius / migration (lambda poisson)
abil=2
abilI=2

## Effectif de migrants (lambda poisson)
flyers=2
flyersI=2

## Mortalité / generation, rbinom(ntrials=nb indivs / cell, proba succès=mort)
mort=0.05
mortI=0.05

OUT=NULL
for(rep in 5:32){
  origA=c(13,1)
  origB=c(13,25)

  oldmap=data.frame(x=rep(1:25,25),
		    y=rep(1:25,each=25),
		    A=0,
		    B=0,
		    nmax=20,
		    nsurv=1) 

  oldmap[oldmap$x==origA[1]&oldmap$y==origA[2],3]=15
  oldmap[oldmap$x==origB[1]&oldmap$y==origB[2],4]=15

  oldmapI=data.frame(x=rep(1:25,25),
		    y=rep(1:25,each=25),
		    A=0,
		    B=0,
		    nmax=20,
		    nsurv=1) 

  oldmapI[oldmapI$x==origA[1]&oldmapI$y==origA[2],3]=15
  oldmapI[oldmapI$x==origB[1]&oldmapI$y==origB[2],4]=15


  ## Start Simulations
  gen=0
  pops=NULL
  corrs=NULL
  glob_pop=NULL
  sht='off'
  count=0
  grph=0
  while(gen<100){
    disperse(oldmap,tx.migr,abil,flyers,col=3,mort)
    disperse(oldmap,tx.migr,abil,flyers,col=4,mort)
    disperseI(oldmap,oldmapI,tx.migrI,abilI,flyersI,col=4,mortI)
    disperseI(oldmap,oldmapI,tx.migrI,abilI,flyersI,col=3,mortI)
    pops=rbind(pops,c(colSums(oldmap[,3:4]),colSums(oldmapI[,3:4])))
   
    glob_pop=sum(oldmapI[,3:4])

    if(glob_pop>11000) sht='on'
    if(sht=='on') gen=gen+1
    count=count+1
      if(gen==25){
     OUT=data.frame(rep,station=1:nrow(oldmap),count,x=oldmap$x,y=oldmap$y,A=oldmap$A,B=oldmap$B,AI=oldmapI$A,BI=oldmapI$B)
	write.table(OUT,paste("SIM25_rep",rep
		    ,"_tx",tx.migr,"x",tx.migrI
		    ,"_ab",abil,"x",abilI
		    ,"_fl",flyers,"x",flyersI
		    ,"_mt",mort,"x",mortI,".txt",sep=''),sep='\t')
      }
      
      if(gen==50){
        OUT=data.frame(rep,station=1:nrow(oldmap),count,x=oldmap$x,y=oldmap$y,A=oldmap$A,B=oldmap$B,AI=oldmapI$A,BI=oldmapI$B)
	write.table(OUT,paste("SIM50_rep",rep
		    ,"_tx",tx.migr,"x",tx.migrI
		    ,"_ab",abil,"x",abilI
		    ,"_fl",flyers,"x",flyersI
		    ,"_mt",mort,"x",mortI,".txt",sep=''),sep='\t')
      }
    }

  final=data.frame(rep,station=1:nrow(oldmap),count,x=oldmap$x,y=oldmap$y,A=oldmap$A,B=oldmap$B,AI=oldmapI$A,BI=oldmapI$B)
  write.table(final,paste("SIM100_rep",rep
	      ,"_tx",tx.migr,"x",tx.migrI
	      ,"_ab",abil,"x",abilI
	      ,"_fl",flyers,"x",flyersI
	      ,"_mt",mort,"x",mortI,".txt",sep=''),sep='\t')
  write.table(pops,paste("POP_rep",rep
	      ,"tx",tx.migr,"x",tx.migrI
	      ,"_ab",abil,"x",abilI
	      ,"_fl",flyers,"x",flyersI
	      ,"_mt",mort,"x",mortI,".txt",sep=''),sep='\t')
  }

