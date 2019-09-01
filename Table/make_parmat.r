parmat <- data.frame(NPar = 0, Animal = c('Multimammate_Rat',
                                          'Raccoon','Raccoon',
                                          'Fox', 'Fox',
                                          'Possum', 'Possum',
                                          'Badger',
                                          'Deer', 'Deer',
                                          'Elk',
                                          'Pig', 'Pig',
                                          'Sheep', 'Sheep',
                                          'Buffalo', 'Buffalo',
                                          'Bison', 'Bison',
                                          'Fruit Bat'
                                          ), tb = NA, LS = NA, d = NA,
			      gamv = 1/14, Nv = 500, T = 365, NPeak = 1000, tvstar = 0,
                     	      downtv90 = 0, uptv90 = 0,
                     	      downtv95 = 0, uptv95 = 0,
		     	      vavgstar = 0, vnull = 0)

wi <- parmat$Animal=='Multimammate_Rat'
parmat[wi,'tb'] <- 4*30.4
parmat[wi, 'LS'] <- 1*365

wi <- parmat$Animal=='Raccoon'
parmat[wi,'tb'] <- c(1.5,5.5)*30.4
parmat[wi, 'LS'] <- 2.5*365

wi <- parmat$Animal=='Fox'
parmat[wi,'tb'] <- c(2,3)*30.4
parmat[wi, 'LS'] <- 2.5*365

wi <- parmat$Animal=='Possum'
parmat[wi,'tb'] <- c(2,4)*30.4
parmat[wi, 'LS'] <- 6.5*365

wi <- parmat$Animal=='Badger'
parmat[wi,'tb'] <- 2*30.4
parmat[wi, 'LS'] <- 4*365

wi <- parmat$Animal=='Deer'
parmat[wi,'tb'] <- c(0.5,5)*30.4
parmat[wi, 'LS'] <- 5*365

wi <- parmat$Animal=='Elk'
parmat[wi,'tb'] <- 2*30.4
parmat[wi, 'LS'] <- 10*365

wi <- parmat$Animal=='Pig'
parmat[wi,'tb'] <- c(5,12)*30.4
parmat[wi, 'LS'] <- 10*365

wi <- parmat$Animal=='Sheep'
parmat[wi,'tb'] <- c(2,3)*30.4
parmat[wi, 'LS'] <- 6.5*365

wi <- parmat$Animal=='Buffalo'
parmat[wi,'tb'] <- c(4,12)*30.4
parmat[wi, 'LS'] <- 18*365

wi <- parmat$Animal=='Bison'
parmat[wi,'tb'] <- c(1,4)*30.4
parmat[wi, 'LS'] <- 10*365

wi <- parmat$Animal=='Fruit Bat'
parmat[wi,'tb'] <- c(2)*30.4
parmat[wi, 'LS'] <- 10*365


## Calculate rate of mortality
parmat$d <- 1/parmat$LS

## Remove duplicated rows
parmat_unique <- unique(parmat)
parmat <- parmat_unique

