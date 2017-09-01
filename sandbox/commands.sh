#-----------------------------------------
# commands for scp'ing sce and cent over
#-----------------------------------------
cd ~/Dropbox/R/drinf/sandbox
scp  cent.R sce.sh makeData.R dbenkese@snail.fhcrc.org:~/drinf

ssh dbenkese@snail.fhcrc.org
cd drinf
scp cent.R sce.sh makeData.R dbenkese@rhino.fhcrc.org:~/drinf

ssh dbenkese@rhino.fhcrc.org
cd drinf
chmod +x cent.R sce.sh 
./sce.sh ./cent.R run_V1

#-----------------------------------------
# commands to get into rhino and load R
#-----------------------------------------
ssh dbenkese@snail.fhcrc.org
ssh dbenkese@rhino.fhcrc.org
 # enter password
ml R/3.2.0
R
# module avail
#-----------------------------------------
# scp results from rhino to local machine
#-----------------------------------------
# from rhino
cd drinf/out
scp allOut.RData dbenkese@snail.fhcrc.org:~/drinf
	# enter snail password
 	# ctrl + shift + t to open up new term
# scp to snail
cd ~/Dropbox/Emory/drinfSieve/simulation
scp dbenkese@snail.fhcrc.org:~/drinf/allOut.RData . 

#-----------------------------------------
# misc commands 
#-----------------------------------------
squeue -u dbenkese
# scancel `seq 51239645 51239655`
