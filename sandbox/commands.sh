#-----------------------------------------
# commands for scp'ing sce and cent over
#-----------------------------------------
cd ~/Dropbox/R/drinf/sandbox
scp  cent* sce* makeData.R dbenkese@snail.fhcrc.org:~/drinf

ssh dbenkese@snail.fhcrc.org
cd drinf
scp cent* sce* makeData.R dbenkese@rhino.fhcrc.org:~/drinf

ssh dbenkese@rhino.fhcrc.org
cd drinf
chmod +x cent* sce*
./sce_1000.sh ./cent_1000.R cvrun_1000_V3
./sce_5000.sh ./cent_5000.R cvrun_5000_V2
./sce_500.sh ./cent_500.R sub_test_2

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
