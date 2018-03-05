#-----------------------------------------
# commands for scp'ing sce and cent over
#-----------------------------------------
cd ~/Dropbox/R/drinf/sandbox
scp cent_all_norestart1.R sce_all_norestart.sh makeData.R dbenkese@snail.fhcrc.org:~/drinf

ssh dbenkese@snail.fhcrc.org
cd drinf
scp cent_all_norestart1.R sce_all_norestart.sh makeData.R dbenkese@rhino.fhcrc.org:~/drinf

ssh dbenkese@rhino.fhcrc.org
cd drinf
chmod +x cent* sce*
# ./sce_9000.sh ./cent_9000.R fix_9000_v1
# ./sce_1000.sh ./cent_1000.R fix_1000_v1
# ./sce_5000.sh ./cent_5000.R fix_5000_v1
# ./sce_500.sh ./cent_500.R fix_500_v1
# ./sce_all.sh ./cent_all.R noboot_v4
./sce_all_norestart.sh ./cent_all_norestart1.R run_6

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
