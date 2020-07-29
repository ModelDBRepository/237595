
function la_run {
	echo $LAPARAMS

	#qsub -v "LAPARAMS=$LAPARAMS" submit_lamodel.sh
	./lamodel $LAPARAMS 
}


for nltype in 0 1 2 3; do
	for clustered in 0 1; do
		for run in {0..20}  ; do

			# sublinear SOM

			LAPARAMS=" -P 1 -T 180 -S 1980$run -o nlTypePV=${nltype} -o nlTypeSOM=1 -o INClustered=${clustered} -s single_${nltype}_1_${clustered}_G_${run} -G"
			la_run

			# supralinear SOM
			LAPARAMS=" -P 1 -T 180 -S 1980$run -o nlTypePV=${nltype} -o nlTypeSOM=0 -o INClustered=${clustered} -s single_${nltype}_0_${clustered}_G_${run} -G"
			la_run

			# linear SOM
			LAPARAMS=" -P 1 -T 180 -S 1980$run -o nlTypePV=${nltype} -o nlTypeSOM=2 -o INClustered=${clustered} -s single_${nltype}_2_${clustered}_G_${run} -G"
			la_run

			# mixed SOM
			LAPARAMS=" -P 1 -T 180 -S 1980$run -o nlTypePV=${nltype} -o nlTypeSOM=3 -o INClustered=${clustered} -s single_${nltype}_3_${clustered}_G_${run} -G"
			la_run

		done
	done
done


for nltype in 0 1 2 3; do
	for clustered in 0 1; do
		for run in {0..9}  ; do


			LAPARAMS=" -P 2 -T 1440 -S 1980$run -o nlTypePV=${nltype}  -o nlTypeSOM=1 -o INClustered=${clustered} -s two1440_${nltype}_${clustered}_G_${run} -G"
			la_run

			LAPARAMS=" -P 2 -T 60 -S 1980$run -o nlTypePV=${nltype} -o nlTypeSOM=1 -o INClustered=${clustered} -s two_${nltype}_${clustered}_G_${run} -G"
			la_run


		done
	done
done


