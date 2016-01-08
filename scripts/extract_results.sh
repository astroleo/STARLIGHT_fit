##
## bash helper script to extract synthetic (best-fit) spectrum from STARLIGHT output file
##
sldir="$HOME/STARLIGHT"

i=0
for outfile in `ls $sldir/out_dir/*.out`; do
	##
	## uncomment next line to only run the routine for first source
#	if [ $i -ge 1 ]; then continue; fi

	file_synspec="${outfile}.synspec"
#	sed -n '217,$p' < $outfile > $file_synspec
	sed -n '532,$p' < $outfile > $file_synspec

	file_popvec="${outfile}.popvec"
#	sed -n '64,108p' < $outfile > $file_popvec
	sed -n '64,213p' < $outfile > $file_popvec
	
	file_mask="$sldir/masks/$(grep 'arq_masks' $outfile | awk '{print $1}')"
	file_mask_out="${outfile}.masks"
#	echo $file_mask
	nmasks=$(cat $file_mask | head -n 1)
	nmasks=$((nmasks + 1))
	sed -n "2,${nmasks}p" < $file_mask | awk '{print $1 " " $2}' > $file_mask_out

	i=$((i+1))
done
