start=1
end=9
len=$(cat kmer.fa)
len=${#len}
numk=$(($len - 9 + 1))

for r in $(seq $start $numk); do
	k=$(cut -c$start-$end kmer.fa)

	echo ========================
	echo current kmer: $k
	echo ========================
	
	echo $k >> kmers.txt
	((start++))
	((end++))
done