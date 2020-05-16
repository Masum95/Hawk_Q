files=$(find . -type d | cut -c 3-)

for file in $files
do
	found=0
	for fl in `cat links.txt`
	do
		if [ $fl == $file ]; then 
			found=1
			break
		fi
		
	done;
	if [ $found == 0 ]; then
		rm -rf $file 
		echo $file >> 'deletedFiles'
		echo {${file},'not found'}
	fi
done;
