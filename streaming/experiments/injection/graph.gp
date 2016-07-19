set term png
FILES = system("ls -1 *.dat")
do for [data in FILES] {
	set output sprintf("%s.png",data)
	unset key
	set xlabel "Message Index"
	set ylabel "Message Gap"
	plot data u 0:($3<0.001 ? $3 : 0/1)
}
