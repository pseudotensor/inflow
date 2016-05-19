pl	0	#
		da t00
		erase
		read {r 1 rho 2 ur 3 up 4 Br 5 Bp 6 Ma 7}
		#
		pla rho 1 2
		pla ur 2 2
		pla up 1 1
		pla Bp 2 1
		#
		window 1 1 1 1
pla	3	#
		window 2 2 $2 $3
		limits r $1
		box
		connect r $1
		xla r
		yla $1
