python ../../../spectraldensity/getS2.py -i rotacf-100ns-0000to1250ns.xvg > S2tol01-0000to1250ns.dat
python ../../../spectraldensity/getS2.py -i rotacf-100ns-1250to2500ns.xvg > S2tol01-1250to2500ns.dat
python ../../../spectraldensity/avgS2.py S2tol01-0000to1250ns.dat S2tol01-1250to2500ns.dat > S2tol01-2blocks.dat
