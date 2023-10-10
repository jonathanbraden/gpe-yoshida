mencoder mf://'*.png' -mf w=528:h=332:fps=30:type=png -ovc lavc -lavcopts vcodec=mpeg4:mbd=2:trell -oac copy -o output.avi
