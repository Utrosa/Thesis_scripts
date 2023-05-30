import numpy as np 
from time import sleep

sbj = np.arange(1000)
length_bar = 10
bar = '-'*length_bar
replace_no = len(sbj)/length_bar

for i in sbj:
	num_squares = int(i/replace_no)
	bar = '#'*num_squares+'-'*(length_bar-num_squares)
	percentage = round(i/len(sbj) * 100)
	sleep(.1)
	print(bar, percentage, '%', end='\r')



	# print(f'{i}/{sbj.shape[0]}', end='\r')
	#sleep(.1)

