# liveRadar.jli
## Overview
liveRadar is a program to use with the MIT Coffee Can Radar (Candar).
This program has functions including:
- offline doppler waterfall
- online (by serial) doppler waterfall
- offline ranging waterfall 
- online ranging waterfall
## Function Calls
Calls are simple, use the module liveRadar and the function
```
	liveRadar.offlineDoppler(MATFILE)
	liveRadar.liveDoppler(PORT,BAUD,SAVE_FILE)
	liveRadar.offlineRanging(MATFILE)
	liveRadar.liveRanding(PORT,BAUD,SAVE_FILE)
```


