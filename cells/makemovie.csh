# 1 = run file, 2 = delay (60)
./simcell $argv[1].run $argv[2] &
tcsh movie.csh $argv[3] -1 2
avconv -r 15 -i "frame1%04d.png" -b:v 20000k $argv[1].mp4
mplayer $argv[1].mp4
