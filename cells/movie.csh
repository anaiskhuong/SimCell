# after argv[1] seconds, grab argv[2] screen dumps every argv[3] seconds (or until "done" exists)
# the movie is saved in argv[4].avi and played back at argv[5] fps
rm done
rm *scrot.png
#sleep 60
sleep $argv[1]
echo Action
@ n = 0
if ( $argv[2] < 0 ) then
	@ end = 99999
else
	@ end = $argv[2]
endif
while ( $n < $end )
	if ( -e done ) break
	@ n++
	scrot -q 100 -u
	set snap = `ls -ltr *scrot.png | tail -1 | awk '{print $9}'`
	@ num = 10000 + $n
	mv $snap frame$num.png
	sleep $argv[3]
end
#mencoder "mf://*scrot.png" -mf type=png:fps=$argv[5] -ovc lavc -lavcopts vcodec=mpeg4:mbd=2:trell -o $argv[4].avi
#mencoder "mf://*scrot.png" -mf type=png:fps=$argv[5] -ovc lavc -o $argv[4].avi
#mplayer $argv[4].avi
