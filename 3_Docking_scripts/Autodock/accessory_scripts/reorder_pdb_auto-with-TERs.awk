BEGIN{i=0;j=0}
{
    if ($1 ~ /CRYST1/) {print $0;i=-1;j=0;next}
    if ($1 ~ /^TER/) {
	j=0
	print "TER"
    } else if ($1 ~ /^END/) {
	print "END"
    } else {
        if(substr($0,23,4) != prev){
	    prev=substr($0,23,4);j++
	    print "TER"
	}
	printf"%s",substr($0,0,4)
	printf "%7d",NR+i
	printf"%s",substr($0,12,10)
	printf "%5d",j
	printf"%-s\n",substr($0,27)
    }
}
