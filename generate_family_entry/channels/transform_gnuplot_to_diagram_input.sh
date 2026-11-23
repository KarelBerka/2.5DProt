#/usr/bin/bash
for i in `ls generateds`; do 
    echo "---1"  $i
    j=`echo $i | sed "s|generated-||" | sed "s|-.*||"`;
    echo "---2" $i $j;
    grep -n "Id" generateds/$i/annotation/*channel.json | sed "s|generateds/$i/annotation/||" | sed 's|-channel.json:.*Id":"| |' | sed 's|",||' 2>/dev/null >families/$j.list;
    echo "---3" families/$j.list;
    k=`ls generateds/$i/annotation/*annotated.sses.json| sed "s|generateds/$i/annotation/||" | sed s"|-annotated.sses.json||"| head -n1`;
    #./transform_sses_to_gnuplot.py --directory generateds/$i/annotation/ --protein $k 
    echo "k je " $k
    ./transform_sses_to_gnuplot.py --directory generateds/$i/annotation/ --protein $k > families/$j.list.sses;
    echo "---4"  $i
    ./transform_channels_to_gnuplot.py --directory generateds/$i/annotation/ --unknown_file_list families/$j.list >families/$j.list.gnuplot;
    ./transform_channels_to_gnuplot.py --directory generateds/$i/annotation/ --unknown_file_list families/$j.list --near >families/$j.list.gnuplot_near;
    echo "---5"  $i
    mkdir -p gnuplots/gnuplots-$j
    mv gnuplot-* gnuplots/gnuplots-$j
    echo "families/$j.list.gnuplot 0.2 1,0,0" > families/ch_$j
    echo "families/$j.list.sses 0.5 0,0,0" >> families/ch_$j
    echo "families/$j.list.gnuplot_near 0.2 1,0,0" > families/ch_near_$j
    echo "families/$j.list.sses 0.5 0,0,0" >> families/ch_near_$j
done
