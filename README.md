Usage
=====

Show the matrix of distances of helices:

    ./show input/stejne/1og5 10

Show the matrix of distances of helices:

 * at first it is the distance of average points
 * distances of the line segment from the first atom to the last athom
 * distances of the aproximation line segment

    ./matrix.py input/stejne/log5

New usage
---------

    ./cleanup-full.sh; export bilkoviny_template_selection_method=ubertemplate; export bilkoviny_ubertemplate_count=100; export bilkoviny_domains_count=20; for f in 1.10.630.10; do ./generate-new.sh cathid:$f; done

Links
=====

Posilam informace k anotacim tech odlisnych cytochromu.
Zde je tabulka s anotacemi:
https://drive.google.com/drive/folders/0B8mglzpeWWb-ZnNUeDRwcnZ0VXM
A zde jsou pak i dalsi informace:
https://drive.google.com/drive/folders/0B8mglzpeWWb-ZXNOR2U4ZHQ2bXc?usp=sharing

Requirements
============

gnuplot-py
----------
Best to do this in virtalenv:

    wget http://downloads.sourceforge.net/project/gnuplot-py/Gnuplot-py/1.8/gnuplot-py-1.8.tar.gz
    tar xzf gnuplot-py-1.8.tar.gz
    cd gnuplot-py-1.8/
    python setup.py install

Packages
--------
In Fedora, for rest of dependencies just run:

    dnf -y install python2-biopython python2-numpy python2-tabulate python2-scikit-learn
