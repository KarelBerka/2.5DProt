#!/bin/bash
ddir=$1
update=$2


echo "tune multiple opacity layout, s35 use family average structure - regenerated after 2021-04-26"
#grep "all" $update | sed "s|all .*||"

for i in `grep "all" $update | sed "s|all .*||"`; do
#    echo $i, `ls $ddir `
    j=`ls $ddir | grep generated-$i- | head -n1| sed "s|generated-$i-||" | sed "s|T.*||" | tr -d "-"`
    if [ "${j}0" -le "202104250" ]; then
#       k=`grep "$i " $update | grep "all " | sed "s|.* all ||" | sed "s| .*||"`
#       if [ $k -ge 3000 ]; then
           grep "$i " $update | grep "all " 
#       fi
    fi
    # TODO: preclusters test should be fixed
done  # | wc -l
#exit


#echo "preclusters missing - regenerated after 2021-04-05"
#for i in `grep "all" $update | sed "s|all .*||"`; do
#
#    j=`ls $ddir | grep generated-$i- | head -n1| sed "s|generated-$i-||" | sed "s|T.*||" | tr -d "-"`
#    if [ "${j}0" -le "202104050" ]; then 
#       k=`grep "$i " $update | grep "all " | sed "s|.* all ||" | sed "s| .*||"`
#
#       if [ $k -ge 3 ]; then
#           grep "$i " $update | grep "all " 
#       fi
#    fi
#    # TODO: still present preclusters test should be fixed
#done # | wc -l

echo "STATISTIC SUMMARIES =========================="
p=`./cath-family-updates.py --show-last-updated-per-family | wc -l`
echo $p, "Sum of all families"


p=`find $ddir -name "imageCircle-multiple.png" | wc -l`
echo $p, "Sum of all multiple layouts"


l=`./cath-family-updates.py --show-last-updated-per-family | sed "s| .*||"`

p=`for i in  $l; do
   ls  $ddir/generated-$i-* >/dev/null 2>/dev/null
   if [ $? -ne 0 ]; then 
      echo $i
   fi
done | wc -l`
#for i in `find $ddir -name "generated-*`"
echo $p, "Empty families"

p=`for i in  $l; do
   ls  $ddir/generated-$i-*/imageCircle-multiple.png >/dev/null 2>/dev/null
   if [ $? -ne 0 ]; then 
      ls  $ddir/generated-$i-* >/dev/null 2>/dev/null
      if [ $? -eq 0 ]; then
          ls  $ddir/generated-$i-*/*fail* >/dev/null 2>/dev/null
          if [ $? -eq 0 ]; then 
#            echo $i
             grep all $update |grep $i
          fi
      fi
   fi
done `
echo `echo $p | wc -w`, "some fails listed below"#, $p


p=`for i in  $l; do
   ls  $ddir/generated-$i-*/imageCircle-multiple.png >/dev/null 2>/dev/null
   if [ $? -ne 0 ]; then 
      ls  $ddir/generated-$i-* >/dev/null 2>/dev/null
      if [ $? -eq 0 ]; then
          ls  $ddir/generated-$i-*/*fail* >/dev/null 2>/dev/null
          if [ $? -ne 0 ]; then 
             #echo $i
             grep all $update |grep $i
          fi
      fi
   fi
done `
echo `echo $p | wc -w`, "crucial fail ", $p


p=`find $ddir -name "layout*" | wc -l`
echo $p, "Sum of all single layouts"

echo "ERRORS: ======================================"

p=`find $ddir -name "*fail*" |wc -l`
q=`find $ddir -name "*fail*" | sed "s|.*generated-||" | sed "s|-.*||"| sort | uniq ` # | wc -l` 
echo "all errors ...................", $p, $q

p=`find $ddir -name "get_domain*fail*" | wc -l `
find $ddir -name "get_domain*fail*" | sed "s|[^-]*-||" | sed "s|-.*||" | sort | uniq >list_gdb
q=`find $ddir -name "get_domain*fail*" | sed "s|[^-]*-||" | sed "s|-.*||" | sort | uniq | wc -l`
echo "get_domain_boundaries fail ...", $p, $q

p=`find $ddir -name "range*fail*" | wc -l`
echo "domain_range fails ...........", $p

pom1=`find $ddir -name "annotated_and_detected*fail*"`
echo "annotated_and_detected celkem ", `echo $pom1 | wc -w`

pom=`find $ddir -name "*SecStr*fail*"`
l=`echo $pom | wc -w`
echo "SecStr celkem ................", $l

p=`find $ddir -name "generate*layout*" | sed "s|[^-]*-||" | sed "s|-.*||" | wc -l`
echo "generate_layout fails ........", $p

p=`find $ddir -name "*2dlayout*fail*" | wc -l`
echo "show-2dlayout fails ..........", $p

p=`find $ddir -name "ubertemplate-fail*" | wc -l`
echo "ubertemplate fails ..........", $p

p=`find $ddir -name "*fail*" ! -name "get_domain*fail*" ! -name "range*fail*" ! -name "annotated_and*fail*" ! -name "*SecStr*fail*" ! -name "generate*layout*" ! -name "*2dlayout*fail*" ! -name "ubertemplate-fail*"`
echo "The rest .....................", $p

echo "SECSTR details ==============================="
echo -e "stderr: WARNING: Loaded structure is of type 'unknown' (-R) ..........................: \c"
for i in `echo $pom`; do 
    grep  "stderr: WARNING: Loaded structure is of type 'unknown'" $i > /dev/null 2>/dev/null; 
    if [ $? -eq 0 ]; then 
           p=`echo $i | sed "s|[^-]*-||" | sed "s|-.*||" `;
#           echo $p
           grep "$p " update_2021-02-06* | grep "all " 

#        echo $i
    fi;
done | wc -l


echo -e "ERROR: Query protein does not contain chain (-R) .....................................: \c"
for i in `echo $pom`; do 
    grep  "ERROR: Query protein does not contain chain" $i > /dev/null 2>/dev/null; 
    if [ $? -eq 0 ]; then 
#       echo $i
        p=`echo $i | sed "s|[^-]*-||" | sed "s|-.*||" `;
#       echo $p
        grep "$p " update_2021-02-06* | grep "all " 

    fi;
done | wc -l


echo -e "stderr: WARNING: Found some atoms with identical chain ID, residue number (-R) .......: \c"
for i in `echo $pom`; do
    grep   "stderr: WARNING: Found some atoms with identical chain ID, residue number" $i > /dev/null 2>/dev/null;
    if [ $? -eq 0 ]; then
             p=`echo $i | sed "s|[^-]*-||" | sed "s|-.*||" `;
#             echo $p
             grep "$p " update_2021-01-31* | grep "all " 

#       echo $i
    fi;
done | wc -l


echo -e "ERROR: Template protein does not contain chain A (R+) ................................: \c"
for i in `echo $pom`; do 
    grep  "ERROR: Template protein does not contain chain A." $i > /dev/null 2>/dev/null; 
    if [ $? -eq 0 ]; then 
        echo $i
##        p=`echo $i | sed "s|[^-]*-||" | sed "s|-.*||" `;
##        echo $p
##           grep "$p " update_2021-01-10T22* | grep "all " | sed "s|all.*||"

    fi;
done | wc -l

echo -e "Unhandled Exception: System.FormatException: Cannot fit CIF data to PDB data model(R+): \c"
for i in `echo $pom`; do 
    grep  "Cannot fit CIF data to PDB data model: CIF file" $i > /dev/null 2>/dev/null; 
    if [ $? -eq 0 ]; then 
             p=`echo $i | sed "s|[^-]*-||" | sed "s|-.*||" `;
             echo $p
#             grep "$p " update_2021-01-31* | grep "all " 

#        echo $i
    fi;
done | wc -l

echo -e "Unhandled Exception: Cif.CifException: Lexical error on line ........................,: \c"
for i in `echo $pom`; do 
    grep  "Cif.CifException: Lexical error" $i > /dev/null 2>/dev/null; 
    if [ $? -eq 0 ]; then 
        echo $i
    fi;
done | wc -l

echo -e "Unhandled Exception: Cif.CifException: Syntactic error on line (R+) ..................: \c"
for i in `echo $pom`; do 
    grep  "Cif.CifException: Syntactic error" $i > /dev/null 2>/dev/null; 
    if [ $? -eq 0 ]; then 
        echo $i
    fi;
done | wc -l

echo -e "Index was outside the bounds of the array.(R -) ......................................: \c"
for i in `echo $pom`; do 
    grep  "Index was outside the bounds of the array." $i > /dev/null 2>/dev/null; 
    if [ $? -eq 0 ]; then 
        echo $i
    fi;
done | wc -l

echo -e "The rest ............................................................................ : \c"
for i in `echo $pom`; do 
    grep  "stderr: WARNING: Loaded structure is of type 'unknown'" $i >/dev/null 2>/dev/null; 
    if [ $? -ne 0 ]; then 
         grep  "ERROR: Query protein does not contain chain" $i >/dev/null 2>/dev/null; 
         if [ $? -ne 0 ]; then 
               grep  "stderr: WARNING: Found some atoms with identical chain ID, residue number" $i > /dev/null 2>/dev/null; 
               if [ $? -ne 0 ]; then 
                     grep  "ERROR: Template protein does not contain chain" $i > /dev/null 2>/dev/null; 
                     if [ $? -ne 0 ]; then 
                          grep  "Cannot fit CIF data to PDB data model: CIF file" $i > /dev/null 2>/dev/null; 
                          if [ $? -ne 0 ]; then
                              grep  "Cif.CifException: Lexical error" $i > /dev/null 2>/dev/null; 
                              if [ $? -ne 0 ]; then
                                   grep  "Cif.CifException: Syntactic error" $i > /dev/null 2>/dev/null; 
                                   if [ $? -ne 0 ]; then
                                       grep  "Index was outside the bounds of the array." $i > /dev/null 2>/dev/null; 
                                       if [ $? -eq 0 ]; then 
                                          echo $i
                                       fi;
                                   fi;
                              fi;
                          fi;
                     fi;
               fi;
         fi;
    fi;
done | wc  -l

pom1=`find $ddir -name "annotated_and_detected*fail*"`
echo "annotated_and_detected details ==============="


echo -e "error: argument annotated: can't open 'generated- (-R) ...............................: \c"
for i in `echo $pom1`; do 
    grep  "error: argument annotated: can't open 'generated-" $i > /dev/null 2>/dev/null; 
    if [ $? -eq 0 ]; then 
#             p=`echo $i | sed "s|[^-]*-||" | sed "s|-.*||" `;
#            echo $p
#             grep "$p " update_2021-01-10T22* | grep "all "
#             p=`echo $i | sed "s|[^-]*-||" | sed "s|-.*||" `;
#             grep "$p " update_2021-01-31* | grep "all " 
        echo $i
    fi;
done | wc  -l

echo -e "error: argument detected: can't open 'generated- (R+) ................................: \c"
for i in `echo $pom1`; do 
    grep  "error: argument detected: can't open 'generated-" $i > /dev/null 2>/dev/null; 
    if [ $? -eq 0 ]; then 
#             p=`echo $i | sed "s|[^-]*-||" | sed "s|-.*||" `;
#            echo $p
#            grep "$p " update_2021-01-31T22* | grep "all "
        echo $i
    fi;
done | wc -l

echo -e "The rest .............................................................................: \c"
for i in `echo $pom1`; do 
    grep  "error: argument annotated: can't open 'generated-" $i >/dev/null 2>/dev/null; 
    if [ $? -ne 0 ]; then 
          grep  "error: argument detected: can't open 'generated-" $i >/dev/null 2>/dev/null; 
          if [ $? -ne 0 ]; then 
               echo $i
          fi;
    fi;
done  | wc -l

echo "OTHER PROBLEMS ============================"
echo "Missing color tag (HAVE TO BE REGENERATED) "

for i in `find $ddir -name consensus-template.sses.json |grep -v "consensus_A00" `; do 
    grep color $i >/dev/null 2>/dev/null; 
    p=$?;
    if [ $p -ne 0 ]; then
         grep label $i >/dev/null 2>/dev/null; 
         p=$?;
         if [ $p -eq 0 ]; then
             p=`echo $i | sed "s|[^-]*-||" | sed "s|-.*||" `;
#             echo $p
             grep "$p " update_2021-01-31* | grep "all " 
         fi
    fi;
done 
#| wc -l

echo "Missing sync tag (HAVE TO BE REGENERATED) "
for i in `grep all $update | sed "s| .*||"`; do 
    j=`ls $ddir/generated-$i-*/main-*.log 2>/dev/null`
    if [ $? -eq 0 ]; then 
       grep "was synced" $j >/dev/null 2>/dev/null;
       p=$?;
       if [ $p -ne 0 ]; then   
           p=`echo $i | sed "s|[^-]*-||" | sed "s|-.*||" `;
            grep "$p " $update | grep "all " 
        fi;
    fi;
    # TODO sync probles should end with fail log - will be deal as a standart error, not a special one
done  | wc -l


#PRODUCE L1,2,3,SSES

#DIFF L1,2,3 SSES if there is any old timestamp given

#echo $pom1
