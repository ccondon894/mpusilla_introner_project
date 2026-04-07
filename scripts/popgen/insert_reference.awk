#! /usr/bin/awk -f 

BEGIN {
    FS = "\t"
    OFS = FS
}

#print the header info as-is
/^##/ {
    print
    next
}

#add sample named "CCMP1545" to the list of samples
/^#CHROM/ {
    print $0"\tCCMP1545"
    next
}

#add homozygous reference allele to every locus.
{
    print $0"\t0:100,0:100:100:0,1000"
}
