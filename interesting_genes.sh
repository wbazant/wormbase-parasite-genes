#!/usr/bin/bash
# Look at iRAP quantification results in ISL
# - per exon
# - per gene
# - look for cases where they're far off
# Try to find "interesting" gene models


# In:
# Smp_154160:001     9.31
# Smp_154160:002     4.37
# Smp_154160:003     0.18
# Smp_154160:004     0
summarise_exons() {
perl -e '
my ($current_gene,$current_gene_exons, $first_exon_tpm, $last_exon_tpm); 
sub flush {
    print join "\t", $current_gene, 0+$current_gene_exons, $first_exon_tpm, $last_exon_tpm . "\n"
       if $current_gene;
}
<>;
while(<>){
   if(s/^$current_gene://){
      ($current_gene_exons, $last_exon_tpm) = split;
   } else {
     flush();
     /([^:]+):\d+\s+(.*)/;
     $current_gene = $1;
     $current_gene_exons = "001";
     $first_exon_tpm = $2;
     $last_exon_tpm = $2;
   }
}
flush();
' "$@"
}

# Skip lines of the sort: 
# Smp_019920+Smp_133030:005       12.9
summarise_exons_clean(){
  perl -ne 'print unless /^\w+\+\w+/' "$@" \
   | summarise_exons \
   |  sort -k 1b,1
}

join_genes(){
  join -t $'\t' - "$@"
}

# Sauce formula
# Main: much higher expression in whole gene than in the boundary exons
# More exons is better (this is debatable), and has to be more than one exon
# Expression above 0.5 TPM
# Log and round for prettier numbers
score(){
  perl -MList::Util -ne '
    chomp;
    my ($gene_id, $num_exons, $e1, $en, $g) = split;
    printf("$_\t%.1f\n", log (List::Util::max(0.3678, ($num_exons - 1) * ($g -0.5) * (1/($e1+0.01) + 1/($en+0.01)))));
'
}

score_run(){
  path=$1
  summarise_exons_clean $path/$(basename $path).*.exons.tpm.dexseq.irap.tsv \
    | join_genes <( sort -k 1b,1 $path/$(basename $path).*.genes.tpm.kallisto.irap.tsv ) \
    | score
}

evidence_from_run(){
   score_run "$@" | cut -f 1,6
}

join_results(){
   perl -e 'my $cmd = shift; print "$cmd ", shift, join (" ",  map { " | join - <( $cmd $_ )"} @ARGV), "\n"; ' "$@" | source /dev/stdin
}

aggregate_max(){
   perl -MList::Util -lane 'printf("%s\t%.1f\n", shift @F, List::Util::max(@F));'
}
aggregate(){
   perl -MList::Util -lane 'printf("%s\t%.1f\n", shift @F, List::Util::sum(map {$_**2} @F));'
}
# Not useful for finding bad genes, mostly picks up differentially expressed genes
aggregate_variance(){
   perl -MList::Util -lane 'printf("%s\t%.1f\n", shift @F, List::Util::sum(map {$_**2} @F) - (List::Util::sum(@F)**2));'
}
# Works if you're on the EBI network only :) 
ftp_to_ebi_folder(){
   perl -e 'print join " ", map {s/ftp:\/\/ftp.ebi.ac.uk/\/ebi\/ftp/; $_} @ARGV' "$@"
}
genes_by_evidence(){
   join_results evidence_from_run $( ftp_to_ebi_folder "$@" ) | aggregate_max | sort -nr -k2 
}

populate_results(){
   species=$1;shift
   DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null && pwd )"
   [ -s "$DIR/out/${species}.interesting_genes.tsv" ] || genes_by_evidence "$@" > "$DIR/out/${species}.interesting_genes.tsv"
}  

get_isl_metadata_nonrepetitive_order(){
   perl -e 'my %seen ; <> ; while(<>){chomp; my @F = split "\t" ; my @C = @F; print join( "\t",$seen{join "", @C[0], splice (@C, 6) }++ , @F)."\n";}' "$@"\
     | perl -e 'my %seen ; <> ; while(<>){chomp; my @F = split "\t" ; my @C = @F; print join( "\t",$seen{@F[2]}++ , @F)."\n";}' \
     | sort -n | cut -f 3-
}

run_all(){
   d=$(ftp_to_ebi_folder ftp://ftp.ebi.ac.uk/pub/databases/wormbase/parasite/releases/.next/rnaseq)
   for f in $d/* ; do
      species=$(basename $f | sed 's/.tsv//' )
      num_runs=$(wc -l < $f )
      if [ $num_runs -gt 50 ]; then
         files=$(get_isl_metadata_nonrepetitive_order $f | grep -v "under 50%" | grep -v "200k" | cut -f 6 | head -n 15 )
      elif [ $num_runs -gt 15 ]; then
         files=$(get_isl_metadata_nonrepetitive_order $f | cut -f 6 | head -n 15 )
      else
         files=$(tail -n+2 $f | cut -f 6 )
      fi
      echo populate_results $species $files
      populate_results $species $files
   done
}
if [ "$0" = "$BASH_SOURCE" ] ; then
   cd $( dirname $0)
   rm out/*
   mkdir -p out log
   logPath=log/run_all.WBPS${PARASITE_VERSION}.log
   echo populate_results schisto_test | tee /dev/stderr >> $logPath
   populate_results schisto_test \
     ftp://ftp.ebi.ac.uk/pub/databases/arrayexpress/data/atlas/rnaseq/SRR505/002/SRR5054482 \
     ftp://ftp.ebi.ac.uk/pub/databases/arrayexpress/data/atlas/rnaseq/ERR582/ERR582667 \
     ftp://ftp.ebi.ac.uk/pub/databases/arrayexpress/data/atlas/rnaseq/SRR321/008/SRR3211868 \
     ftp://ftp.ebi.ac.uk/pub/databases/arrayexpress/data/atlas/rnaseq/ERR022/ERR022883 \
     ftp://ftp.ebi.ac.uk/pub/databases/arrayexpress/data/atlas/rnaseq/SRR629/SRR629229 \
     ftp://ftp.ebi.ac.uk/pub/databases/arrayexpress/data/atlas/rnaseq/SRR320/007/SRR3209257 \
     ftp://ftp.ebi.ac.uk/pub/databases/arrayexpress/data/atlas/rnaseq/ERR506/ERR506079 \
     | tee /dev/stderr >> $logPath

   run_all | tee /dev/stderr >> $logPath
fi
