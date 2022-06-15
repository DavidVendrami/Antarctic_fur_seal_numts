### Step 1:
### get AFS mt and nuclear genomes and find mt sequences within the genome.
## We will look for the presence of mt genome bits in the nuclear genome with a blastn search as in Lammers et al. 2017.

## First, make a database (from nuclear genome):
makeblastdb -in antarctic_fur_seal_02Jun2018_WAj4l.fasta -parse_seqids -dbtype nucl

## Now do the blast search using a word size of 20:
blastn -db AFS_nucl_newHeaders.fasta -query ../AFS_mitoGenome.fasta -outfmt 7 -word_size 20 -out numt_nucdb

## Now remove '#' from the output of the search:
grep -v '#' numt_nucdb > numt_nucdb.txt

# Let's switch to R and fill in the table "numt_report.txt"
# you basically read in "numt_mtdb.txt", keep only matches equal or longer than 200 (equaling >95% of the cumulative numt length), and fill the table.
# matches within the same nuclear contig separated by no more than 10kb are merged (length do not include gap, gap info under "distance.. column")
# for practicity split the dataframe into a list by Query ID:
numt<-split(data,f=data$Query)

## We will then need to extract sequences to check for insert (make a .bed file to take out only the inserts) and whether primers have been designed (extract 100 bp more from both sides).
## use bedtools getfasta (https://bedtools.readthedocs.io/en/latest/content/tools/getfasta.html)

## Make bed file
## bed file = 3 columns. col 1=chr_ID, col 2=start position (will be excluded), col 3=end position (included)
~/bin/bedtools2/bin/bedtools getfasta -fi ../antarctic_fur_seal_02Jun2018_WAj4l.fasta -bed get_numt.txt -fo ./numts.fasta

## Create primers for every numt (so that to test with PCR if they are real rather than artefacts). One set of primer per numt (select one fragment in case there are multiple). One primer (F) half in and half out of
## the numt and the other primer (R) completely in. PCR product up to 500bp. Use primer3 (web interface). 


# numts phylogenetic trees. Get mtDNA of the other pinniped species (+ the dog) and build phylogenies separately for each relevant numt.


## Now let's try to blast each AFS numt to each mt ref genome
## First you need to self-concatenate the ref mt genomes (so that you dont miss alignments in case they go over the end - mt genome is circular)
for i in *.fasta; do cat $i $i > ../${i%.fasta}_selfConc.fasta; done
## Now open them in nano and concatenate them
## Then make blast database:
for i in *.fasta; do makeblastdb -in $i -parse_seqids -dbtype nucl; done &
## Now extract the first numt (where numts.fasta is a fasta file containing the AFS numt squences):
sed -n '1,2p' numts/numts.fasta > numt_1.fasta
## Let's try to blast  it to some species:
for i in *.fasta; do blastn -db $i -query ./numts/numt_25.fasta -outfmt 7 -word_size 20 -out ${i%.fasta}_numt25_blast; done
## Let's create a bunch of bed files to extract the blasted sequences (*_numt1_blast.txt)
for i in *_blast; do sed -n '6p' $i | cut -d$'\t' -f2,10 > $i.txt; done
for i in *_blast; do sed -n '6p' $i | cut -d$'\t' -f9 > $i.sec; done
for i in *.txt; do paste -d$'\t' $i ${i%txt}sec > ${i%txt}bed; done
## for bedtools to work you need fasta files with sequences stored in line of same length,
## let's then put the seuqnece in a single line (using Kanchon's script "convert_to_one_line_fasta.pl")
for i in *.fasta; do ./convert_to_one_line_fasta.pl $i; done
# However for some unexplainable reasons most one-line fasta files will have ^M as new lines. Correct by opening these files in vi and entering
# :%s/^M//g # type "^M" as Ctrl + vm
# :wq # to save and exit
# or better (put it in for loop):
sed 's/^M//g' Miroun_leo_selfConc.fasta_one_line.fas > Miroun_leo_selfConc_OL.fasta # type "^M" as Ctrl + vm
for i in *fasta_one_line*; do sed 's/^M//g' $i > ${i%.fasta_one_line.fas}_OL.fasta; done
# Now extract the mt sequences correponding to numt_1 from all species
for i in *bed; do ~/bin/bedtools2/bin/bedtools getfasta -fi ../OL_fastas/${i%_numt1_blast.bed}_OL.fasta -bed $i -fo ./${i%_numt1_blast.bed}_numt_1.fasta; done
# and concatenate them in a single fasta file:
cat *.fasta > numt_1_all.fasta


# In case of numts with multiple fragments:
# blast 
for i in *.fasta; do blastn -db $i -query ./numt_2/numt_2.fasta -outfmt 7 -word_size 20 -out ${i%.fasta}_numt2_blast; done
# check they have all fragments:
cat *_numt2_blast | less
# make bed
for i in *_blast; do sed -n '6p' $i | cut -d$'\t' -f2,9 > $i.txt; done
for i in *_blast; do sed -n '6p' $i | cut -d$'\t' -f10 > $i.sec; done
for i in *.txt; do paste -d$'\t' $i ${i%txt}sec > ${i%.txt}_1.bed; done

for i in *_blast; do sed -n '8p' $i | cut -d$'\t' -f2,9 > $i.txt; done
for i in *_blast; do sed -n '8p' $i | cut -d$'\t' -f10 > $i.sec; done
for i in *.txt; do paste -d$'\t' $i ${i%txt}sec > ${i%.txt}_2.bed; done

for i in *_blast; do sed -n '10p' $i | cut -d$'\t' -f2,10 > $i.txt; done
for i in *_blast; do sed -n '10p' $i | cut -d$'\t' -f9 > $i.sec; done
for i in *.txt; do paste -d$'\t' $i ${i%txt}sec > ${i%.txt}_2.bed; done

for i in *_blast; do sed -n '12p' $i | cut -d$'\t' -f2,10 > $i.txt; done
for i in *_blast; do sed -n '12p' $i | cut -d$'\t' -f9 > $i.sec; done
for i in *.txt; do paste -d$'\t' $i ${i%txt}sec > ${i%.txt}_2.bed; done

# bedtools
for i in *_1.bed; do ~/bin/bedtools2/bin/bedtools getfasta -fi ../OL_fastas/${i%_numt2_blast_1.bed}_OL.fasta -bed $i -fo ./${i%_numt2_blast_1.bed}_numt_2.fasta_1; done
for i in *_2.bed; do ~/bin/bedtools2/bin/bedtools getfasta -fi ../OL_fastas/${i%_numt2_blast_2.bed}_OL.fasta -bed $i -fo ./${i%_numt2_blast_2.bed}_numt_2.fas_2; done
for i in *_3.bed; do ~/bin/bedtools2/bin/bedtools getfasta -fi ../OL_fastas/${i%_numt2_blast_3.bed}_OL.fasta -bed $i -fo ./${i%_numt2_blast_3.bed}_numt_2.fas_3; done
for i in *_4.bed; do ~/bin/bedtools2/bin/bedtools getfasta -fi ../OL_fastas/${i%_numt2_blast_4.bed}_OL.fasta -bed $i -fo ./${i%_numt1_blast_4.bed}_numt_2.fas_4; done
# manipulate a bit and put the entire numt in 1 fasta file with 1 header
for i in *.fas_2; do grep -v '>' $i > ${i%fas_2}fasta_2; done
for i in *.fas_3; do grep -v '>' $i > ${i%fas_3}fasta_3; done
for i in *.fas_4; do grep -v '>' $i > ${i%fas_4}fasta_4; done

for i in *.fasta_1; do cat $i ${i%_1}_2 ${i%_1}_3 ${i%_1}_4 > ${i%fasta_1}almost; done

mv numt_2.fasta numt_2_complete.almost

for i in *almost; do ../convert_to_one_line_fasta.pl $i; done
for i in *almost_one_line.fas; do mv "$i" "$(echo $i | sed -e s/.almost_one_line.fas/complete.fasta/)"; done

cat *complete.fasta > numt_2_all.almost_fasta
sed 's/>/\n>/g' numt_2_all.almost_fasta | sed '1d' > numt_2_all.fasta

### In case of revcomp:
# If its only one fragment, you can revcomp in mega
# Otherwise:
for i in *selfConc*fasta_1; do ../revcomp.txt $i; done
for i in *selfConc*fas_2; do ../revcomp.txt $i; done
for i in *selfConc*fas_3; do ../revcomp.txt $i; done
for i in *selfConc*fas_4; do ../revcomp.txt $i; done

for i in *.fas_2revcomp.fas; do grep -v '>' $i > ${i%fas_2revcomp.fas}fasta_2revcomp.fasta; done
for i in *.fas_3revcomp.fas; do grep -v '>' $i > ${i%fas_3revcomp.fas}fasta_3revcomp.fasta; done
for i in *.fas_4revcomp.fas; do grep -v '>' $i > ${i%fas_4revcomp.fas}fasta_4revcomp.fasta; done
for i in *.fasta_1revcomp.fas; do mv "$i" "$(echo $i | sed -e s/fasta_1revcomp.fas/fasta_1revcomp.fasta/)"; done

for i in *.fasta_1revcomp.fasta; do cat $i ${i%_1revcomp.fasta}_2revcomp.fasta ${i%_1revcomp.fasta}_3revcomp.fasta ${i%_1revcomp.fasta}_4revcomp.fasta > ${i%fasta_1revcomp.fasta}revcomp.fastaalmost; done

mv numt_22.fasta numt_22_completerevcomp.fastaalmost

for i in *revcomp.fastaalmost; do ../convert_to_one_line_fasta.pl $i; done
for i in *revcomp.fastaalmost_one_line.fas; do mv "$i" "$(echo $i | sed -e s/.fastaalmost_one_line.fas/complete.fasta/)"; done

cat *revcompcomplete.fasta > numt_22_all.revcompalmost_fasta
sed 's/>/\n>/g' numt_22_all.revcompalmost_fasta | sed '1d' > numt_22_revcomp_all.fasta

# Use MEGA to align and trim sequences
# For phylogenetic trees use MEGA 
