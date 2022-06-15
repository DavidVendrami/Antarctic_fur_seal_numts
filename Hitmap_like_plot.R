### Need to redo the carnival figure for each fragment separately.
### Run blast on each species 
blastn ...

# Now parse the outputs, remove mitochondrail hits, and retrieve information (numt_frag, species, query_cov, perc_identity)

# For one species:

# Fields: query acc.ver, subject acc.ver, % identity, alignment length, mismatches, gap opens, q. start, q. end, s. start, s. end, evalue, bit score

grep -v '#' M_leonina.fna10kbclust_blast_out > temp_files/Mleo_blast_grepped.txt
cat temp_files/Header.txt temp_files/Mleo_blast_grepped.txt > Mleo_blast_rinput.txt

R
fai<-read.table("../All_numts.fasta.fai",h=F)
data<-read.table("../R_input/",h=T)
spec<-"Z_californicus" # Put here species name


data[1:50,1:5]


# bu<-data
data<-data[data$evalue<0.01,]
data<-data[data$Hit!="NC_008416.1",] # !!! Here put the mtDNA !!!
data_list<-split(data, f = data$Query)

# Create a df with the best hit for each numt
df<-data[FALSE,]
for (i in 1:length(data_list)){
bla<-data_list[[i]][order(data_list[[i]]$evalue),]
if (dim(bla)[1]!=0){
df[i,]<-bla[1,]
}}

# Output importnat info for the species
numt<-c()
sp<-c()
qc<-c()
perci<-c()

for (i in 1:25){
if (length(df$Align_length[which(fai[i,1]==df$Query)])!=0){
qc[i]<-(df$Align_length[which(fai[i,1]==df$Query)]-df$Gap_opens[which(fai[i,1]==df$Query)])/fai[i,2]
}

qc[qc>1]<-1

numt[i]<-as.character(fai[i,1])

sp[i]<-spec

if (length(df$Align_length[which(fai[i,1]==df$Query)])!=0){
perci[i]<-df$Perc_id[which(fai[i,1]==df$Query)]
}}

out<-data.frame(Numt=numt,Species=sp,Query_Cov=qc,Perc_ID=perci)
write.table(out,paste(spec,"_Sep_plot.txt",sep=""),col.names=T,row.names=F,quote=F,sep="\t")


# for loop:
# create file list and loop through it to import the various files
# Create also a character in R with the IDs of the mtDNA (same length of file list, so you can easily loop thrugh it with the same 
# loop variable)

# Once you have all the files, cat them:
mv C_ursinus_Sep_plot.txt C_ursinus_Sep_plot_Head.txt
cp C_ursinus_plot_Head.txt Final_all_sp.txt

for i in *plot.txt; do tail -25 $i | cat >> Final_all_sp.txt; done

# PLOT

library(ggplot2)
library(gridExtra)

df<-read.table("Final_all_sp_Clust_numts.txt",h=T)
#int<-c(1.5,2.5,3.5,4.5,5.5,8.5,10.5,11.5,13.5,14.5,15.5,16.5,24.5)
levs<-c("Z_californicus","E_jubatus","C_ursinus","O_rosmarus","H_gryous","P_vitulina","L_weddelli","M_angustirostris","M_leonina","M_schauislandi")
levs<-levs[length(levs):1]
bla<-factor(df$Species,levels=levs)
df$Species2<-bla
mne<-factor(df$Numt,levels=c("Numt_1","Numt_2","Numt_3","Numt_4","Numt_5","Numt_6","Numt_7","Numt_8","Numt_9","Numt_10","Numt_11","Numt_12","Numt_13","Numt_14","Numt_15","Numt_16","Numt_17","Numt_18","Numt_19","Numt_20","Numt_21","Numt_22","Numt_23","Numt_24","Numt_25"))
df$Numt2<-mne

tiff("Carnival2.tiff",width=11.9, height=4.85, units= 'in', res=600, pointsize=1/600)
ggplot(data=df, 
aes(x=Numt2,y=Species2, size=Query_Cov, color=Perc_ID)) +
geom_point(shape=15) +
scale_colour_gradientn(colours=rainbow(3,rev=T)) +
#theme_classic() +
scale_x_discrete(expand=c(0,1)) +
scale_y_discrete(expand=c(0,2)) +
geom_vline(xintercept=c(int),linetype="dashed",size=0.5,color="grey50")
dev.off()


# Separate fragments (a):

df<-read.table("Final_all_sp_Sep_sp.txt",h=T)
int<-c(1.5,5.5,7.5,8.5,10.5,11.5,13.5,17.5,18.5,22.5,23.5)
levs<-c("Z_californicus","E_jubatus","C_ursinus","O_rosmarus","H_grypus","P_vitulina","L_weddelli","M_angustirostris","M_leonina","M_schauislandi")
levs<-levs[length(levs):1]
bla<-factor(df$Species,levels=levs)
df$Species2<-bla
mne<-factor(df$Numt,levels=c("Numt_1","Numt_2_fragment_1","Numt_2_fragment_2","Numt_2_fragment_3","Numt_2_fragment_4","Numt_3_fragment_1","Numt_3_fragment_2","Numt_4",
"Numt_5_fragment_1","Numt_5_fragment_2","Numt_6","Numt_7_fragment_1","Numt_7_fragment_2","Numt_8_fragment_1","Numt_8_fragment_2","Numt_8_fragment_3","Numt_8_fragment_4",
"Numt_9","Numt_10_fragment_1","Numt_10_fragment_2","Numt_10_fragment_3","Numt_10_fragment_4","Numt_11","Numt_12_fragment_1","Numt_12_fragment_2","Numt_12_fragment_3",
"Numt_12_fragment_4","Numt_13","Numt_14","Numt_15_fragment_1","Numt_15_fragment_2","Numt_15_fragment_3","Numt_16_fragment_1","Numt_16_fragment_2","Numt_17_fragment_1",
"Numt_17_fragment_2","Numt_17_fragment_3","Numt_18_fragment_1","Numt_18_fragment_2","Numt_18_fragment_3","Numt_19_fragment_1","Numt_19_fragment_2","Numt_19_fragment_3",
"Numt_20_fragment_1","Numt_20_fragment_2","Numt_20_fragment_3","Numt_21_fragment_1","Numt_21_fragment_2","Numt_21_fragment_3","Numt_22_fragment_1","Numt_22_fragment_2",
"Numt_22_fragment_3","Numt_23_fragment_1","Numt_23_fragment_2","Numt_24_fragment_1","Numt_24_fragment_2","Numt_25"))
df$Numt2<-mne

df2<-df[order(df$Numt2),]

df<-df2[1:250,]



tiff("Carnival_sepA.tiff",width=11.9, height=4.85, units= 'in', res=600, pointsize=1/600)
ggplot(data=df, 
aes(x=Numt2,y=Species2, size=Query_Cov, color=Perc_ID)) +
geom_point(shape=15) +
scale_colour_gradientn(colours=rainbow(3,rev=T)) +
#theme_classic() +
scale_x_discrete(expand=c(0,1)) +
scale_y_discrete(expand=c(0,2)) +
geom_vline(xintercept=c(int),linetype="dashed",size=0.5,color="grey50")
dev.off()


# Separate fragments (b):

df<-read.table("Final_all_sp_Sep_sp.txt",h=T)
int<-c(2.5,3.5,4.5,7.5,9.5,12.5,15.5,18.5,21.5,24.5)
levs<-c("Z_californicus","E_jubatus","C_ursinus","O_rosmarus","H_grypus","P_vitulina","L_weddelli","M_angustirostris","M_leonina","M_schauislandi")
levs<-levs[length(levs):1]
bla<-factor(df$Species,levels=levs)
df$Species2<-bla
mne<-factor(df$Numt,levels=c("Numt_1","Numt_2_fragment_1","Numt_2_fragment_2","Numt_2_fragment_3","Numt_2_fragment_4","Numt_3_fragment_1","Numt_3_fragment_2","Numt_4",
"Numt_5_fragment_1","Numt_5_fragment_2","Numt_6","Numt_7_fragment_1","Numt_7_fragment_2","Numt_8_fragment_1","Numt_8_fragment_2","Numt_8_fragment_3","Numt_8_fragment_4",
"Numt_9","Numt_10_fragment_1","Numt_10_fragment_2","Numt_10_fragment_3","Numt_10_fragment_4","Numt_11","Numt_12_fragment_1","Numt_12_fragment_2","Numt_12_fragment_3",
"Numt_12_fragment_4","Numt_13","Numt_14","Numt_15_fragment_1","Numt_15_fragment_2","Numt_15_fragment_3","Numt_16_fragment_1","Numt_16_fragment_2","Numt_17_fragment_1",
"Numt_17_fragment_2","Numt_17_fragment_3","Numt_18_fragment_1","Numt_18_fragment_2","Numt_18_fragment_3","Numt_19_fragment_1","Numt_19_fragment_2","Numt_19_fragment_3",
"Numt_20_fragment_1","Numt_20_fragment_2","Numt_20_fragment_3","Numt_21_fragment_1","Numt_21_fragment_2","Numt_21_fragment_3","Numt_22_fragment_1","Numt_22_fragment_2",
"Numt_22_fragment_3","Numt_23_fragment_1","Numt_23_fragment_2","Numt_24_fragment_1","Numt_24_fragment_2","Numt_25"))
df$Numt2<-mne

df2<-df[order(df$Numt2),]

df<-df2[251:500,]



tiff("Carnival_sepB.tiff",width=11.9, height=4.85, units= 'in', res=600, pointsize=1/600)
ggplot(data=df, 
aes(x=Numt2,y=Species2, size=Query_Cov, color=Perc_ID)) +
geom_point(shape=15) +
scale_colour_gradientn(colours=rainbow(3,rev=T)) +
#theme_classic() +
scale_x_discrete(expand=c(0,1)) +
scale_y_discrete(expand=c(0,2)) +
geom_vline(xintercept=c(int),linetype="dashed",size=0.5,color="grey50")
dev.off()



# Separate fragments (c):

df<-read.table("Final_all_sp_Sep_sp.txt",h=T)
int<-c(2.5,4.5,6.5)
levs<-c("Z_californicus","E_jubatus","C_ursinus","O_rosmarus","H_grypus","P_vitulina","L_weddelli","M_angustirostris","M_leonina","M_schauislandi")
levs<-levs[length(levs):1]
bla<-factor(df$Species,levels=levs)
df$Species2<-bla
mne<-factor(df$Numt,levels=c("Numt_1","Numt_2_fragment_1","Numt_2_fragment_2","Numt_2_fragment_3","Numt_2_fragment_4","Numt_3_fragment_1","Numt_3_fragment_2","Numt_4",
"Numt_5_fragment_1","Numt_5_fragment_2","Numt_6","Numt_7_fragment_1","Numt_7_fragment_2","Numt_8_fragment_1","Numt_8_fragment_2","Numt_8_fragment_3","Numt_8_fragment_4",
"Numt_9","Numt_10_fragment_1","Numt_10_fragment_2","Numt_10_fragment_3","Numt_10_fragment_4","Numt_11","Numt_12_fragment_1","Numt_12_fragment_2","Numt_12_fragment_3",
"Numt_12_fragment_4","Numt_13","Numt_14","Numt_15_fragment_1","Numt_15_fragment_2","Numt_15_fragment_3","Numt_16_fragment_1","Numt_16_fragment_2","Numt_17_fragment_1",
"Numt_17_fragment_2","Numt_17_fragment_3","Numt_18_fragment_1","Numt_18_fragment_2","Numt_18_fragment_3","Numt_19_fragment_1","Numt_19_fragment_2","Numt_19_fragment_3",
"Numt_20_fragment_1","Numt_20_fragment_2","Numt_20_fragment_3","Numt_21_fragment_1","Numt_21_fragment_2","Numt_21_fragment_3","Numt_22_fragment_1","Numt_22_fragment_2",
"Numt_22_fragment_3","Numt_23_fragment_1","Numt_23_fragment_2","Numt_24_fragment_1","Numt_24_fragment_2","Numt_25"))
df$Numt2<-mne

df2<-df[order(df$Numt2),]

df<-df2[501:570,]



tiff("Carnival_sepC.tiff",width=5, height=4.85, units= 'in', res=600, pointsize=1/600)
ggplot(data=df, 
aes(x=Numt2,y=Species2, size=Query_Cov, color=Perc_ID)) +
geom_point(shape=15) +
scale_colour_gradientn(colours=rainbow(3,rev=T)) +
#theme_classic() +
scale_x_discrete(expand=c(0,1)) +
scale_y_discrete(expand=c(0,2)) +
geom_vline(xintercept=c(int),linetype="dashed",size=0.5,color="grey50")
dev.off()
