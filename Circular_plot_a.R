data<-read.table("mtDNA_map.txt",h=T) #
nst<-data$start-1
nsto<-data$stop
stt<-(360*nst)/16156
stot<-(360*nsto)/16156
tdat<-data.frame(Circ_start=stt,Circ_stop=stot,Color=as.character(data$col))

numts<-read.table("numts.txt",h=T)
bla<-gsub("&","#",numts$Color)
numts$Color<-bla
nust<-numts$Start-1
nusto<-numts$Stop
nstt<-(360*nust)/16156
nstot<-(360*nusto)/16156
num_t<-data.frame(Circ_start=nstt,Circ_stop=nstot,Color=as.character(numts$Color))

tiff("Figure_1_Circ.tiff",width=7, height=7, units= 'in', res=600, pointsize=1/600)
par(mar = c(2, 2, 2, 2))
plot(c(-1.75, 1.75), c(-1.85, 1.85), type = "n", axes = FALSE, ann = FALSE, asp = 1)
for (i in 1:dim(tdat)[1]){
draw.sector(tdat[i,2], tdat[i,1], rou1 = 0.7, rou2 = 0.5, col = tdat[i,3])
}

draw.sector(num_t[1,2], num_t[1,1], rou1 = 0.8, rou2 = 0.75, col = num_t[1,3],lty=1,border=NA)

draw.sector(num_t[4,2], num_t[2,1], rou1 = 0.9, rou2 = 0.85, col = num_t[2,3],lty=2)
draw.sector(num_t[2,2], num_t[4,1], rou1 = 0.9, rou2 = 0.85, col = num_t[4,3],lty=2)
draw.sector(num_t[3,2], num_t[3,1], rou1 = 0.9, rou2 = 0.85, col = num_t[3,3],lty=1,border=NA)
draw.sector(num_t[5,2], num_t[5,1], rou1 = 0.9, rou2 = 0.85, col = num_t[5,3],lty=1,border=NA)

draw.sector(num_t[6,2], num_t[6,1], rou1 = 1, rou2 = 0.95, col = num_t[6,3],lty=1,border=NA)
draw.sector(num_t[7,2], num_t[7,1], rou1 = 1, rou2 = 0.95, col = num_t[7,3],lty=2)

draw.sector(num_t[8,2], num_t[8,1], rou1 = 0.8, rou2 = 0.75, col = num_t[8,3],lty=1,border=NA)

draw.sector(num_t[10,2], num_t[9,1], rou1 = 1.1, rou2 = 1.05, col = num_t[9,3],lty=2)
draw.sector(num_t[9,2], num_t[10,1], rou1 = 1.1, rou2 = 1.05, col = num_t[10,3],lty=2)

draw.sector(num_t[11,2], num_t[11,1], rou1 = 0.8, rou2 = 0.75, col = num_t[11,3],lty=1,border=NA)

draw.sector(num_t[12,2], num_t[13,1], rou1 = 1.2, rou2 = 1.15, col = num_t[12,3],lty=1,border=NA)

draw.sector(num_t[17,2], num_t[14,1], rou1 = 1.3, rou2 = 1.25, col = num_t[14,3],lty=2)
draw.sector(num_t[15,2], num_t[17,1], rou1 = 1.3, rou2 = 1.25, col = num_t[14,3],lty=2)
draw.sector(num_t[16,2], num_t[16,1], rou1 = 1.3, rou2 = 1.25, col = num_t[14,3],lty=1,border=NA)

draw.sector(num_t[18,2], num_t[18,1], rou1 = 0.8, rou2 = 0.75, col = num_t[18,3],lty=1,border=NA)

draw.sector(num_t[19,2], num_t[19,1], rou1 = 1.4, rou2 = 1.35, col = num_t[19,3],lty=1,border=NA)
draw.sector(num_t[20,2], num_t[20,1], rou1 = 1.4, rou2 = 1.35, col = num_t[20,3],lty=1,border=NA)
draw.sector(num_t[21,2], num_t[21,1], rou1 = 1.4, rou2 = 1.35, col = num_t[21,3],lty=1,border=NA)
draw.sector(num_t[22,2], num_t[22,1], rou1 = 1.4, rou2 = 1.35, col = num_t[22,3],lty=1,border=NA)

draw.sector(num_t[23,2], num_t[23,1], rou1 = 0.8, rou2 = 0.75, col = num_t[23,3],lty=1,border=NA)

draw.sector(num_t[24,2], num_t[24,1], rou1 = 1.4, rou2 = 1.35, col = num_t[24,3],lty=1,border=NA)
draw.sector(num_t[25,2], num_t[25,1], rou1 = 1.4, rou2 = 1.35, col = num_t[25,3],lty=1,border=NA)
draw.sector(num_t[26,2], num_t[26,1], rou1 = 1.4, rou2 = 1.35, col = num_t[26,3],lty=1,border=NA)
draw.sector(num_t[27,2], num_t[27,1], rou1 = 1.4, rou2 = 1.35, col = num_t[27,3],lty=1,border=NA)

draw.sector(num_t[28,2], num_t[28,1], rou1 = 0.8, rou2 = 0.75, col = num_t[28,3],lty=1,border=NA)

draw.sector(num_t[29,2], num_t[29,1], rou1 = 0.8, rou2 = 0.75, col = num_t[29,3],lty=1,border=NA)

draw.sector(num_t[31,2], num_t[32,1], rou1 = 1.1, rou2 = 1.05, col = num_t[33,3],lty=2)
draw.sector(num_t[32,2], num_t[31,1], rou1 = 1.1, rou2 = 1.05, col = num_t[33,3],lty=2)
draw.sector(num_t[30,2], num_t[30,1], rou1 = 1.1, rou2 = 1.05, col = num_t[33,3],lty=1,border=NA)

draw.sector(num_t[33,2], num_t[33,1], rou1 = 0.8, rou2 = 0.75, col = num_t[33,3],lty=1,border=NA)
draw.sector(num_t[34,2], num_t[34,1], rou1 = 0.8, rou2 = 0.75, col = num_t[33,3],lty=1,border=NA)

draw.sector(num_t[35,2], num_t[35,1], rou1 = 1.5, rou2 = 1.45, col = num_t[35,3],lty=1,border=NA)
draw.sector(num_t[36,2], num_t[36,1], rou1 = 1.5, rou2 = 1.45, col = num_t[36,3],lty=1,border=NA)
draw.sector(num_t[37,2], num_t[37,1], rou1 = 1.5, rou2 = 1.45, col = num_t[37,3],lty=1,border=NA)

draw.sector(num_t[40,2], num_t[38,1], rou1 = 1.6, rou2 = 1.55, col = num_t[38,3],lty=2)
draw.sector(num_t[38,2], num_t[40,1], rou1 = 1.6, rou2 = 1.55, col = num_t[40,3],lty=2)
draw.sector(num_t[39,2], num_t[39,1], rou1 = 1.6, rou2 = 1.55, col = num_t[39,3],lty=1,border=NA)

draw.sector(num_t[41,2], num_t[41,1], rou1 = 1.7, rou2 = 1.65, col = num_t[41,3],lty=1,border=NA)
draw.sector(num_t[42,2], num_t[42,1], rou1 = 1.7, rou2 = 1.65, col = num_t[42,3],lty=1,border=NA)
draw.sector(num_t[43,2], num_t[43,1], rou1 = 1.7, rou2 = 1.65, col = num_t[43,3],lty=1,border=NA)

draw.sector(num_t[46,2], num_t[44,1], rou1 = 1.8, rou2 = 1.75, col = num_t[44,3],lty=2)
draw.sector(num_t[44,2], num_t[46,1], rou1 = 1.8, rou2 = 1.75, col = num_t[46,3],lty=2)
draw.sector(num_t[45,2], num_t[45,1], rou1 = 1.8, rou2 = 1.75, col = num_t[45,3],lty=1,border=NA)

draw.sector(num_t[49,2], num_t[49,1], rou1 = 1.9, rou2 = 1.85, col = num_t[46,3],lty=1,border=NA)
draw.sector(num_t[47,2], num_t[47,1], rou1 = 1.9, rou2 = 1.85, col = num_t[47,3],lty=1,border=NA)
draw.sector(num_t[48,2], num_t[48,1], rou1 = 1.9, rou2 = 1.85, col = num_t[48,3],lty=1,border=NA)

draw.sector(num_t[52,2], num_t[50,1], rou1 = 2, rou2 = 1.95, col = num_t[49,3],lty=2)
draw.sector(num_t[50,2], num_t[52,1], rou1 = 2, rou2 = 1.95, col = num_t[50,3],lty=2)
draw.sector(num_t[51,2], num_t[51,1], rou1 = 2, rou2 = 1.95, col = num_t[51,3],lty=1,border=NA)

draw.sector(num_t[53,2], num_t[53,1], rou1 = 1, rou2 = 0.95, col = num_t[53,3],lty=1,border=NA)
draw.sector(num_t[54,2], num_t[54,1], rou1 = 1, rou2 = 0.95, col = num_t[54,3],lty=1,border=NA)

draw.sector(num_t[55,2], num_t[55,1], rou1 = 1.2, rou2 = 1.15, col = num_t[55,3],lty=1,border=NA)
draw.sector(num_t[56,2], num_t[56,1], rou1 = 1.2, rou2 = 1.15, col = num_t[56,3],lty=1,border=NA)

draw.sector(num_t[57,2], num_t[57,1], rou1 = 1.5, rou2 = 1.55, col = num_t[33,3],lty=1,border=NA)
dev.off()

tiff("Legend_fig1_Circ.tiff",width=7, height=7, units= 'in', res=600, pointsize=1/600)
plot(1,1,col="white")
legend("topright",c("1","2","3","6","7","8","9","10","11","12","14","Unplaced scaffolds"),
col=c(num_t[29,3],num_t[8,3],num_t[23,3],num_t[24,3],num_t[11,3],num_t[2,3],num_t[1,3],num_t[10,3],num_t[19,3],num_t[55,3],num_t[6,3],num_t[33,3]),
fill=c("#a6cee3","#1f78b4","#b2df8a","#33a02c","#fb9a99","#e31a1c","#fdbf6f","#ff7f00","#cab2d6","#6a3d9a","#ffff99","#b15928","blue2","hotpink2"),
title="AFS chromosome")
dev.off()
