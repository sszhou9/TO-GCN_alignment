#TO-GCN alignment
library(parallel)
library(dplyr)
uvb <- read.table("UV_B_TO_GCN.txt")
uvc <- read.table("UV_C_TO_GCN.txt")
blast <- read.table("UV_B_UV_C_blast.txt",header=F)
bait <- read.table("bait.txt")
nbait <- nrow(bait)
auvb <- read.table("uvb_array.txt")
auvc <- read.table("uvc_array.txt")
auvb <- arrange(auvb,desc(V2))
auvc <- arrange(auvc,desc(V2))
p <- c("UV_B_UV_C_network_alignment")
p1 <- c("UV_B_UV_C_network_alignment")
for(baiti in 1:nbait){
uvbx <- uvb[!uvb[,1] %in% bait[,1],]
uvcx <- uvc[!uvc[,1] %in% bait[,1],]
b <- as.character(bait[baiti,][1])
bb <- auvb[auvb[,1] %in% b,][1,][2]
cc <- auvc[auvc[,1] %in% b,][1,][2]
bbb <- filter(auvb,V2 == (as.character(bb))|V2 == (as.character(bb+1))|V2 == (as.character(bb+2)))
ccc <- filter(auvc,V2 == (as.character(cc))|V2 == (as.character(cc+1))|V2 == (as.character(cc+2)))
uvb1 <- uvb[(uvb[,1] %in% b)&(uvb[,2] %in% bbb[,1]),]
uvc1 <- uvc[(uvc[,1] %in% b)&(uvc[,2] %in% ccc[,1]),]
nuvb1 <- nrow(uvb1)
nuvc1 <- nrow(uvc1)
for(nuvb1i in 1:nuvb1){
  for(nuvc1i in 1:nuvc1){
b1 <- as.character(uvb1[nuvb1i,][2])
c1 <- as.character(uvc1[nuvc1i,][2])  
blast1 <- filter(blast,((V1 == b1&V2 == c1)|(V2 == b1&V1 == c1)))
if(b1 == c1|!is.na(blast1[1,][1]))
{
uvb2 <- filter(uvbx,V2 == b1)
uvc2 <- filter(uvcx,V2 == c1)
uvbx <- uvbx[!uvbx[,2] %in% b1,]
uvcx <- uvcx[!uvcx[,2] %in% c1,]
bb1 <- auvb[auvb[,1] %in% b1,][1,][2]
cc1 <- auvc[auvc[,1] %in% c1,][1,][2]
bbb1 <- filter(auvb,V2 == (as.character(bb1))|V2 == (as.character(bb1+1))|V2 == (as.character(bb1+2)))
ccc1 <- filter(auvc,V2 == (as.character(cc1))|V2 == (as.character(cc1+1))|V2 == (as.character(cc1+2)))
uvb2 <- uvb2[(uvb2[,2] %in% b1)&(uvb2[,1] %in% bbb1[,1]),]
uvc2 <- uvc2[(uvc2[,2] %in% c1)&(uvc2[,1] %in% ccc1[,1]),]
nuvb2 <- nrow(uvb2)
nuvc2 <- nrow(uvc2)
for(nuvb2i in 1:nuvb2){
  for(nuvc2i in 1:nuvc2){
b2 <- as.character(uvb2[nuvb2i,][1])
c2 <- as.character(uvc2[nuvc2i,][1])
if(!is.na(b2)&!is.na(c2))
{
blast2 <- filter(blast,((V1 == b2&V2 == c2)|(V2 == b2&V1 == c2)))
if(b2 == c2|!is.na(blast2[1,][1]))
{
uvb3 <- filter(uvbx,V1 == b2)
uvc3 <- filter(uvcx,V1 == c2)
uvbx <- uvbx[!uvbx[,1] %in% b2,]
uvcx <- uvcx[!uvcx[,1] %in% c2,]
bb2 <- auvb[auvb[,1] %in% b2,][1,][2]
cc2 <- auvc[auvc[,1] %in% c2,][1,][2]
bbb2 <- filter(auvb,V2 == (as.character(bb2))|V2 == (as.character(bb2+1))|V2 == (as.character(bb2+2)))
ccc2 <- filter(auvc,V2 == (as.character(cc2))|V2 == (as.character(cc2+1))|V2 == (as.character(cc2+2)))
uvb3 <- uvb3[(uvb3[,2] %in% bbb2[,1]),]
uvc3 <- uvc3[(uvc3[,2] %in% ccc2[,1]),]
nuvb3 <- nrow(uvb3)
nuvc3 <- nrow(uvc3)
for(nuvb3i in 1:nuvb3){
  for(nuvc3i in 1:nuvc3){
b3 <- as.character(uvb3[nuvb3i,][2])
c3 <- as.character(uvc3[nuvc3i,][2])
if(!is.na(b3)&!is.na(c3))
{
blast3 <- filter(blast,((V1 == b3&V2 == c3)|(V2 == b3&V1 == c3)))
if(b3 == c3|!is.na(blast3[1,][1]))
{
m <- paste0("UV-B_",b,",UV-B_",b1,",UV-B_",b2,",UV-B_",b3,",UV-C_",b,",UV-C_",c1,",UV-C_",c2,",UV-C_",c3)
print(m)
p <- c(p,m)
}}}}}}}}}}}}
write.table(p,file="network_alignment.csv",quote = F,row.names = F,col.names = F)
