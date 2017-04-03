setwd("/Users/pracz/Work/NZILBB/maori")
load("maori_functions.rda")

library(ndl)
library(plyr)
library(stringr)
library(reshape2)
library(tidyr)

# files
cfeat <- read.delim("~/cons-matrix.txt", sep="\t")
vfeat <- read.delim("~/vowels-matrix.txt", sep="\t")
baselist <- read.delim("/maori_dict_head_deriv_4_cv.txt", sep=",", header=F)
# todo
# need to test outputs
########################################################################################################################
# list of functions
# filter1										-- takes baselist as input and filters it  bit
# filter1NoInitRedupWords		-- takes baselist as input and removes words that have an initial repeated CV
# filter1NoRedup						-- takes baselist as input and removes repeated sequences of any length (but leaves the rest of the word)

# filter2ForC								-- takes output of filter1, removes vowels so you can then do cons digrams
# filter2ForV								-- takes output of filter1, removes cons so you can then do vowel digrams

# allDigrams								-- takes output of filter2, creates digrams from whatever you hand it to, calculates o, e, o/e
# posDigrams								-- takes output of filter2 and output of alldigrams, creates positional digrams from c/v, calculates o, e, o/e. obviously if you are doing consonants, you will need the base list filtered to have consonants as well as consonant digrams (this script is a massive hack and will make R scream bloody murder but it works)
# posUnigrams								-- takes output of filter2, counts segments in different positions

# cPlace										-- takes output of a digram maker, adds place features for consonants
# cPlace2										-- exactly like cPlace, except this one takes the output of a unigram maker

# posRatios									-- takes output of a positional digram maker, calculates positional o over e and o minus e (which are moderately useful)

# vPlace										-- same as cPlace but for vowels
# vPlace2										-- same as cPlace2 but for vowels


########################################################################################################################
# so how do you make the datasets?
# all consonantal digrams:
# dat1 <- filter1(baselist)
# dat1 <- filter2ForC(dat1)
# dat1 <- allDigrams(dat1)
# dat1 <- cPlace(dat1,cfeat)
# consonantal digrams from non-repeating chunks:
# dat2 <- filter1NoRedup(baselist)
# dat2 <- filter2ForC(dat2)
# dat2 <- allDigrams(dat2)
# dat2 <- cPlace(dat2,cfeat)
# positional consonantal digrams:
# dat3 <- filter1(baselist)
# dat3 <- filter2ForC(dat3)
# dat3 <- posDigrams(dat3,dat1)
# dat3 <- cPlace(dat3,cfeat)
# dat3 <- posRatios(dat3)
# consonant unigrams per position:
# dat4 <- filter1(baselist)
# dat4 <- filter2ForC(dat4)
# dat4 <- posUnigrams(dat4)
# dat4 <- cPlace2(dat4,cfeat)
# that's all folks!

########################################################################################################################
########################################################################################################################
# get base lists

########################################################################################################################
# filter list 


filter1 <- function(dat){
# 1 data from cvcv words
dat <- dat[,1]
dat <- as.character(dat)

# cleaning up. 

# getting rid of first three rows which have QUESTION MARKS for some reason
dat <- dat[4:8950]
return(dat)}

########################################################################################################################
# filter list to remove reduplicated chunks

filter1NoRedup <- function(dat){
	
# 1 data from cvcv words
dat <- dat[,1]
dat <- as.character(dat)

# replace anythingx2 by _

dat2 <- sub("([WhknmNprtw][a-zA-z]{1,})\\1", "_", dat) 
dat2 <- sub("([WhknmNprtw][a-zA-z]{1,})\\1", "_", dat2) 
dat2 <- sub("([WhknmNprtw][a-zA-z]{1,})\\1", "_", dat2) 
# dat2 <- sub("([WhknmNprtw][a-zA-z]{1,})\\1", "_", dat2)  # makes no diff

dat3 <- strsplit(dat2, "_")
dat3 <- unlist(dat3)
# we end up with 9724 bits, some of them no doubt empty
dat3 <- dat3[dat3!=""]
# final count: 8473 bits
return(dat3)
}

# modified version: throws out the word if it has a repeated anything the size of or larger than a cv syllable

filter1NoRedupKillWord <- function(dat){
	
# 1 data from cvcv words
dat <- dat[,1]
dat <- as.character(dat)

# kill words with anythingx2

dat2 <- grep("([WhknmNprtw][a-zA-z]{1,})\\1", dat, invert=T, value=T) 

# final count: 5820 bits
return(dat2)
}



########################################################################################################################
# filter list to remove initial repeated CV


filter1NoInitRedupWords <- function(dat){
	
# 1 data from cvcv words
dat <- dat[,1]
dat <- as.character(dat)

# drop ^cvcv

hits <- !grepl("^(.[aeiouAEIOU])\\1", dat, perl=T)
dat <- cbind(dat,hits)
dat <- subset(dat, hits==T)
dat <- dat[,1]

# original: 8950 words. w/o initial reduplication: 8357 words.

# cleaning up. 

# getting rid of first three rows which have QUESTION MARKS for some reason
dat <- dat[4:8950]
return(dat)
}

########################################################################################################################
########################################################################################################################
# narrow it down

########################################################################################################################
# narrow down to consonants

filter2ForC <- function(dat){
# getting rid of vowels (since we want to have CC transitional probabilities)
for (i in 1:10){
dat <- sub("[aeoiuAEIOU]","", dat)
}

for (i in 1:10){
dat <- sub("W","f", dat)
}
return(dat)
}

########################################################################################################################
# narrow down to vowels

filter2ForV <- function(dat){

# getting rid of CONSONANTS
for (i in 1:10){
dat <- sub("[vbfjsglptkhmnNWwr?]","", dat) # all sorts of weird characters in there
}

return(dat)
}


########################################################################################################################
########################################################################################################################
# make digram lists

########################################################################################################################
# digrams from all positions

# input: output of a filter function for consonants
# same for vowels
allDigrams <- function(baselist){
# building digrams using Harald Baayen's digrammifier
digrams <- as.data.frame(orthoCoding(baselist, g=2))
# getting TP-s across the board
digrams <- as.character(digrams[,1])
# split it across "_"
digrams <- strsplit(digrams, "_")
# now it's a list
digrams <- unlist(digrams)
# there we go.

# we don't care about `#'

digrams <- sub("#.","",digrams)
digrams <- sub(".#","",digrams)
digrams <- digrams[digrams != ""]

# character frequencies NOT from bigrams
letters2 <-  as.data.frame(orthoCoding(baselist, g=1))
letters2 <- as.character(letters2[,1])
letters2 <- strsplit(letters2, "_")
letters2 <- unlist(letters2)
letters2 <- count(letters2)

# digram frequencies
digrams <- count(digrams)
# total frequency of every character
tot.freq.uni <- sum(letters2$freq)
# total frequency of every digram
tot.freq.bi <- sum(digrams$freq)

# calculate relative frequency for each character
letters2$rel.freq <- letters2$freq/tot.freq.uni

names(digrams) <- c("digram","observed")
# get first segment for digrams
digrams$seg1 <- substr(digrams$digram,1,1)
# get second segment for digrams
digrams$seg2 <- substr(digrams$digram,2,2)
# get relative freq for seg1
digrams <- merge(digrams,letters2,by.x="seg1",by.y="x")
names(digrams) <- c("seg1","digram","observed","seg2","seg1.freq","seg1.rel.freq")
# get relative freq for seg2
digrams <- merge(digrams,letters2,by.x="seg2",by.y="x")
names(digrams) <- c("seg2", "seg1","digram","observed","seg1.freq","seg1.rel.freq","seg2.freq","seg2.rel.freq")
# get expected frequency
digrams$expected <- digrams$seg1.rel.freq * digrams$seg2.rel.freq * tot.freq.bi
# get observed over expected
digrams$o.over.e <- digrams$observed/digrams$expected

# get rid of hapaxes
digrams <- subset(digrams, observed>1)

# groovy

return(digrams)
}

### JEN's VERSION WHICH CALCULATES EXPECTED IN POSITION SPECIFIC WAY, AS PER PIERREHUMBERT AND OTHERS....
allDigramsJen <- function(baselist){
  # building digrams using Harald Baayen's digrammifier
  digrams <- as.data.frame(orthoCoding(baselist, g=2))
 
  # getting TP-s across the board
  digrams <- as.character(digrams[,1])
  # split it across "_"
  digrams <- strsplit(digrams, "_")
  # now it's a  list
  digrams <- unlist(digrams)
  # there we go.
 
  # we don't care about `#'
 
  digrams <- sub("#.","",digrams)
  digrams <- sub(".#","",digrams)
  digrams <- digrams[digrams != ""]
 
  
#   letters2 <-  as.data.frame(orthoCoding(baselist, g=1))
#    letters2 <- as.character(letters2[,1])
#    letters2 <- strsplit(letters2, "_")
#    letters2 <- unlist(letters2)
#    letters2 <- count(letters2)
#  
  # do this instead:
  # digram frequencies
  digrams <- count(digrams)
 
  names(digrams) <- c("digram","observed")
  # get first segment for digrams
  digrams$seg1 <- substr(digrams$digram,1,1)
  # get second segment for digrams
  digrams$seg2 <- substr(digrams$digram,2,2)
  # get relative freq for seg1
  seg1freq= aggregate(observed ~ seg1, sum, dat=digrams)
  names(seg1freq) <- c("seg1", "seg1o")
  seg2freq = aggregate(observed ~ seg2, sum, dat=digrams)
  names(seg2freq) <- c("seg2", "seg2o")
  digramsa = merge(digrams, seg1freq, by.x="seg1", by.y="seg1")
  digramsb = merge(digramsa, seg2freq, by.x="seg2", by.y="seg2")
 
  digramsb$expected =digramsb$seg1o * digramsb$seg2o / sum(digramsb$observed)
  digramsb$o.over.e = digramsb$observed/digramsb$expected
 
  # get rid of hapaxes
  digramsb <- subset(digramsb, observed>1)
 
  return(digramsb)
}

########################################################################################################################
# digrams w/ pos encoded

# this is a massive hack. so R will complain about the unspeakable things it makes it do.
# NOT sensitive to max word length in the baselist, ha

# input: output of a filter function for consonants (preferably the vanilla one), output of allConDigrams,
# same for vowels
posDigrams <- function(baselist,digrams){
	
baselist <- baselist[!is.element(baselist, "")]
posgrams <- as.data.frame(orthoCoding(baselist, g=2))
posgrams <- sub("^...", "", as.character(posgrams[,1]))
posgrams <- sub("...$", "", posgrams)
# the bits that were too short are still there, like `##' and `h#'
posgrams <- posgrams[nchar(posgrams) > 2]
# getting rid of #. this effectively gets rid of words that have less than 3 consonants in them. but we want to mark digram position not cons position (we have the unigram positional thing for that.)
posgrams <- strsplit(posgrams, "_") # sfsg
# so here, when you mark positions, for the word xyz the first ngram will be xy, and then both x and y are 1nd. but yz is the second ngram, so y and z are 2nd.
# turn this LIST into a df which remembers the position of digrams within the words

# http://www.r-bloggers.com/converting-a-list-to-a-data-frame/
# install.packages("devtools")
require(devtools)
source_gist(4676064)
bumm <- as.data.frame(posgrams)
# The length of vectors are not the same and do not are not named, the results may not be correct!
names(bumm) <- paste("pos", 1:ncol(bumm),sep="")
bumm2 <- gather(bumm)
names(bumm2) <- c("position","digram")
posgrams <- bumm2
posgrams$digram <- as.character(posgrams$digram)
posgrams.list <- split(posgrams, posgrams$position)	
posgrams.list2 <- as.list(NULL)

for (i in 1:length(posgrams.list)){
this.pos <- posgrams.list[[i]]
c.pos <- count(this.pos, var="digram")
names(c.pos)[2] <- "observed"
c.pos$position <- paste("pos",i,sep="")
posgrams.list2[[i]] <- c.pos
}

c.posgrams <- do.call("rbind",posgrams.list2)

# it has all the NA-s now since I turned it into char above
c.posgrams <- subset(c.posgrams, digram != "<NA>")
####################################################################################
# c.posgrams: 509 rows
digrams.merge <- digrams[,c("seg1","seg2","digram","expected")]
digrams.merge$digram <- as.character(digrams.merge$digram)
c.posgrams[!(c.posgrams[,1] %in% digrams.merge$digram),]
mega.digrams <- merge(c.posgrams,digrams.merge)
# this is all sorts of weird stuff that's not actually Maori segments, like d and l and g
mega.digrams$num.position <- substr(mega.digrams$position,4,4)
mega.digrams$position <- as.numeric(mega.digrams$num.position)
mega.digrams$num.position <- NULL
mega.digrams$o.over.e <- mega.digrams$observed/mega.digrams$expected
return(mega.digrams)
}

########################################################################################################################
# unigrams across position

# input: output of a filter function for consonants (preferably the vanilla one)
# same for vowels
posUnigrams <- function(baselist){
pos.list <- as.list(NULL)	

for ( i in 1:max(nchar(baselist)) ){
pos <- substr(baselist, i,i)
pos.list[[i]] <- pos
}	

pos.list2 <- as.list(NULL)	
# count things in positions
for (i in 1:length(pos.list)){	
c.pos <- count(pos.list[[i]])
names(c.pos) <- c("segment", "observed")
c.pos$position <- paste("pos",i,sep="")
pos.list2[[i]] <- c.pos
}
c.posunigrams <- do.call("rbind",pos.list2)

return(c.posunigrams)
}


########################################################################################################################
########################################################################################################################
# adding place and stuff

########################################################################################################################
# place features (consonants)

cPlace <- function(dat,cfeatures){
# get features for sounds

# consonants
cfeatures[cfeatures$X=="r",]$labial <- F
cfeatures[cfeatures$X=="r",]$bilabial <- F

c.seg1 <- cfeatures
names(c.seg1) <- c("seg1","seg1.labial","seg1.coronal","seg1.dorsal","seg1.bilabial","seg1.dental","seg1.alveolar","seg1.palatal","seg1.velar","seg1.obstruent","seg1.sonorant","seg1.stop","seg1.continuant","seg1.glide","seg1.consonantal","seg1.oral","seg1.affricate","seg1.strident","seg1.distributed","seg1.lateral","seg1.rhotic","seg1.nasal","seg1.voice","seg1.voiceless","seg1.spread.glottis")
c.seg2 <- cfeatures
names(c.seg2) <- c("seg2","seg2.labial","seg2.coronal","seg2.dorsal","seg2.bilabial","seg2.dental","seg2.alveolar","seg2.palatal","seg2.velar","seg2.obstruent","seg2.sonorant","seg2.stop","seg2.continuant","seg2.glide","seg2.consonantal","seg2.oral","seg2.affricate","seg2.strident","seg2.distributed","seg2.lateral","seg2.rhotic","seg2.nasal","seg2.voice","seg2.voiceless","seg2.spread.glottis")
#######################################################
# go back for digrams

dat <- merge(dat,c.seg1)
dat <- merge(dat,c.seg2)

dat$has.nasal <- ifelse(dat$seg1.nasal==T | dat$seg2.nasal==T, T, F)
dat$has.stop <- ifelse(dat$seg1.stop==T | dat$seg2.stop==T, T, F)
dat$has.labial <- ifelse(dat$seg1.labial==T | dat$seg2.labial==T, T, F)
dat$has.coronal <- ifelse(dat$seg1.coronal==T | dat$seg2.coronal==T, T, F)
dat$has.velar <- ifelse(dat$seg1.velar==T | dat$seg2.velar==T, T, F)
dat$place <- "diff.place"
# if same place
dat$place <- ifelse( ( dat$seg1.labial==T & dat$seg2.labial==T ) | ( dat$seg1.coronal==T & dat$seg2.coronal==T )  | ( dat$seg1.velar==T & dat$seg2.velar==T ), "same.place", dat$place)
# if identical
dat$place <- ifelse(as.character(dat$seg1)==as.character(dat$seg2),"ident",dat$place)
dat$place <- as.factor(dat$place)
dat$o.minus.e <- dat$observed - dat$expected

dat$seg1.place<-NA
dat[dat$seg1.labial==T,]$seg1.place<-"labial"
dat[dat$seg1.coronal==T,]$seg1.place<-"coronal"
dat[dat$seg1.velar==T,]$seg1.place<-"velar"
dat$seg2.place<-NA
dat[dat$seg2.labial==T,]$seg2.place<-"labial"
dat[dat$seg2.coronal==T,]$seg2.place<-"coronal"
dat[dat$seg2.velar==T,]$seg2.place<-"velar"	
	
return(dat)
}

########################################################################################################################
# place features vowels

vPlace <- function(dat,vfeatures){

# adding features
v.seg1 <- vfeatures
names(v.seg1) <- c("seg1","seg1.long",	"seg1.front",	"seg1.back",	"seg1.close",	"seg1.open",	"seg1.round")
v.seg2 <- vfeatures
names(v.seg2) <- c("seg2","seg2.long",	"seg2.front",	"seg2.back",	"seg2.close",	"seg2.open",	"seg2.round")

dat <- merge(dat,v.seg1)
dat <- merge(dat,v.seg2)
	
	
dat$has.front <- ifelse(dat$seg1.front==T | dat$seg2.front==T,T,F)
dat$has.back <- ifelse(dat$seg1.back==T | dat$seg2.back==T,T,F)
dat$has.open <- ifelse(dat$seg1.open==T | dat$seg2.open==T,T,F)
dat$has.close <- ifelse(dat$seg1.close==T | dat$seg2.close==T,T,F)
dat$has.long <- ifelse(dat$seg1.long==T | dat$seg2.long==T,T,F)

dat$ident <- ifelse(as.character(dat$seg1)==as.character(dat$seg2),T,F) 
return(dat)
}


########################################################################################################################
# this one's for the output of the unigram pos thing

cPlace2 <- function(dat,cfeatures){
dat2 <- dat	
names(cfeatures) <- c("seg1","labial","coronal","dorsal","bilabial","dental","alveolar","palatal","velar","obstruent","sonorant","stop","continuant","glide","consonantal","oral","affricate","strident","distributed","lateral","rhotic","nasal","voice","voiceless","spread.glottis")
names(cfeatures)[1] <- "segment"
cfeatures$segment <- as.character(cfeatures$segment)
dat2$segment <- as.character(dat2$segment)
dat2 <- merge(dat2,cfeatures)
dat2$num.position <- substr(dat2$position,4,4)
dat2$position <- as.numeric(dat2$num.position)
dat2$num.position <- NULL
return(dat2)
}

vPlace2 <- function(dat,vfeatures){
names(vfeatures)[1] <- "segment"
dat <- merge(dat,vfeatures)
dat$num.position <- substr(dat$position,4,4)
dat$position <- as.numeric(dat$num.position)
dat$num.position <- NULL
return(dat)
}


########################################################################################################################
# pos o/e ratios


posRatios <- function(dat){
# I need the frequency of the two segments per position. well I can get it from the digrams in that position
seg1fr <- aggregate(observed ~ seg1 + position, dat, sum)
seg2fr <- aggregate(observed ~ seg2 + position, dat, sum)
# and then divide that by the total frequency of letters in that position (actually, in the digram in that position)
segsum <- aggregate(observed ~ position, dat, sum)
names(segsum)[2] <- "pos.total"
seg1fr <- merge(segsum,seg1fr)
seg2fr <- merge(segsum,seg2fr)
seg1fr$seg1.pos.rel.freq <- seg1fr$observed/seg1fr$pos.total
seg2fr$seg2.pos.rel.freq <- seg2fr$observed/seg2fr$pos.total
head(seg1fr,1)
seg1fr <- seg1fr[,c("seg1","position","seg1.pos.rel.freq")]
seg2fr <- seg2fr[,c("seg2","position","seg2.pos.rel.freq")]
dat <- merge(dat,seg1fr)
dat <- merge(dat,seg2fr)
dat <- merge(dat,segsum)

dat$pos.expected <- dat$seg1.pos.rel.freq * dat$seg2.pos.rel.freq * dat$pos.total

dat$pos.o.over.e <- dat$observed / dat$pos.expected
dat$pos.o.minus.e <- dat$observed - dat$pos.expected
return(dat)
}



# save.image("maori_functions.rda")
