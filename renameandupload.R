# script to rename raw EEG files to preserve participant anonymity
# also re-orders by AQ score, and uploads data to OSF
# DHB 23/5/22

packagelist <- c('remotes','tictoc','utils')
missingpackages <- packagelist[!packagelist %in% installed.packages()[,1]]
if (length(missingpackages)>0){install.packages(missingpackages)}
if (packagelist[!'osfr' %in% installed.packages()[,1]]){remotes::install_github("centerforopenscience/osfr")}
packagelist <- c(packagelist,'osfr')
toinstall <- packagelist[which(!packagelist %in% (.packages()))]
invisible(lapply(toinstall,library,character.only=TRUE))


osf_auth(token = '3dEYuhZNmwWbG3xuhRIiVRo0T2oniOkEP3Ip8i1LPG4PNjeTRln54eGNAG7oyTT9xozWwJ')
osfproject <- osf_retrieve_node("ab3yv")
componentlist <- osf_ls_nodes(osfproject)
EEGID <- pmatch('Raw EEG data 1',as.character(unlist(componentlist[,1])))
EEGtokenLow <- componentlist[EEGID,2]
EEGfilesLow <- osf_ls_files(EEGtokenLow,n_max=300)
EEGID <- pmatch('Raw EEG data 2',as.character(unlist(componentlist[,1])))
EEGtokenHi <- componentlist[EEGID,2]
EEGfilesHi <- osf_ls_files(EEGtokenHi,n_max=300)

rawdata <- read.csv('~/Google Drive/Current work/Teaching/BSc projects 2021/ASD/BSc project questionnaires 2021 (Responses) - Form Responses 1.csv')

IDnos <- readMat('~/Google Drive/Current work/Research/normalization/IDnos.mat')
uniqueID <- unlist(IDnos$allU)
load('~/Google Drive/Current work/Research/normalization/questionnairedata.RData')
QIDs <- as.character(qdata$ID)

# fix duplicate ID numbers
uniqueID[12] <- '2001b'
QIDs[12] <- '2001b'
uniqueID[46] <- '2001c'
QIDs[46] <- '2001c'
uniqueID[56] <- '1357b'
QIDs[55] <- '1357b'
uniqueID[70] <- '1234b'
QIDs[68] <- '1234b'
uniqueID[77] <- '1234c'
QIDs[75] <- '1234c'
uniqueID[85] <- '2705b'
QIDs[83] <- '2705b'
QIDs[53] <- '5417'

IDmatching <- data.frame(uniqueID,QIDs)

# these are the indices in the questionnaire data for each EEG dataset
Qindices <- NULL
for (n in 1:nrow(IDmatching)){
  Qindices[n] <- match(IDmatching[n,1],IDmatching[,2])
}

Eindices <- which(!is.na(Qindices))
orderedAQ <- qdata$AQ[Qindices]

subjorder <- sort(orderedAQ,index.return=TRUE)


for (s in 1:length(Eindices)){
  if (s!=57){    # subject 57 did not complete the experiment so is omitted
   
  thisIDnumber <- uniqueID[s]  
  QIDnumber <- match(thisIDnumber,QIDs)  
  qresponses <- rawdata[QIDnumber,]
  qscores <- qdata[QIDnumber,]
  
  targetparticipant <- subjorder$ix[s]
  
  
  thisAQ <- orderedAQ[targetparticipant]
  qdata$AQ[Eindices[targetparticipant]]
  IDmatching[targetparticipant,]
  
  
  
  }
}





