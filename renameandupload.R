# script to rename raw EEG files to preserve participant anonymity
# also re-orders by AQ score, and uploads data to OSF
# DHB 23/5/22

packagelist <- c('remotes','tictoc','utils','R.matlab')
missingpackages <- packagelist[!packagelist %in% installed.packages()[,1]]
if (length(missingpackages)>0){install.packages(missingpackages)}
if (packagelist[!'osfr' %in% installed.packages()[,1]]){remotes::install_github("centerforopenscience/osfr")}
packagelist <- c(packagelist,'osfr')
toinstall <- packagelist[which(!packagelist %in% (.packages()))]
invisible(lapply(toinstall,library,character.only=TRUE))

my_wd <- getwd()

osf_auth(token = '')
osfproject <- osf_retrieve_node("ab3yv")
componentlist <- osf_ls_nodes(osfproject)
EEGID <- pmatch('Raw EEG data 1',as.character(unlist(componentlist[,1])))
EEGtokenLow <- componentlist[EEGID,2]
EEGfilesLow <- osf_ls_files(EEGtokenLow,n_max=300)
EEGID <- pmatch('Raw EEG data 2',as.character(unlist(componentlist[,1])))
EEGtokenHi <- componentlist[EEGID,2]
EEGfilesHi <- osf_ls_files(EEGtokenHi,n_max=300)

rawdatadir <- '~/Desktop/UG2021/'

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

subjcount <- 0
AQinsubjorder <- NULL
for (s in 1:length(uniqueID)){
  if (s!=58){    # subject 58 did not complete the questionnaire so is omitted
  subjcount <- subjcount + 1
  thisIDnumber <- uniqueID[s]  
  QIDnumber <- match(thisIDnumber,QIDs) 
  qscores <- qdata[QIDnumber,]
  AQinsubjorder[subjcount] <- qscores$AQ
  }
}

newsubjorder <- sort(AQinsubjorder,index.return=TRUE)

subjcount <- 0
for (s in 1:length(uniqueID)){
  if (s!=58){    # subject 58 did not complete the experiment so is omitted
      print(s)
      
      setwd(my_wd)
      subjcount <- subjcount + 1
      thisIDnumber <- uniqueID[s]  
      QIDnumber <- match(thisIDnumber,QIDs)  
      qresponses <- rawdata[QIDnumber,3:148]
      
      newpno <- which(newsubjorder$ix==subjcount)
      
      csvfiles <- dir(paste0(rawdatadir,'P',s+200),pattern='csv.gz',full.names=TRUE)
      
      newdir <- paste0('temp/raw/S',newpno)
      dir.create(newdir)
      for (f in 1:length(csvfiles)){
      file.copy(csvfiles[f],paste0(newdir,'/Block',f,'.csv.gz'))
      }
      write.csv(qresponses,file=paste0(newdir,'/questionnairedata.csv'),row.names=FALSE)
      
      setwd('temp/raw/')
      d <- dir(paste0('S',newpno))
      zip(paste0('S',newpno,'.zip'),files=paste0('S',newpno,'/',d))
      setwd(my_wd)
      
  }
}

# upload first 49 participants to the low AQ folder
for (s in 1:49){
  print(s)
  osf_upload(EEGtokenLow,paste0('temp/raw/S',s,'.zip'),progress=TRUE)
}

# upload remaining 51 participants to the high AQ folder
for (s in 50:100){
  print(s)
  osf_upload(EEGtokenHi,paste0('temp/raw/S',s,'.zip'),progress=TRUE)
}

for (s in 1:100){
  print(s)
  file.remove(paste0('temp/raw/S',s,'.zip'))
}
