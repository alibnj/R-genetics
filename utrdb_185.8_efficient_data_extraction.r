

m.185.8 <- read.csv('./185 (8) Match_Hits.csv', header=TRUE)
m.185.8$X <- NULL
m.185.8$rRNA <- NULL
m.185.8 <- unique(m.185.8)

#Removing the mRNAs that were found before:
p.10 <- read.csv('./185 mRNA List-Peak-(10).csv', header=TRUE)
np.10 <- read.csv('./185 mRNA List-Not Peak-(10).csv', header=TRUE)
m.185.8.rem <- as.data.frame(m.185.8[!(m.185.8$mRNA %in% p.10$t.185.10.mRNA),])
colnames(m.185.8.rem) <- c('mRNA')
m.185.8.rem <- as.data.frame(m.185.8.rem[!(m.185.8.rem$mRNA %in% np.10$mRNA),])
colnames(m.185.8.rem) <- c('Locus')
write.csv(m.185.8.rem, './m.185.8.rem.csv')

#Searching UTRdb
m.185.8.rem <- read.csv('./m.185.8.rem.csv', header=TRUE)
m.185.8.rem$X <- NULL
library(RSelenium)

for (j in 20:32)
  {
  checkForServer()
  startServer()
  remDrv <- remoteDriver()
  remDrv$open()
  C <- data.frame(Gene = character(), Organism = character(), Locus = character(), Region = character(), Length = numeric())
  filename <- paste("./UTRdb.185.8.rem", j, ".csv", sep = " ")
  for (i in ((j*1500)+1):((j+1)*1500))
  {
    remDrv$navigate('http://utrdb.ba.itb.cnr.it/search')
    remDrv$findElement(using = "xpath", "//select[@id = 'utr_db']/option[@value = '1']")$clickElement()
    remDrv$findElement(using = "xpath", "//input[@id ='ids']")$sendKeysToElement(list(m.185.8.rem[i,1], "\uE007"))
    doc <- htmlParse(remDrv$getPageSource()[[1]])
    tab <- readHTMLTable(doc)
    tab <- as.data.frame(tab)
    C <- rbind(C,tab)
    Sys.sleep(1)
  }
  write.csv(C, filename)
  remDrv$closeWindow()
  remDrv$quit()
  remDrv$closeServer()
  Sys.sleep(5)
}
#------------------------------------------ 
checkForServer()
startServer()
remDrv <- remoteDriver()
remDrv$open()
C <- data.frame(Gene = character(), Organism = character(), Locus = character(), Region = character(), Length = numeric())

for (i in 49501:49706)
{
  remDrv$navigate('http://utrdb.ba.itb.cnr.it/search')
  remDrv$findElement(using = "xpath", "//select[@id = 'utr_db']/option[@value = '1']")$clickElement()
  remDrv$findElement(using = "xpath", "//input[@id ='ids']")$sendKeysToElement(list(m.185.8.rem[i,1], "\uE007"))
  doc <- htmlParse(remDrv$getPageSource()[[1]])
  tab <- readHTMLTable(doc)
  tab <- as.data.frame(tab)
  C <- rbind(C,tab)
  Sys.sleep(1)
}
write.csv(C, './UTRdb.185.8.rem 33 .csv')
remDrv$closeWindow()
remDrv$quit()
remDrv$closeServer()
Sys.sleep(5)

# - - - - - - - - Binding Rows - - - - - - - - - - - - - - 

C <- read.csv('./UTRdb.185.8.rem 0 .csv')
for (j in 1:33)
{
  filename <- paste("./UTRdb.185.8.rem", j, ".csv", sep = " ")
  temp <- read.csv(filename, header=TRUE)
  C <- rbind(C, temp)
}
write.csv(C, './UTRdb.185.8.rem COMPLETE.csv')

# - - - - - - - - Adding missing data to rem file - - - - - - - - - - - - - - 

m.hit185.8 <- read.csv('./185 (8) Match_Hits.csv', header=TRUE)
m.hit185.8$X <- NULL
m.hit185.8$rRNA <- NULL
m.hit185.8 <- unique(m.hit185.8)

mr1 <- read.csv('./UTRdb.185.8.rem COMPLETE.csv', header=TRUE)
mr2 <- read.csv('./185 UTRdb-Not Peak-(10).csv', header=TRUE)
mr3 <- read.csv('./185 UTRdb-Peak-(10).csv', header=TRUE)

m.hit185.8$gene <- mr1$NULL.Gene[match(m.hit185.8$mRNA, mr1$NULL.Locus)]
m1 <- m.hit185.8[complete.cases(m.hit185.8$gene),]
m.hit185.8$gene <- mr2$Gene[match(m.hit185.8$mRNA, mr2$Locus)]
m2 <- m.hit185.8[complete.cases(m.hit185.8$gene),]
m.hit185.8$gene <- mr3$NULL.Gene[match(m.hit185.8$mRNA, mr3$NULL.Locus)]
m3 <- m.hit185.8[complete.cases(m.hit185.8$gene),]

UTR.185.8 <- rbind(m1,m2)
UTR.185.8 <- rbind(UTR.185.8,m3)

write.csv(UTR.185.8, './UTRdb.185(8).complete.csv')
