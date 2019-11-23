#-----------------------------------------------------------
#--------READING AND CLEANING-------------------------->>>>>
#-----------------------------------------------------------

GAF <- readLines('C:/Users/nslavov/Desktop/GSDB/goa_mouse.gaf')
GAF <- GAF[!startsWith(GAF, "!")]
writeLines(GAF, con = 'C:/Users/nslavov/Desktop/GSDB/goa_mouse-only.table.txt')
GAF <- read.delim('C:/Users/nslavov/Desktop/GSDB/goa_mouse-only.table.txt', header = FALSE, stringsAsFactors = FALSE, sep = '\t')
GAF <- GAF[c('V3', 'V5')]
colnames(GAF) <- c('gene', 'GO')

go.names <- read.delim('C:/Users/nslavov/Desktop/GSDB/GO.names.csv', header = TRUE, sep = ' ')
GAF$description <- go.names$name[match(GAF$GO, go.names$identifier)]
GAF$description <- paste("GO", GAF$description, sep = "_")
GAF$description <- gsub(" ", "_", GAF$description)

#-----------------------------------------------------------
#--------CONVERTING------------------------------------>>>>>
#-----------------------------------------------------------

GOs <- as.character(sort(unique(GAF$description)))
dimension <- aggregate(gene~GO,GAF,length)
GMT <- data.frame(matrix(ncol=max(dimension$gene)))

for (i in 1:length(GOs)) {
  GMT[i,1] <- GOs[i]
  GMT[i,2] <- 'NA'
  D <- unique(with(GAF, gene[description==GOs[i]]))
  GMT[i,3:(2+length(D))] <- D
}

write.table(GMT, 'C:/Users/nslavov/Desktop/GSDB/goa_mouse.complete.gmt', na="", row.names = FALSE, col.names = FALSE, sep = '\t')
#write.table(GOs, 'C:/Users/nslavov/Desktop/GSDB/GOs.txt', na="", row.names = FALSE, col.names = FALSE, sep = '\t')
