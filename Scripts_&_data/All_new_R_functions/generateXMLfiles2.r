generateXMLfiles2 = function(templateFile, newXMLFile, prefix, fastaFile, jitter) {

	fasta = scan(file=fastaFile, what="", sep="\n", quiet=TRUE)
	sequenceIDs = c(); sequences = c()
	for (i in 1:length(fasta))
		{
			if (grepl(">",fasta[i])) sequenceIDs = c(sequenceIDs, gsub(">","",fasta[i]))
			if (!grepl(">",fasta[i])) sequences = c(sequences, fasta[i])
		}
	metadata = matrix(nrow=length(sequenceIDs), ncol=3); row.names(metadata) = sequenceIDs
	colnames(metadata) = c("collection_date","latitude","longitude")
	for (i in 1:dim(metadata)[1])
		{
			sequenceID = unlist(strsplit(sequenceIDs[i],"_"))
			metadata[i,"collection_date"] = sequenceID[2]
			metadata[i,"latitude"] = sequenceID[3]
			metadata[i,"longitude"] = sequenceID[4]
		}
	template = scan(templateFile, what="", sep="\n", quiet=T, blank.lines.skip=F)
	sink(file=newXMLFile); clusters = list(); clusters[[1]] = metadata
	for (i in 1:length(template))
		{
			if (grepl("TEMPLATE",template[i]))
				{
					template[i] = gsub("TEMPLATE", prefix, template[i])
				}
			if (grepl("JITTER",template[i]))
				{
					template[i] = gsub("JITTER", jitter, template[i])
				}
			cat(template[i],"\n")
			if (grepl("Insert taxa blocks",template[i]))
				{
					for (j in 1:length(clusters))
						{
							if ((!is.null(dim(clusters[[j]])))&&(dim(clusters[[j]])[1] >= 3)&(sum(!is.na(clusters[[j]][,"longitude"])) >= 3))
								{
									cat(paste0("\t<taxa id=\"taxa_",j,"\">","\n"))
									for (k in 1:dim(clusters[[j]])[1])
										{
											if (!is.na(clusters[[j]][k,"longitude"]))
												{
													collectionDate = round(as.numeric(clusters[[j]][k,"collection_date"]),6)
													latitude = round(as.numeric(clusters[[j]][k,"latitude"]),6)
													longitude = round(as.numeric(clusters[[j]][k,"longitude"]),6)												
													cat(paste0("\t\t<taxon id=\"",row.names(clusters[[j]])[k],"\">","\n"))
													cat(paste0("\t\t\t<date value=\"",collectionDate,"\" direction=\"forwards\" units=\"years\"/>","\n"))
													cat("\t\t\t<attr name=\"latitude\">\n")
													cat(paste0("\t\t\t\t",latitude,"\n"))
													cat("\t\t\t</attr>\n")
													cat("\t\t\t<attr name=\"longitude\">\n")
													cat(paste0("\t\t\t\t",longitude,"\n"))
													cat("\t\t\t</attr>\n")
													cat("\t\t\t<attr name=\"location\">\n")
													cat(paste0("\t\t\t\t",latitude," ",longitude,"\n"))
													cat("\t\t\t</attr>\n")
													cat("\t\t</taxon>\n")
												}
										}
									cat("\t</taxa>","\n")
								}
						}
				}
			if (grepl("Insert alignment blocks",template[i]))
				{
					for (j in 1:length(clusters))
						{
							if ((!is.null(dim(clusters[[j]])))&&(dim(clusters[[j]])[1] >= 3)&(sum(!is.na(clusters[[j]][,"longitude"])) >= 3))
								{
									cat(paste0("\t<alignment id=\"alignment_",j,"\" dataType=\"nucleotide\">","\n"))
									for (k in 1:dim(clusters[[j]])[1])
										{
											if (!is.na(clusters[[j]][k,"longitude"]))
												{
													cat("\t\t<sequence>\n")
													cat(paste0("\t\t\t<taxon idref=\"",row.names(clusters[[j]])[k],"\"/>","\n"))
													cat(paste0("\t\t\t",sequences[k],"\n"))
													cat("\t\t</sequence>\n")
												}
										}
									cat("\t</alignment>","\n")
								}
						}
				}
		}
	sink(NULL)
}
