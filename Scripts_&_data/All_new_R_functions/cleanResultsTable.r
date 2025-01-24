cleanResultsTable = function(tab1) {
	
	addZeroDigits = function(tt, zeroDigits)
		{
			if (grepl("\\.",tt))
				{
					tt1 = unlist(strsplit(tt,"\\."))[1]
					tt2 = unlist(strsplit(tt,"\\."))[2]
				}	else	{
					tt1 = tt; tt2 = ""
				}
			if (nchar(tt2) < zeroDigits)
				{
					if ((zeroDigits-nchar(tt2)) == 3) tt2 = paste0(tt2,"000")
					if ((zeroDigits-nchar(tt2)) == 2) tt2 = paste0(tt2,"00")
					if ((zeroDigits-nchar(tt2)) == 1) tt2 = paste0(tt2,"0")
				}
			return(paste0(tt1,".",tt2))
		}
	tab2 = tab1
	for (i in 1:dim(tab2)[2])
		{
			col = colnames(tab2)[i]
			if (grepl("_BF",col))
				{
					for (j in 1:dim(tab2)[1])
						{
							if (as.character(tab1[j,i]) == "Inf")
								{
									tab2[j,i] = ">99"
								}	else	{
									if (tab1[j,i] != "-") tab2[j,i] = addZeroDigits(tt=as.character(tab1[j,i]), zeroDigits=1)
								}
						}
				}
			if (grepl("_pD",col)|grepl("_pQ",col))
				{
					for (j in 1:dim(tab2)[1])
						{
							if (tab1[j,i] != "-") tab2[j,i] = addZeroDigits(tt=as.character(tab1[j,i]), zeroDigits=2)
						}
				}
			if (grepl("_rP2",col)|grepl("_beta",col)|grepl("_R2",col)|grepl("_Q",col)|grepl("_D",col))
				{
					for (j in 1:dim(tab2)[1])
						{
							if (tab1[j,i] != "-")
								{
									if (!grepl("\\[",tab2[j,i]))
										{
											tab2[j,i] = addZeroDigits(as.character(tab1[j,i]),3)
										}
									if (grepl("\\[",tab2[j,i]))
										{
											t1 = unlist(strsplit(tab1[j,i]," \\["))[1]
											t = unlist(strsplit(tab1[j,i]," \\["))[2]
											t2 = unlist(strsplit(t,", "))[1]
											t3 = gsub("\\]","",unlist(strsplit(t,", "))[2])
											t1 = addZeroDigits(tt=t1, zeroDigits=3)
											t2 = addZeroDigits(tt=t2, zeroDigits=3)
											t3 = addZeroDigits(tt=t3, zeroDigits=3)
											tab2[j,i] = paste0(t1," [",t2,", ",t3,"]")
										}
								}
						}
				}
		}
	return(tab2)
}