generateXMLfiles1 = function(templateFile, newXMLFile, prefix, treeFile, metadata, jitter) {

	template = scan(templateFile, what="", sep="\n", quiet=T, blank.lines.skip=F)
	sink(file=newXMLFile); clusters = list(); clusters[[1]] = metadata
	for (i in 1:length(template))
		{
			if (grepl("TEMPLATE",template[i]))
				{
					template[i] = gsub("TEMPLATE", prefix, template[i])
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
													collectionDate = round(clusters[[j]][k,"collection_date"],6)
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
													cat("\t\t\tNNNN\n")
													cat("\t\t</sequence>\n")
												}
										}
									cat("\t</alignment>","\n")
								}
						}
				}
			if (grepl("Insert pattern blocks",template[i]))
				{
					for (j in 1:length(clusters))
						{
							if ((!is.null(dim(clusters[[j]])))&&(dim(clusters[[j]])[1] >= 3)&(sum(!is.na(clusters[[j]][,"longitude"])) >= 3))
								{
									cat(paste0("\t<patterns id=\"patterns_",j,"\" from=\"1\" strip=\"false\">","\n"))
									cat(paste0("\t\t<alignment idref=\"alignment_",j,"\"/>","\n"))
									cat("\t</patterns>","\n")
								}
						}
				}
			if (grepl("<newick id=\"startingTree\">",template[i]))
				{
					for (j in 1:length(clusters))
						{
							if ((!is.null(dim(clusters[[j]])))&&(dim(clusters[[j]])[1] >= 3)&(sum(!is.na(clusters[[j]][,"longitude"])) >= 3))
								{
									tre = scan(treeFile, what="", sep="\n", quiet=T)
									cat(paste0("\t\t",tre,"\n"))
								}
						}
				}
			if (grepl("Insert tree model blocks",template[i]))
				{
					for (j in 1:length(clusters))
						{
							if ((!is.null(dim(clusters[[j]])))&&(dim(clusters[[j]])[1] >= 3)&(sum(!is.na(clusters[[j]][,"longitude"])) >= 3))
								{
									cat(paste0("\t<treeModel id=\"treeModel_",j,"\">","\n"))
									cat(paste0("\t\t<coalescentTree idref=\"startingTree_",j,"\"/>","\n"))
									cat("\t\t<rootHeight>","\n")
									cat(paste0("\t\t\t<parameter id=\"treeModel.rootHeight_",j,"\"/>","\n"))
									cat("\t\t</rootHeight>","\n")
									cat("\t\t<nodeHeights internalNodes=\"true\">","\n")
									cat(paste0("\t\t\t<parameter id=\"treeModel.internalNodeHeights_",j,"\"/>","\n"))
									cat("\t\t</nodeHeights>","\n")
									cat("\t\t<nodeHeights internalNodes=\"true\" rootNode=\"true\">","\n")
									cat(paste0("\t\t\t<parameter id=\"treeModel.allInternalNodeHeights_",j,"\"/>","\n"))
									cat("\t\t</nodeHeights>","\n")
									cat("\t</treeModel>","\n")
								}
						}
				}
			if (grepl("Insert arbitraryBranchRates blocks",template[i]))
				{
					for (j in 1:length(clusters))
						{
							if ((!is.null(dim(clusters[[j]])))&&(dim(clusters[[j]])[1] >= 3)&(sum(!is.na(clusters[[j]][,"longitude"])) >= 3))
								{
									cat(paste0("\t<arbitraryBranchRates id=\"location.diffusion.branchRates",j,"\">","\n"))
									cat(paste0("\t\t<treeModel idref=\"treeModel_",j,"\"/>","\n"))
									cat("\t\t<rates>","\n")
									cat(paste0("\t\t\t<parameter id=\"location.diffusion.rates",j,"\" lower=\"0.0\"/>","\n"))
									cat("\t\t</rates>","\n")
									cat("\t</arbitraryBranchRates>","\n")
								}
						}
				}
			if (grepl("Insert distributionLikelihood blocks 1",template[i]))
				{
					for (j in 1:length(clusters))
						{
							if ((!is.null(dim(clusters[[j]])))&&(dim(clusters[[j]])[1] >= 3)&(sum(!is.na(clusters[[j]][,"longitude"])) >= 3))
								{
									cat(paste0("\t<distributionLikelihood id=\"location.diffusion.prior",j,"\">","\n"))
									cat("\t\t<data>","\n")
									cat(paste0("\t\t\t<parameter idref=\"location.diffusion.rates",j,"\"/>","\n"))
									cat("\t\t</data>","\n")
									cat("\t\t<distribution>","\n")
									cat(paste0("\t\t\t<onePGammaDistributionModel>","\n"))
									cat("\t\t\t\t<shape>","\n")
									cat("\t\t\t\t\t<parameter value=\"0.5\"/>","\n")
									cat("\t\t\t\t</shape>","\n")
									cat("\t\t\t</onePGammaDistributionModel>","\n")
									cat("\t\t</distribution>","\n")
									cat("\t</distributionLikelihood>","\n")
								}
						}
				}
			if (grepl("Insert location.traitLikelihood blocks",template[i]))
				{
					for (j in 1:length(clusters))
						{
							if ((!is.null(dim(clusters[[j]])))&&(dim(clusters[[j]])[1] >= 3)&(sum(!is.na(clusters[[j]][,"longitude"])) >= 3))
								{
									cat(paste0("\t<multivariateTraitLikelihood id=\"location.traitLikelihood",j,"\" traitName=\"location\" useTreeLength=\"true\" scaleByTime=\"true\" reportAsMultivariate=\"true\" reciprocalRates=\"true\" integrateInternalTraits=\"true\">","\n"))
									cat("\t\t<multivariateDiffusionModel idref=\"location.diffusionModel\"/>","\n")
									cat(paste0("\t\t<treeModel idref=\"treeModel_",j,"\"/>","\n"))
									cat("\t\t<traitParameter>","\n")
									cat(paste0("\t\t\t<parameter id=\"leaf.location",j,"\"/>","\n"))
									cat("\t\t</traitParameter>","\n")
									cat(paste0("\t\t<jitter window=\"",jitter," ",jitter,"\" duplicatesOnly=\"true\">","\n"))
									cat(paste0("\t\t\t<parameter idref=\"leaf.location",j,"\"/>","\n"))
									cat("\t\t</jitter>","\n")
									cat("\t\t<conjugateRootPrior>","\n")
									cat("\t\t\t<meanParameter>","\n")
									cat("\t\t\t\t<parameter value=\"0.0 0.0\"/>","\n")
									cat("\t\t\t</meanParameter>","\n")
									cat("\t\t\t<priorSampleSize>","\n")
									cat("\t\t\t\t<parameter value=\"0.000001\"/>","\n")
									cat("\t\t\t</priorSampleSize>","\n")
									cat("\t\t</conjugateRootPrior>","\n")
									cat(paste0("\t\t<arbitraryBranchRates idref=\"location.diffusion.branchRates",j,"\"/>","\n"))
									cat("\t</multivariateTraitLikelihood>","\n")
								}
						}
				}
			if (grepl("Insert continuousDiffusionStatistic blocks 1",template[i]))
				{
					for (j in 1:length(clusters))
						{
							if ((!is.null(dim(clusters[[j]])))&&(dim(clusters[[j]])[1] >= 3)&(sum(!is.na(clusters[[j]][,"longitude"])) >= 3))
								{
									cat(paste0("\t<continuousDiffusionStatistic id=\"location.diffusionRate",j,"\" greatCircleDistance=\"true\">","\n"))
									cat(paste0("\t\t<multivariateTraitLikelihood idref=\"location.traitLikelihood",j,"\"/>","\n"))
									cat("\t</continuousDiffusionStatistic>","\n")
								}
						}
				}
			if (grepl("Insert scaleOperator blocks",template[i]))
				{
					for (j in 1:length(clusters))
						{
							if ((!is.null(dim(clusters[[j]])))&&(dim(clusters[[j]])[1] >= 3)&(sum(!is.na(clusters[[j]][,"longitude"])) >= 3))
								{
									cat(paste0("\t\t<scaleOperator scaleFactor=\"0.75\" weight=\"30\">","\n"))
									cat(paste0("\t\t\t<parameter idref=\"location.diffusion.rates",j,"\"/>","\n"))
									cat("\t\t</scaleOperator>","\n")
								}
						}
				}
			if (grepl("Insert precisionGibbsOperator blocks",template[i]))
				{
					for (j in 1:length(clusters))
						{
							if ((!is.null(dim(clusters[[j]])))&&(dim(clusters[[j]])[1] >= 3)&(sum(!is.na(clusters[[j]][,"longitude"])) >= 3))
								{
									cat(paste0("\t\t<precisionGibbsOperator weight=\"2\">","\n"))
									cat(paste0("\t\t\t<multivariateTraitLikelihood idref=\"location.traitLikelihood",j,"\"/>","\n"))
									cat("\t\t\t<multivariateWishartPrior idref=\"location.precisionPrior\"/>","\n")
									cat("\t\t</precisionGibbsOperator>","\n")
								}
						}
				}
			if (grepl("Insert distributionLikelihood blocks 2",template[i]))
				{
					for (j in 1:length(clusters))
						{
							if ((!is.null(dim(clusters[[j]])))&&(dim(clusters[[j]])[1] >= 3)&(sum(!is.na(clusters[[j]][,"longitude"])) >= 3))
								{
									cat(paste0("\t\t\t\t<distributionLikelihood idref=\"location.diffusion.prior",j,"\"/>","\n"))
								}
						}
				}
			if (grepl("Insert multivariateTraitLikelihood blocks 1",template[i]))
				{
					for (j in 1:length(clusters))
						{
							if ((!is.null(dim(clusters[[j]])))&&(dim(clusters[[j]])[1] >= 3)&(sum(!is.na(clusters[[j]][,"longitude"])) >= 3))
								{
									cat(paste0("\t\t\t\t<multivariateTraitLikelihood idref=\"location.traitLikelihood",j,"\"/>","\n"))
								}
						}
				}
			if (grepl("Insert continuousDiffusionStatistic blocks 2",template[i]))
				{
					for (j in 1:length(clusters))
						{
							if ((!is.null(dim(clusters[[j]])))&&(dim(clusters[[j]])[1] >= 3)&(sum(!is.na(clusters[[j]][,"longitude"])) >= 3))
								{
									cat(paste0("\t\t\t\t<continuousDiffusionStatistic idref=\"location.diffusionRate",j,"\"/>","\n"))
								}
						}
				}
			if (grepl("Insert multivariateTraitLikelihood blocks 2",template[i]))
				{
					for (j in 1:length(clusters))
						{
							if ((!is.null(dim(clusters[[j]])))&&(dim(clusters[[j]])[1] >= 3)&(sum(!is.na(clusters[[j]][,"longitude"])) >= 3))
								{
									cat(paste0("\t\t\t\t<multivariateTraitLikelihood idref=\"location.traitLikelihood",j,"\"/>","\n"))
								}
						}
				}
			if (grepl("<!-- Insert logTree blocks -->",template[i]))
				{
					for (j in 1:length(clusters))
						{
							if ((!is.null(dim(clusters[[j]])))&&(dim(clusters[[j]])[1] >= 3)&(sum(!is.na(clusters[[j]][,"longitude"])) >= 3))
								{
									cat(paste0("\t\t<logTree id=\"treeFileLog",j,"\" logEvery=\"100000\" nexusFormat=\"true\" fileName=\"",prefix,".trees\" sortTranslationTable=\"true\">","\n"))
									cat(paste0("\t\t\t<treeModel idref=\"treeModel_",j,"\"/>","\n"))
									cat("\t\t\t<joint idref=\"joint\"/>","\n")
									cat("\t\t\t<trait name=\"location\" tag=\"location\">","\n")
									cat(paste0("\t\t\t\t<multivariateTraitLikelihood idref=\"location.traitLikelihood",j,"\"/>","\n"))
									cat("\t\t\t</trait>","\n")
									cat("\t\t\t<multivariateDiffusionModel idref=\"location.diffusionModel\"/>","\n")
									cat("\t\t\t<trait name=\"rate\" tag=\"location.rate\">","\n")
									cat(paste0("\t\t\t\t<arbitraryBranchRates idref=\"location.diffusion.branchRates",j,"\"/>","\n"))
									cat("\t\t\t</trait>","\n")
									cat("\t\t</logTree>","\n")
								}
						}
				}
		}
	sink(NULL)
}
