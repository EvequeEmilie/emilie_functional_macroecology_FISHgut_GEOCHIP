
levels_possible <- c("Genbank.ID", "Genbank.id", "Gene.name","Gene.category", "Sub.category1","Sub.category2", "Sub.category3","Sub.category4","Gene","Gene_category","GeneCategory","subcategory1", "subcategory2","Subcategory1", "Subcategory2","Organism", "Lineage","Annotation","Genbank.GI","Gene.Name", "annotation","ProbeCategory","Probe.Category","X","Species","OTU_ID","OTU","Superkingdom","Kingdom","Domain","Phylum","Class","Order","Family","Genus","Species", "Sk","Kg","Ph","Cl","Or","Fa","Ge","Sp","superkingdom","kingdom","phylum","class","order","family","genus","species","TraitType")

L = function(data)
{
labs=c("Genbank.ID", "Genbank.id", "Gene.name","Gene.category", "Sub.category1","Sub.category2", "Sub.category3","Sub.category4","Gene","Gene_category","GeneCategory","subcategory1", "subcategory2","Subcategory1", "Subcategory2","Organism", "Lineage","Annotation","Genbank.GI","Gene.Name", "annotation","ProbeCategory","Probe.Category","X","Species","OTU_ID","OTU","Superkingdom","Kingdom","Domain","Phylum","Class","Order","Family","Genus","Species", "Sk","Kg","Ph","Cl","Or","Fa","Ge","Sp","superkingdom","kingdom","phylum","class","order","family","genus","species","TraitType")

	dataL = data[,colnames(data) %in% labs]		#data labels
return(dataL)
}

S = function(data)
{	
labs=c("Genbank.ID", "Genbank.id", "Gene.name","Gene.category", "Sub.category1","Sub.category2", "Sub.category3","Sub.category4","Gene","Gene_category","GeneCategory","subcategory1", "subcategory2","Subcategory1", "Subcategory2","Organism", "Lineage","Annotation","Genbank.GI","Gene.Name", "annotation","ProbeCategory","Probe.Category","X","Species","OTU_ID","OTU","Superkingdom","Kingdom","Domain","Phylum","Class","Order","Family","Genus","Species", "Sk","Kg","Ph","Cl","Or","Fa","Ge","Sp","superkingdom","kingdom","phylum","class","order","family","genus","species","TraitType")

	dataS = data[,!colnames(data) %in% labs]	#data Samples
return(dataS)
}



