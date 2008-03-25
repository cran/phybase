`name2node` <-
function(treestr,name="")
{
    str<-treestr    
    if(length(name)<2)
	speciesname<-sort(species.name(str))
    else
	speciesname<-name
    
    new<-paste(1:length(speciesname),":",sep="")
    old<-paste(speciesname,":",sep="")
    strlength<-nchar(old)
    inputorder<-order(rank(strlength,ties.method="random"),decreasing=TRUE)
    for(i in 1:length(speciesname))	
	str<-gsub(old[inputorder[i]],new[inputorder[i]],str)

    str

}