`species.name` <-
function(str)
{
    string<-gsub("[ ]","",str)
    string <- unlist(strsplit(string, NULL))      
    leftpar<-which(string=="(" | string==",")
    rightpar<-which(string==")")
    colon<-which(string==":")
    comma<-which(string==",")
    nspecies<-length(comma)+1
    name<-as.character(1:nspecies)
    iname<-1
    for(i in 1:length(colon)){
	x1<-leftpar[sum(leftpar<colon[i])]
 	x2<-sum(rightpar<colon[i])
	if(x2 > 0)
		x2<-rightpar[x2]
 	if(x1>x2){
		name[iname]<-paste(string[(x1+1):(colon[i]-1)],sep="",collapse="")
		iname<-iname+1
	}
    }
    return(name)
    
}

