`write.subtree` <-
function(inode,nodes,nspecies,inodeindex)
{
    if(inode<=nspecies){
	if(nodes[inode,4] >= 0){
       		str<-paste(inode,":",sep="")
       		str<-paste(str,nodes[inode,4],sep="")
	}
    	if(nodes[inode,4] < 0)
		str<-paste(inode)
	if(nodes[inode,5] > 0){
		str<-paste(str,"#",sep="")
		str<-paste(str,nodes[inode,5],sep="")
	}
    }
    if(inode>nspecies){
        son1<-nodes[inode,2]
        son2<-nodes[inode,3]
	str1<-write.subtree(son1,nodes,nspecies,inodeindex)
	str2<-write.subtree(son2,nodes,nspecies,inodeindex)
	str<-paste("(",str1,sep="")
 	str<-paste(str,",",sep="")
	str<-paste(str,str2,sep="")
	if(nodes[inode,1] == -8){
		son3<-as.integer(nodes[inode,4])
		str3<-write.subtree(son3,nodes,nspecies,inodeindex)
		str<-paste(str,",",sep="")
		str<-paste(str,str3,sep="")
	}
	str<-paste(str,")",sep="")
	if(nodes[inode,4] >= 0 & nodes[inode,1]>0 & inode != inodeindex){
		str<-paste(str,":",sep="")
		str<-paste(str,round(nodes[inode,4],6),sep="")
	}
	if(nodes[inode,5] >= 0){
		str<-paste(str,"#",sep="")
		str<-paste(str,round(nodes[inode,5],6),sep="")
	}
    }
    if(inode == inodeindex)
	str<-paste(str,";",sep="")
    return(str)
}
