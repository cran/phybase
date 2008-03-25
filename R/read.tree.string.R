`read.tree.string` <-
function (file = "",format="nexus") 
{
    
   X <- scan(file = file, what = "", sep = "\n", quiet = TRUE)

    i1 <- grep("BEGIN TREES;", X, ignore.case = TRUE)
    i2 <- grep("TRANSLATE", X, ignore.case = TRUE)
    translation <- FALSE
    if (length(i2) == 1 && length(i1) == 1){
	if(i2 > i1) translation<-TRUE
    }
    if (translation) {
       	semico<-grep(";",X) 
        end <- semico[semico > i2][1]
        x <- X[(i2 + 1):end]
        x <- unlist(strsplit(x, "[,; \t]"))
        x <- x[nzchar(x)]
        TRANS <- matrix(x, ncol = 2, byrow = TRUE)
        TRANS[, 2] <- gsub("['\"]", "", TRANS[, 2])
        nspecies <- dim(TRANS)[1]
	speciesname<-TRANS[,2]
    }

    X<- gsub("[]]", "", X)
    X<- gsub("[[]", "", X)

    X<-paste(X,collapse="")
    X<-unlist(strsplit(X,split=";"))
    X<-paste(X,";",sep="")

    if(length(grep("nexus",format,ignore.case=TRUE))>0)
    	X<-unlist(strsplit(X,split="="))
    
    
    tree<-rep("", length(X))
    itree<-1
    for(i in 1:length(X)){
  	string <- unlist(strsplit(X[i], NULL))
	leftpar<-which(string=="(")
    	rightpar<-which(string==")")
     	semicolon<-which(string==":")
    	if(length(leftpar)>1 & length(leftpar)>1 & length(semicolon)>1 ){
		left<-which(string=="(")[1]
		right<-which(string==";")[1]
		tree[itree]<-paste(string[left:right],sep="",collapse="")
		itree<-itree+1
	}
    }

    tree<-tree[tree != ""]

    if(!translation)
    {   

	speciesname<-species.name(tree[1])

    }

    string <- unlist(strsplit(tree[1], NULL))
    leftpar<-which(string=="(")  
    rightpar<-which(string==")") 
    comma<-which(string==",")
    if(length(leftpar) != length(leftpar))
     	stop("The number of left parenthesis is NOT equal to the number of right  parenthesis")


    if(length(leftpar) == length(comma))
	rooted<-TRUE   else if(length(leftpar) == (length(comma)-1))
	rooted<-FALSE  else
        stop("Wrong number of comma in the tree string!")

    z <- list(tree="", names ="", root=TRUE)

    z$tree<-tree
    z$names<-speciesname
    z$root<-rooted
    z
}

