#' HTGM
#'
#' @import minimalistGODB
#' @import GoMiner
#' @import grDevices
#' @import stats
#' @importFrom gplots heatmap.2
#' @import vprint
#'
#' @description driver to invoke GoMiner for multiple studies, and integrate the results
#'  in a categories versus study hyperlinked heatmap
#'
#' @param title character string descriptive title
#' @param dir character string full pathname to the directory acting as result repository
#' @param sampleLists list of character vector of user-supplied genes of interest
#' @param GOGOA3 return value of subsetGOGOA()
#' @param ONT c("molecular_function","cellular_component","biological_process") 
#' @param enrichThresh  numerical acceptance threshold for enrichment passed to GoMiner
#' @param countThresh numerical acceptance threshold for gene count passed to GoMiner
#' @param fdrThresh numerical acceptance threshold for fdr passed to GoMiner
#' @param nrand integer number of randomizations passed to GoMiner
#' @param mn integer param passed to trimGOGOA3, min size threshold for a category
#' @param mx integer param passed to trimGOGOA3, max size threshold for a category
#' @param opt integer 0:1 parameter used to select randomization method
#' @param verbose integer parameter passed to vprint()
#'
#' @examples
#' \dontrun{
#' # GOGOA3.RData is too large to include in the R package
#' # so I need to load it from a file that is not in the package.
#' # Since this is in a file in my own file system, I could not
#' # include this as a regular example in the package.
#' # you can generate it using the package 'minimalistGODB'
#' # or you can retrieve it from https://github.com/barryzee/GO/tree/main/databases
#' load("/Users/barryzeeberg/personal/GODB_RDATA/goa_human/GOGOA3_goa_human.RData")
#' 
#' # load("data/Housekeeping_Genes.RData")
#' sampleList<-unique(as.matrix(Housekeeping_Genes[,"Gene.name"]))
#' n<-nrow(sampleList)
#' sampleLists<-list()
#' # test the effect of random sampling of the entire gene set
#' # this can give an idea of the quality of the GoMiner results
#' # when the complete gene set is yet to be determined
#' sampleLists[["1"]]<-sampleList[sample(n,n/2)]
#' sampleLists[["2"]]<-sampleList[sample(n,n/2)]
#' sampleLists[["3"]]<-sampleList[sample(n,n/2)]
#' sampleLists[["4"]]<-sampleList[sample(n,n/2)]
#' sampleLists[["5"]]<-sampleList[sample(n,n/2)]
#' sampleLists[["ALL"]]<-sampleList
#' m<-HTGM(title=NULL,dir=tempdir(),sampleLists,GOGOA3,ONT="biological_process",
#'  enrichThresh=2,countThresh=5,fdrThresh=0.10,nrand=100)
#' }
#' 
#' @return returns the matrix of significant categories versus study
#' 
#' @export
HTGM<-
  function(title=NULL,dir=tempdir(),sampleLists,GOGOA3,ONT,enrichThresh=2,countThresh=5,
           fdrThresh=0.10,nrand=100,mn=2,mx=200,opt=0,verbose=1) {
    stamp<-gsub(":","_",format(Sys.time(), "%a_%b_%d_%Y_%X"))
    if(is.null(title))
      title<-"HTGMresults"
    subd<-sprintf("%s/%s_%s",dir,title,stamp)
    dir.create(subd)
    
    l<-list()
    if(verbose)
      class<-"none"
    else
      class<-"message"
    for(id in names(sampleLists))
      #l[[id]]<-suppressMessages(GoMiner(title=id,subd,sampleLists[[id]],GOGOA3,
       # ONT,enrichThresh,countThresh,fdrThresh,nrand),classes=class)
      l[[id]]<-suppressMessages(GoMiner(title=id,subd,sampleLists[[id]],GOGOA3,
        ONT,enrichThresh,countThresh,fdrThresh,nrand,mn,mx,opt,verbose),classes=class)
    
    #x_l<-l
    #save(x_l,file="data/x_l.RData")
    # set up the matrix of significant categories versus gene sample list title
    m<-htgmM(l,fdrThresh)
    
    # set up the genes mapping to each significant category per gene sample list title
    hyperGenes(l,subd)
    
    file<-sprintf("%s/htgm.svg",subd)
    #svgWidth<-(0.375+0.025517241*ncol(m)) * 1.059602649
    #svgHeight<-(8.0 + 0.5*nrow(m)) * 0.526
    
    svgWidth<-(0.375+0.025517241*ncol(m)) * 1.059602649 * 4
    svgHeight<-(8.0 + 0.5*nrow(m)) * 0.526 * 4

    svg(filename=file,width=svgWidth,height=svgHeight*1.5)
    # trick - use row for 'key' to get more space for long category names, but suppress key
    
    hm<-heatmap.2(m,col = heat.colors(n=100,rev=FALSE),trace="none",lmat=rbind( c(0, 3),
        c(2,1), c(0,4) ),lhei=c(1,4,15),lwid=c(1,50),key=FALSE,margins = c(1, 10))
    
    dev.off()
    #file.copy(file,"inst/extdata/x_htgm.svg")
    #x_m<-m
    #save(x_m,file="data/x_m.RData")
    
    ff<-hyperlinks(file,rownames(m[hm$rowInd,]),colnames(m[,hm$colInd]))
    
    return(m)
  }

#' hyperlinks
#' 
#' @description driver to add gene list hyperlinks to the HTGM heatmap
#' 
#' @param s character path name of the file containing the HTGM svg
#' @param rownames character vector of row names
#' @param colnames character vector of column names
#' 
#' @examples
#' #load("data/x_rn.RData")
#' #load("data/x_cn.RData")
#' #load("data/x_svg.RData")
#' s<-system.file("extdata","x_htgm.svg",package="HTGM")
#' # need to avoid writing to "extdata"
#' dir<-tempdir()
#' file.copy(from=s, to=dir)
#' hyperlinkedFileName<-hyperlinks(sprintf("%s/%s",dir,"x_htgm.svg"),x_rn,x_cn)
#' print("hyperlinkedFileName")
#' print(hyperlinkedFileName)
#' 
#' @return returns the path name of the file containing the hyperlinked HTGM svg
#' 
#' @export
hyperlinks<-
  function(s,rownames,colnames) {
    
    SVG<-scan(s,what="character",sep="\n")
    
    #x_rn<-rownames
    #save(x_rn,file="data/x_rn.RData")
    #x_cn<-colnames
    #save(x_cn,file="data/x_cn.RData")
    #x_svg<-SVG
    #save(x_svg,file="data/x_svg.RData")
    
    nr<-length(rownames)
    nc<-length(colnames)
    
    n<-0
    row<-0
    col<-1
    for(r in 1:length(SVG)) {
      #gr1<-grep("<path style=",SVG[r])
      gr1<-grep("<path fill-rule=",SVG[r])
      gr2<-grep("nonzero",SVG[r])
      if(length(gr1)>0 && length(gr2)>0) {
        n<-n+1
        if(row==nr) {
          row<-0
          col<-col+1
        }
        row<-row+1
        #x_svgr<-SVG[r]
        #save(x_svgr,file="data/x_svgr.RData")
        #x_rnr<-rownames[row]
        #save(x_rnr,file="data/x_rnr.RData")
        #x_cnc<-colnames[col]
        #save(x_cnc,file="data/x_cnc.RData")
        SVG[r]<-pasteHyperlinks(SVG[r],rownames[row],colnames[col])
      }
      if(n==nr*nc) break
    }
   
    fname<-sprintf("%s/%s",dirname(s),"hyperlink.svg")
    writeLines(SVG,fname)
    
    return (fname)
  }

#' pasteHyperlinks
#' 
#' @description add gene list hyperlinks to the HTGM heatmap
#' 
#' @param str character a line from the svg that is to have a hyperlink inserted
#' @param c1 character list of row names
#' @param c2 character list of column names
#' 
#' @examples
#' #load("data/x_svgr.RData")
#' #load("data/x_rnr.RData")
#' #load("data/x_cnc.RData")
#' hl<-pasteHyperlinks(x_svgr,x_rnr,x_cnc)
#'
#' @return returns a line of code to insert into svg
#' 
#' @export
pasteHyperlinks<-
  function(str,c1,c2) {
    # allow this to work for both HTGM and HTGM2D
    if(substr(c1,1,3)=="GO_")
      c1<-substr(c1,1,10)
    pre <- sprintf("%s%s%s__%s%s%s","<a xlink:href= '","./hyperGenes/",c1,substr(c2,1,10),
              ".txt","' target='new'>")
    post <- "</a>"
    
    return (paste(pre,str,post))
  }

#' hyperGenes
#' 
#' @description populate subdirectory of hyperlinked gene lists
#' 
#' @param l return value of GoMiner()
#' @param dir character string containing path name of results directory
#' 
#' @examples
#' # x_l<-load("data/x_l.RData")
#' dir<-tempdir()
#' print(dir)
#' hyperGenes(x_l,dir)
#'  
#' @return returns no value but has side effect of populating subdirectory of
#'  hyperlinked gene lists
#' 
#' @export
hyperGenes<-
  function(l,dir) {
    # pathname to subdirectory of hyperlinked gene lists
    
    # debug for ngene = 1 cat
    #save(l,file="~/gce.RData")

    subdir<-sprintf("%s/%s",dir,"hyperGenes")
    if(!dir.exists(subdir))
      dir.create(subdir)
    for(id in names(l)) {
      thresh<-l[[id]][["thresh"]]
      gce<-l[[id]][["gce"]]
      for(i in 1:nrow(thresh)) {
        cat<-thresh[i,"Row.names"]
        w<-which(gce[,"GO_NAME"]==cat)
        genes<-gce[w,"HGNC"]
        ##label<-sprintf("%s__%s.txt",id,cat)
        label<-sprintf("%s__%s.txt",id,substr(cat,1,10))
        writeLines(sort(unique(genes)),sprintf("%s/%s",subdir,label))
      }
    }
  }

#' htgmM
#' 
#' @description generate FDR matrix of id versus cat
#' 
#' @param l list of return values of GoMiner()
#' @param fdrThresh numerical acceptance threshold for fdr
#' 
#' @examples
#' # load("data/x_l.RData")
#' m<-htgmM(x_l,.1)
#' 
#' @return returns numeric matrix m containing FDR values
#'
#' @export     
htgmM<-
  function(l,fdrThresh) {
    # list of all uniqued categories in thresh's
    v<-vector("character",0)
    for(id in names(l)) {
      thresh<-l[[id]][["thresh"]]
      v<-unique(c(v,thresh[,"Row.names"]))
    }
    
    # m is FDR matrix of id versus cat
    m<-matrix(fdrThresh,nrow=length(names(l)),ncol=length(v))
    rownames(m)<-names(l)
    colnames(m)<-v
    
    for(id in names(l)) {
      thresh<-l[[id]][["thresh"]]
      for(i in 1:nrow(thresh))
        m[id,thresh[i,"Row.names"]]<-thresh[i,"FDR"]
    }
 
    return(m)
  }
