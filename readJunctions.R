## this function takes a vector of splice junction positions, a reference BED file(reference genome
## in ncbi or uscs format), and a target chromosome.
## it returns a dataframe table with columns of accession number (geneid), position of splices, 
## and # of each splice position

library(GenomicRanges)
library(rtracklayer)
library(bedr)
library(Biostrings)
library(Rsamtools)
readJunctions= function(spliceJunction, ref_file, ref_chromosome){

y =data.table::fread(ref_file, data.table=FALSE, header=FALSE,select =c(1,2,3,4))
  
for (pos in spliceJunctions){
  ##subsets y for genes that are splicejunct'd
  spliced.genes = subset(y,V2<=pos & pos <=V3 & V1 == ref_chromosome)
  ##selects genes as to always be the first in the list, purges repeats,
  unique.genes = c(spliced.genes[1,1], spliced.genes[1,4])
  ##turns the names into a vector. 
  names= as.vector(y[4])
  ##adds counts column to track gene counts and pos.
  names['pos'] = vector('numeric')
  names['counts'] = vector('numeric')
  
  ##establishes the row in which the gene is found in names
  name.row=which(names[1:length(names[[1]]),1] %in% unique.genes[2])
  ##adds to names a count.
  if (is.na(names[name.row[1],2])){
    names[name.row[1],2]= pos
    names[name.row[1],3]=1
  }
  ##adds 1 to count if junction is counted already.
  if (!(is.na(names[name.row[1],2])) & pos %in% names[name.row[1],2])
{
  i = match(pos, names[name.row[1],2])
  c=names[name.row[1],3]
  c[i]+1
}
## if it's a new position, adds to vector in position slot.
if(!is.na(names[name.row[1],2]) & !(pos %in% names[name.row[1],2])){
  names[name.row[1],2]= c(names[name.row[1],2], pos)
  names[name.row[1],3]= c(names[name.row[1],3], 1)
  
}
}
} 
class(readJunctions)

















