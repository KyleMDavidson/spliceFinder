##function takes a chromosome, startrange and endrange, and a BAM file.
##function returns a numeric vector with the positions of all splice junctions found, or if there is no
##splicing in the transcripts (no cigar strings along a sequence that have gaps, so the definitiveness of this
##is dependent upon high coverage of the given area, it returns a numeric vector with no entries.


spliceFinder = function(chromosome,startRange,endRange,filepath){

##each numeral is 64 bits of memory, or 8 bytes. Max local memory:3gb
##creates a parameter for read ranges
##startRange= 13633501
##endRange = 13633680
##filepath = 'dataset6-8.bam'
#which <- RangesList(chromosome=IRanges(startRange,endRange))
which = GRanges(seqnames=chromosome,IRanges(startRange,endRange))
 ##specifies what we would like to extract from BAM
what <- c("pos", "seq")
##gives which ranges and what of them we want from BAM into param
param <- ScanBamParam(what=what, which=which)
ga = GenomicAlignments::readGAlignments(bamFile, filepath, FALSE, param,TRUE)
##check for splice candidacy
cigars=cigar(ga)
##builds vector from the cigar strands that are splice candidates
spliced= cigars[cigars!='40M']
##returns the indices that definte splice candidates in ga
splice.indices = match(spliced, cigars)
##uses candidate indices to subset GA to only reads with breaks
splice.candidates = ga[splice.indices]
length(splice.candidates)


##searches out the number of bases before a splice and adds them to a corresponding start
##sets up start positions
splice.starts = c(start(splice.candidates))
splice.additions = vector('numeric',length = length(splice.candidates))
##manages a correspondingly indexed vector in which to place the splice positions
splice.junctions = vector('numeric',length(splice.candidates))
cig = cigar(splice.candidates)


## inorder substr arrays. Gets indexes that identify positional mods that need to be handled in one dig or 2
oneval.splice.mods = substr(cig,2,2)
oneval.indices = which( oneval.splice.mods %in% 'M')
twoval.indices = which((!(oneval.splice.mods %in% 'M')))
cig[twoval.indices]
cig[oneval.indices]
##sets up vectors to hold the numerals
splice.additions[oneval.indices] =  (substr(cig[oneval.indices], 1,1))
splice.additions[twoval.indices] =  (substr(cig[twoval.indices], 1,2))


##determines splice junctions and fills in the splicejunction vector.
splice.addition = vector('numeric',length(splice.candidates))
splice.addition = as.numeric(splice.additions)

splice.junctions[oneval.indices] = splice.starts[oneval.indices]+splice.addition[oneval.indices]
splice.junctions[twoval.indices] = splice.starts[twoval.indices]+splice.addition[twoval.indices]

}
class (spliceFinder)





