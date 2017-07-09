## this function takes a BED file and extracts the chromosome, start, and stop coordinates.

readBeds=function(read_file){
    x= data.table::fread(read_file, data.table=FALSE, header=FALSE,select =c(1,2,3))
}
class (readBeds)