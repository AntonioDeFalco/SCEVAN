sortData <- function(dataset, output_file_name = "")
{

    if (!file.exists(dataset)) {
	message("ERROR: Specified Input File does not Exist")
	return(FALSE)
    }
    if( output_file_name == "" || 
	    substr(output_file_name, 
		nchar(output_file_name), nchar(output_file_name)) == "/" ){
        message("ERROR: You must provide the name of the file",
                " in which sorted data are saved")
	return(FALSE)
    }
    res <- .C("sort_data", as.character(dataset),
              as.character(output_file_name))
}
