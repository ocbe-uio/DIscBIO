#' @title Plotting the network.
#' @description This function uses STRING-api to plot the network. 
#' @export
#' @param data A gene list.
#' @param FileName A string vector showing the name to be used to save the resulted network.
#' @param species The taxonomy name/id. Default is "9606" for Homo sapiens.
#' @importFrom httr GET status_code
#' @importFrom utils download.file
#' @importFrom png readPNG
#' @importFrom graphics plot rasterImage
Networking<-function(data,FileName,species="9606"){
	if (length(data)>600){
		print("Your gene list is too big")
			
	}else{
		string_api_url <- "https://string-db.org/api/"
		output_format <- "image"
		method <- "network"
		your_identifiers <- ""
		optional_parameters <- ""

		# Construct API request
		genes <- data
		repos <- GET(url = paste0(string_api_url,output_format,'/',method,'?identifiers=',
                          paste(as.character(data), collapse = "%0d"),"&", "species=",species))
                cat("Examine response components =",status_code(repos),"\t","200 means successful","\n")
		y = repos$request$ url
		download.file(y,paste0("network",FileName,".png"), mode = 'wb')
		Network<- readPNG(paste0("network",FileName,".png"), native = TRUE)
		plot(0:1,0:1,type="n",ann=FALSE,axes=FALSE)
		rasterImage(Network,0,0,1,1)	
	}
}