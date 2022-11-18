#' Generates lagrangePHYLIP file for using in BioGeoBEARS
#' 
#' It prepares and save the file to be used as geographical areas in `BioGeoBEARS`
#' scripts. It was based on `save_tipranges_to_LagrangePHYLIP` function of
#' `BioGEOBEARS` package.
#' 
#' @param tip_range data frame. Tip names in rows, geographical areas in columns. 
#'   Please provide the tip names as `row.names`. See example. 
#' @param filename filename to store the data. By default it is 
#'   `"lagrange_area_data_file.data"`. 
#' @param areanames the names of the areas. By defaut the area names are generated
#'  as sequencial upper case letters. 
#'  
#' @importFrom utils write.table
#'  
#' @examples
#' 
#' \dontrun{
#' set.seed(082622)
#' geo.areas <- data.frame(
#'   a = rbinom(5, 1, 0.5),
#'   b = rbinom(5, 1, 0.5),
#'   c = rbinom(5, 1, 0.5)
#' )
#' 
#' # set the row names
#' row.names(geo.areas) <- paste0("tip_", 1:nrow(geo.areas))
#' 
#' # it saves the data at 'filename' path.
#' tipranges_to_BioGeoBEARS(
#'   geo.areas, 
#'   filename = "lagrange_area_data_file.data", 
#'   areanames = NULL)
#' )
#' 
#' # file content must be like this:
#' # 5	3	(A B C)
#' # tip_1	111
#' # tip_2	001
#' # tip_3	010
#' # tip_4	101
#' # tip_5	010
#'  }
#'  
#' @export

tipranges_to_BioGeoBEARS <- function (
    tip_range, 
    filename = "lagrange_area_data_file.data", 
    areanames = NULL) 
{
  tipranges_df = tip_range
  
  
  if ((is.null(areanames)) || (length(areanames) != ncol(tipranges_df))) {
    new_areanames = LETTERS[1:ncol(tipranges_df)]
    areanames_txt = paste(new_areanames, collapse = " ")
    warning("\nNote: assigning '", areanames_txt, "' as area names.\n", 
        sep = "")
  }
  else {
    areanames_txt = paste(areanames, collapse = " ")
  }
  
  ranges_strings = apply(X = tipranges_df, MARGIN = 1, FUN = paste, 
                         collapse = "")
  taxon_names = row.names(tipranges_df)
  ranges_table = cbind(taxon_names, ranges_strings)
  ranges_table_txt = apply(X = ranges_table, MARGIN = 1, FUN = paste, 
                           collapse = "\t")
  header_string = paste(nrow(tipranges_df), "\t", ncol(tipranges_df), 
                        "\t(", areanames_txt, ")", sep = "")
  write.table(x = header_string, file = filename, append = FALSE, 
              quote = FALSE, sep = "\n", row.names = FALSE, col.names = FALSE)
  write.table(x = ranges_table_txt, file = filename, append = TRUE, 
              quote = FALSE, sep = "\n", row.names = FALSE, col.names = FALSE)
  return(filename)
}
