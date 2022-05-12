## To Supress Note
utils::globalVariables(c(".", "..mycols","..tcga_pheno_columns","tempdir","heat.colors","topo.colors"))

#' Function to extract the clinical annotations from TCGA
#' @description This function download and extract the clinical
#' annotations from TCGA.
#' @author Mayank Tandon, Ashish Jain
#' @param cancerCodes A character vector of TCGA cancer codes
#' @param outputFolder path of the file containing the clinical
#' annotations from TCGA
#' @param plotdata Flag to plot the annotations
#' @export
#' @return A list containing the TCGA clinical annotations
#'
#' @importFrom TCGAbiolinks GDCquery_clinic
#'
#' @examples
#' library(MAFDash)
#' cancerCodes <- "ACC"
#' #clinicalData <- getTCGAClinicalAnnotation(cancerCodes = cancerCodes)
#'
getTCGAClinicalAnnotation <- function(cancerCodes="ACC",outputFolder=tempdir(), plotdata=NULL) {
  cancerCodes <- ensurer::ensure_that(cancerCodes,
                                      !is.null(.) && (class(.) == "character") && all(. %in% TCGAbiolinks::getGDCprojects()$tumor),
                                      err_desc = "Please enter the correct TCGA cancer code. See TCGAbiolinks::getGDCprojects().")
  plotdata <- ensurer::ensure_that(plotdata,
                                   is.null(.) || (class(.) == "data.frame"),
                                   err_desc = "Please enter the data to be plotted in a correct format.")

  tcgaClinicalAnnotationDataList<-list()
  for(cancerCode in cancerCodes)
  {
    tcga_clinical_file=file.path(outputFolder,paste0("TCGA_",cancerCode,".clinical.txt"))
    if (! file.exists(tcga_clinical_file)) {
      if (!dir.exists(dirname(tcga_clinical_file))) {dir.create(dirname(tcga_clinical_file), recursive = TRUE)}
      tcga_clinical <- TCGAbiolinks::GDCquery_clinic(project = paste0("TCGA-",cancerCode), type = "clinical")
      write.table(tcga_clinical, file=tcga_clinical_file, quote=TRUE, sep="\t", row.names = FALSE, col.names = TRUE)
    }
    tcga_clin_data <- read.table(tcga_clinical_file, sep="\t",header = TRUE,stringsAsFactors = FALSE)
    tcga_clin_data$Tumor_Sample_Barcode <- tcga_clin_data$bcr_patient_barcode
    #tcga_clin_data <- tcga_maf_obj@clinical.data
    tcga_pheno_columns <- c("Tumor_Sample_Barcode","ajcc_pathologic_stage","age_at_diagnosis","gender","race","vital_status","tissue_or_organ_of_origin")
    matched_order=1:nrow(tcga_clin_data)
    if (!is.null(plotdata)) {
      matched_order=match(colnames(plotdata), tcga_clin_data$Tumor_Sample_Barcode, nomatch=0)
    }
    tcga_anno_data <- tcga_clin_data[matched_order,tcga_pheno_columns]
    cancerCode <- paste0(unique(tcga_clin_data$disease), collapse=",")
    tcga_anno_data$Dataset <- cancerCode
    tcgaClinicalAnnotationDataList[[cancerCode]]<-tcga_anno_data
    #anno_data <- tcga_anno_data
  }
  anno_data<-do.call("rbind",tcgaClinicalAnnotationDataList)
  stages=sort(unique(anno_data$ajcc_pathologic_stage))
  #stage_colors <- setNames(brewer.pal(n = length(stages), name = "Reds"), stages)
  stage_colors <- setNames(rainbow(length(stages)), stages)

  anno_data$age_at_diagnosis <- as.numeric(as.character(anno_data$age_at_diagnosis))
  age_range=round(range(anno_data$age_at_diagnosis, na.rm = TRUE),-1)
  age_color_length=10
  age_breaks=round(seq(age_range[1], age_range[2], length.out=age_color_length),0)
  age_color_vals=colorRampPalette(c("lightblue1","royalblue1","navy"))(age_color_length)
  age_colors=colorRamp2(age_breaks, age_color_vals)

  gender_colors=c(female="hotpink", male="cornflowerblue")

  races=sort(unique(anno_data$race))
  #race_colors <- setNames(rev(brewer.pal(n = length(races), name = "Set1")), races)
  race_colors <- setNames(rev(heat.colors(length(races))), races)

  statuses=sort(unique(anno_data$vital_status))
  vitstat_colors <- c(Alive="darkgreen",Dead="darkred")

  tissues=sort(unique(anno_data$tissue_or_organ_of_origin))
  #tissue_colors <- setNames(brewer.pal(n = length(tissues), name = "Dark2"), tissues)
  tissue_colors <- setNames(topo.colors(length(tissues)), tissues)

  # dataset_colors <- setNames(c("mediumorchid1","darkolivegreen1"),
  dataset_colors <- setNames(c("grey30","darkolivegreen1"),
                             c(cancerCode, "Other"))

  anno_colors <- setNames(list(stage_colors, age_colors, gender_colors, race_colors, vitstat_colors, tissue_colors, dataset_colors),
                          setdiff(colnames(anno_data),"Tumor_Sample_Barcode"))


  mycols <- which(!colnames(anno_data) %in% c("Tumor_Sample_Barcode"))
  anno_df <- anno_data[,mycols]
  myanno <- ComplexHeatmap::HeatmapAnnotation(df=anno_df,col = anno_colors)

  return(list(colorList=anno_colors, annodata=anno_data, HManno=myanno))

}

#' Makes reasonable colors for some TCGA clinical annoations
#' @description This will return a list of colors that can be used with TCGA clinical annotations
#' @author Mayank Tandon, Ashish Jain
#' @param ageRange The range of patient's age to generate color vectors
#' @export
#' @return A list containing the TCGA clinical annotations
#'
#' @examples
#' library(MAFDash)
#' colorList <- getTCGAClinicalColors()
getTCGAClinicalColors <- function(ageRange=c(0,100)) {
  #require(RColorBrewer)
  #suppressPackageStartupMessages(require(circlize))
  # tcga_pheno_columns <- c("Tumor_Sample_Barcode","ajcc_pathologic_stage","age_at_diagnosis","gender","race","vital_status","tissue_or_organ_of_origin")
  tcga_pheno_columns <- c("ajcc_pathologic_stage","age_at_diagnosis","gender","race","vital_status")

  # if (is.null(stages)) {
    stages <- paste("Stage",c("X",paste(rep(c("I","II","III","IV"), each=4),c("","A","B","C"), sep="")))
  # }
  #stage_colors <- setNames(colorRampPalette(brewer.pal(n = 5, name = "Reds"))(length(stages)), stages)
  stage_colors <- setNames(colorRampPalette(rainbow(5))(length(stages)), stages)

  # browser()
  ageRange=sort(ageRange)
  age_range=round(c(max(c(0,ageRange[1]-4)),ageRange[2]+4),-1)
  age_color_length=10
  age_breaks=round(seq(age_range[1], age_range[2], length.out=age_color_length),0)
  age_color_vals=colorRampPalette(c("lightblue1","royalblue1","navy"))(age_color_length)
  age_colors=colorRamp2(age_breaks, age_color_vals)

  gender_colors=c(female="hotpink", male="cornflowerblue")
  vitstat_colors <- c(Alive="darkgreen",Dead="darkred","Not Reported"="grey80")

  races=c("american indian or alaska native","asian","black or african american","native hawaiian or other pacific islander","white","not allowed to collect","not reported","other","unknown")
  #race_colors <- setNames(brewer.pal(n = length(races), name = "Set1"), races)
  race_colors <- setNames(heat.colors(length(races)), races)


  # tissues=sort(unique(anno_data$tissue_or_organ_of_origin))
  # tissue_colors <- setNames(brewer.pal(n = length(tissues), name = "Dark2"), tissues)

  # dataset_colors <- setNames(c("mediumorchid1","darkolivegreen1"),
  # dataset_colors <- setNames(c("grey30","darkolivegreen1"),
  #                            c(cancerCode, "Other"))

  anno_colors <- setNames(list(stage_colors, age_colors, gender_colors, race_colors, vitstat_colors),
                          tcga_pheno_columns)

  return(anno_colors)
}
