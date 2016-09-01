#' summaryTable
#'
#' \code{summaryTable} gives basic statistics for patient groupings in IDEA data extract
#'
#' @param data Individual patient records
#' @param dir filename to save output
#' @return dataframe

summaryTable <- function(datat, dir=NULL){
    ## summaryTable(data)

    long.table <- function(x){
        if(is.factor(x)){x <- droplevels(x)}
        d <- as.data.frame(table(x))
        d[,1] <- as.character(d[,1])
        d}

    out <- rbind(
        c("N", nrow(data)),
        c("Age (years)", paste(round(mean(datat$PatientAge, na.rm=T)), " (", round(sqrt(var(datat$PatientAge, na.rm=T))), ")", sep="")),
        c("Gender (male)", sum(datat$Sex=="Male", na.rm=T)),
        c("Ethnic group",""),
        long.table(datat$Ethnclass),
        c("Region",""),
        long.table(datat$Region),
        c("HIV positive", sum(datat$HIVpos==TRUE, na.rm=T)),

        c("Symptoms",""),
        long.table(datat$numSymptoms),
        c("Cough",sum(datat$Cough==T, na.rm=T)),
        c("Fever", sum(datat$Fever==T, na.rm=T)),
        c("Night sweats",sum(datat$Ngtsweat==T, na.rm=T)),
        c("Lethargy",sum(datat$Leth==T, na.rm=T)),
        c("Weight loss",sum(datat$Wghtloss==T, na.rm=T)),
        c("Haemop",sum(datat$Haemop==T, na.rm=T)),

        c("Dosanjh category",""),
        long.table(datat$Dosanjh),
        c("Diagnosis outcome",""),
        # long.table(datat$Diagoutcome),
        long.table(datat$DiagoutcomeGrouped),
        #     as.datat.frame(prop.table(table(datat$Diagoutcome)))

        c("Times",""),
        c("Treatment to culture test", paste(round(mean(datat$DrugCult_diff, na.rm=T)), " (", round(sqrt(var(datat$DrugCult_diff, na.rm=T))), ")", sep="")),
        c("First test to start treatment", paste(round(mean(datat$testDrug_diff, na.rm=T)), " (", round(sqrt(var(datat$testDrug_diff, na.rm=T))), ")", sep="")),
        c("First test to confirmed diagnosis", paste(round(mean(datat$testDiagCon_diff, na.rm=T)), " (", round(sqrt(var(datat$testDiagCon_diff, na.rm=T))), ")", sep="")),
        c("First test to culture", paste(round(mean(datat$testCult_diff, na.rm=T)), " (", round(sqrt(var(datat$testCult_diff, na.rm=T))), ")", sep="")),
        c("Total treatment duration", paste(round(mean(datat$TBDrug_diff, na.rm=T)), " (", round(sqrt(var(datat$TBDrug_diff, na.rm=T))), ")", sep=""))
        #     c("Consent to first treatment", paste(round(mean(datat$testDrug_diff, na.rm=T)), " (", round(sqrt(var(datat$testDrug_diff, na.rm=T))), ")", sep="")),
        #     c("Consent to confirmed diagnosis", paste(round(mean(datat$testDiagCon_diff, na.rm=T)), " (", round(sqrt(var(datat$testDiagCon_diff, na.rm=T))), ")", sep="")),
        #     c("Consent to culture", paste(round(mean(datat$testCult_diff, na.rm=T)), " (", round(sqrt(var(datat$testCult_diff, na.rm=T))), ")", sep=""))

    )

    colnames(out) <- c(" "," ")

    if(!is.null(dir)){
      write.csv(out, file=dir)}

    out
}
