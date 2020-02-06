##
## Batch generate all of the documents and plots for the HTA report
## saves as separate files

library(knitr)
library(pander)
library(rmarkdown)
library(IDEAdectree)

# http://www.r-statistics.com/2013/03/write-ms-word-document-using-r-with-as-little-overhead-as-possible/
# if(!require(installr)) { install.packages("installr"); require(installr)} #load / install+load installr
# install.pandoc()


filenames.full <- list.files(path = "./docs", pattern = "\\.Rmd$", full.names = TRUE)
filenames <- list.files(path = "./docs", pattern = "\\.Rmd$", full.names = FALSE)
filenames <- sub("^([^.]*).*", "\\1", filenames)


# knit2html(filenames.full[2], output=".")    #deprecated
# rmarkdown::render(input=filenames.full[2], output_format="all")

# system(paste0("pandoc -o ", filenames[2], ".docx ", filenames[2], ".md"))

## this way saves figures separately too ##
for (i in 1:length(filenames)) {
  knitr::knit(input = filenames.full[i],
              output = paste0("docs/", filenames[i], ".md"))
  # rmarkdown::render(input=paste0("docs/", filenames[i], ".md"), output_format="all")
}





