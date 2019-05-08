
# 
Created with
```
devtools::create("ratools")
```
Compile with:
```
devtools::document("ratools")
devtools::install("ratools")
#devtools::use_package("dplyr")
```

# Notes

magrittr pipe imported in `utilities.R`:
```
#' @importFrom magrittr %>%
#' @export
magrittr::`%>%`
```
(although the export line is not used)

https://stackoverflow.com/questions/31328404/how-to-call-r-script-from-another-r-script-both-in-same-package
https://stackoverflow.com/questions/27947344/r-use-magrittr-pipe-operator-in-self-written-package