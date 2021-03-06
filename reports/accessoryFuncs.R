## http://www.r-bloggers.com/identifying-records-in-data-frame-a-that-are-not-contained-in-data-frame-b-%E2%80%93-a-comparison/
thisnotinthat <- function(x.1,x.2,...){
  x.1p <- do.call("paste", x.1)
  x.2p <- do.call("paste", x.2)
  x.1[! x.1p %in% x.2p, ]
}

# display datatable function
viewDataTable <- function(dat){
  
  cols2keep <- grep('^transcript_id$|^set$|^family_count$|^db$|^ds$', colnames(dat), invert = T, value = T)
  dat <- dat[,colnames(dat) %in% cols2keep]

  cols2crop <- grep('family_genotypes|ref|alt', colnames(dat))-1
  cols2hide <- grep('Description|Variants|Perc|^chr$|^pos$|rsid|^ref$|^alt$|^gene$|gene_freq|
                    impact_severity|rvis|depths|gt_quals|family_genotypes|denovo_prob|
                    is_exonic|is_splicing|functiongvs|exac_ac_all',
                    colnames(dat), invert = TRUE)-1
  
  DT::datatable(dat,
                rownames = FALSE,
                escape = FALSE,
                extensions = c('Buttons'),
                selection = "single",
                filter = "bottom",
                options = list(
                  columnDefs = list(list(visible = FALSE, targets = cols2hide),
                                    list(targets = cols2crop, render = JS("function(data, type, row, meta) {", 
                                                                          "return type === 'display' && data.length > 10 ?", 
                                                                          "'<span title=\"' + data + '\">' + data.substr(0, 10) + '...</span>' : data;", 
                                                                          "}"))),
                  dom = 'Bfrtip',
                  buttons = list(list(extend = "colvis", 
                                      collectionLayout = "fixed four-column", 
                                      prefixButtons = list('colvisRestore')),
                                 'pageLength', 'copy', 'print',
                                 list(extend = "collection",
                                      buttons = c('csv', 'excel', 'pdf'),
                                      text = 'Download'
                                 )),
                  searchHighlight = TRUE,
                  lengthMenu = list(c(5, 10, 15, 20, 25, -1), c('5', '10', '15', '20', '25', 'All')),
                  initComplete = DT::JS("function(settings, json) {",
                                        "$(this.api().table().header()).css({'background-color': '#003366', 'color': '#fff'});",
                                        "}"),
                  scrollX = TRUE
                ))
}