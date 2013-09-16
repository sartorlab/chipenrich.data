example_datasets = function() {
  message("Example ChIP-seq peak datasets:");
  print(grep("peaks_",data(package="chipenrich.data")$results[,3],value=T));

  message("Use 'data(name of dataset)' to load it.");
}
