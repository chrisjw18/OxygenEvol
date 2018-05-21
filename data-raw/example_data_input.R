#load example datasets to global environment
example_fluorwin_data<-system.file("extdata", "example_fluorwin_data.txt", package = "OxygenEvol")
example_oxygen_data<-system.file("extdata", "example_oxygen_data.txt", package = "OxygenEvol")

devtools::use_data(example_fluorwin_data)
devtools::use_data(example_oxygen_data)
