capture.output(
  utils::sessionInfo(),
  file = file.path(output_dir, "R_sessionInfo.txt")
)

package_names <- c(
  "limma", "DESeq2", "ggplot2", "dplyr", "tidyr",
  "patchwork", "scales"
)
package_versions <- data.frame(
  Package = package_names,
  Installed = vapply(
    package_names,
    requireNamespace,
    logical(1),
    quietly = TRUE
  ),
  Version = vapply(
    package_names,
    function(package) {
      if (!requireNamespace(package, quietly = TRUE)) return(NA_character_)
      as.character(utils::packageVersion(package))
    },
    character(1)
  ),
  stringsAsFactors = FALSE
)
write_result(package_versions, "R_package_versions.csv")

