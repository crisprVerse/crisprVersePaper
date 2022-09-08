guides <- guidesCD46
spacers <- guides$spacer
write.table(spacers,
          quote=FALSE,
          row.names=FALSE,
          col.names=FALSE,
          file="objects/spacers_guides_cd46.txt")

guides  <- guidesCD55
spacers <- guidesCD55$spacer
write.table(spacers,
          quote=FALSE,
          row.names=FALSE,
          col.names=FALSE,
          file="objects/spacers_guides_cd55.txt")
