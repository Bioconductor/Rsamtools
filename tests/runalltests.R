suppressMessages(library("Rsamtools"))
suppressMessages(library("RUnit"))

options(warn=1)

dirs <- 'unit'

testFilePat <- ".*_test\\.R$"

allSuite <- defineTestSuite(name="allSuite",
                            dirs=dirs,
                            testFileRegexp=testFilePat,
                            rngKind="default",
                            rngNormalKind="default")

results <- capture.output(runTestSuite(allSuite))
if (interactive())
    cat(paste(results, collapse="\n"), "\n")

q(runLast=FALSE)
