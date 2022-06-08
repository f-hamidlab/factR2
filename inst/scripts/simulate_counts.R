data('factRsample')

samples <- c("SampleA.Rep1", "SampleA.Rep2",
             "SampleB.Rep1", "SampleB.Rep2")

# simulate library size
set.seed(1987)
lib.sizes <- sample(10000:14000, size = length(samples))

# get list of transcripts
txsData <- featureData(factRsample[["transcript"]])
txs.length <- txsData$width
txs.name <- rownames(txsData)

# simulate mean counts per tx
set.seed(1987)
txs.mean <- rnbinom(length(txs.name), size = 1, prob = 0.001) # randomize expression
txs.mean <- txs.mean * (txs.length/100) # take into consideration tx length

