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
txs.mean <- rnbinom(length(txs.name), size = 1, prob = 0.001)/100 # randomize expression
txs.mean <- txs.mean * (txs.length/100) # take into consideration tx length

# get within and between sd
set.seed(1987)
within.sd.mat <- matrix(rbinom(length(txs.length)*2, size = 100, prob = 0.2)/50, ncol = 2)
set.seed(1987)
btw.sd <- rnbinom(length(txs.length), size = 1, prob = 0.1)
btw.sd <- btw.sd/(2*max(btw.sd))

sim.counts <- do.call(rbind, lapply(seq_len(length(txs.name)), function(x){
    set.seed(1987)
    group.mean <- rnorm(2, mean = txs.mean[x], sd = btw.sd[x]*txs.mean[x])
    do.call(cbind, lapply(seq_len(length(group.mean)), function(y){
        set.seed(1987)
        matrix(rnorm(2, mean = group.mean[y], sd = within.sd.mat[x,y]*group.mean[y]),ncol = 2)
    }))
}))

countsSample <- do.call(cbind, lapply(seq_len(length(samples)), function(x){
    ceiling(lib.sizes[x]*sim.counts[,x]/sum(sim.counts[,x]))
}))
colnames(countsSample) <- samples
rownames(countsSample) <- txs.name
usethis::use_data(countsSample, overwrite = TRUE)

testSamples <- data.frame(samples = samples,
           group = rep(c("SampleA", "SampleB"), each = 2),
           row.names = samples)
