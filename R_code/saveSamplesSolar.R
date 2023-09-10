######################################################################################

# Copyright 2023, Nicolau Andrés-Thió

# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:

# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.

# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.

######################################################################################

# This file is used to save the objective function values of the SOLAR simulator
# for different sample plans and different fidelities.
# This is done so that features can be calculated without resampling the simulator
# multiple times at the same location and for the same fidelity, as this is very
# inefficient. New users should not be required to run this as the objective function
# values are stored in the data folder.


library(Rcpp)
sourceCpp('cpp_code/rFeatureAnalysis.cpp')

args = commandArgs(trailingOnly=TRUE)
seed = as.numeric(args[[1]])

lowFiSizes <- seq(4, 20, by = 4)

# Create dataframe with the results
for(lowFiSize in lowFiSizes){
  # Create data where the values will be saved
  data <- data.frame(matrix(nrow = 0, ncol = 15))
  colnames(data) <- c(paste0("x", 1:5), paste0("fid", c("1.00", "0.90", "0.80", "0.70", "0.60", "0.50", "0.40", "0.30", "0.20", "0.10")))
  lowFiBudget <- lowFiSize * 5
  print(paste0("Working on size ", lowFiBudget, " for seed ", seed))
  fmin = 0
  fmax = 0
  # First get the sample
  samples <- functionSample("SOLAR0.90", seed, lowFiBudget, lowFiBudget, TRUE, fmin, fmax)[[1]]
  # Now want to get the value for each of the fidelities
  nrow <- 0
  for(sample in samples){
    nrow <- nrow + 1
    vals <- sample
    for(fid in c("1.00", "0.90", "0.80", "0.70", "0.60", "0.50", "0.40", "0.30", "0.20", "0.10")){
      print(fid)
      vals <- c(vals, sampleFunction(paste0("SOLAR",fid,seed), sample, 1, TRUE, fmin, fmax))
    }
    data[nrow, ] <- vals
    write.table(data, paste0("data/misc/SOLARn", lowFiBudget, "s",seed,".txt"), quote = FALSE, row.names = FALSE)
  }
  # Have the data, now want to write it to a file for future use
}
 