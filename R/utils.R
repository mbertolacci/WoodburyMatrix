.recycle_vector_to <- function(x, to_length) {
  x[
    ((seq_len(to_length) - 1) %% length(x)) + 1
  ]
}
