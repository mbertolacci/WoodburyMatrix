setValidity('GWoodburyMatrix', function(object) {
  A <- object@A
  B <- object@B
  U <- object@U
  V <- object@V
  if (nrow(A) != ncol(A)) {
    'A is not square'
  } else if (nrow(A) != nrow(U)) {
    'A and U are not compatible'
  } else if (nrow(A) != ncol(V)) {
    'A and V are not compatible'
  } else if (ncol(U) != nrow(B)) {
    'U and B are not compatible'
  } else if (ncol(B) != nrow(V)) {
    'B and V are not compatible'
  } else {
    TRUE
  }
})

setValidity('SWoodburyMatrix', function(object) {
  A <- object@A
  B <- object@B
  X <- object@X
  if (nrow(A) != nrow(X)) {
    'A and X are not compatible'
  } else if (ncol(X) != nrow(B)) {
    'X and B are not compatible'
  } else {
    TRUE
  }
})
