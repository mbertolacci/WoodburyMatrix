#' Virtual class for Woodbury identity matrices
#'
#' The \code{WoodburyMatrix} is a virtual class, contained by both
#' \code{GWoodburyMatrix} (for general matrices) and \code{SWoodburyMatrix}
#' (for symmetric matrices). See \link{WoodburyMatrix} for construction
#' of these classes. The methods available for these classes are described
#' below; see also the \link{solve} methods. This class is itself a subclass of
#' \code{\link[Matrix:Matrix-class]{Matrix}}, so basic matrix methods like
#' \code{nrow}, \code{ncol}, \code{dim} and so on also work.
#'
#' @slot A n x n subclass of \code{\link[Matrix:Matrix-class]{Matrix}}
#' (\code{GWoodburyMatrix}) or \code{\link[Matrix:symmetricMatrix-class]{symmetricMatrix}}
#' (\code{SWoodburyMatrix}).
#' @slot B p x p subclass of \code{\link[Matrix:Matrix-class]{Matrix}}
#' (\code{GWoodburyMatrix}) or \code{\link[Matrix:symmetricMatrix-class]{symmetricMatrix}}
#' (\code{SWoodburyMatrix}).
#' @slot U n x p subclass of \code{\link[Matrix:Matrix-class]{Matrix}} (only for
# \code{GWoodburyMatrix})
#' @slot V p x m subclass of \code{\link[Matrix:Matrix-class]{Matrix}} (only for
# \code{GWoodburyMatrix})
#' @slot X n x p subclass of \code{\link[Matrix:Matrix-class]{Matrix}} (only for
# \code{SWoodburyMatrix})
#' @slot O p x p subclass of \code{\link[Matrix:Matrix-class]{Matrix}}
#' @param x \code{WoodburyMatrix} object
#' @param object \code{WoodburyMatrix} object
#' @param y Matrix or vector
#' @seealso \link{WoodburyMatrix} for object construction,
#' \link[Matrix:Matrix-class]{Matrix} (the parent of this class).
#' @export
setClass(
  'WoodburyMatrix',
  contains = c('Matrix', 'VIRTUAL')
)

#' @rdname WoodburyMatrix-class
#' @description Sub-class representing a generic matrix.
#' @export
setClass(
  'GWoodburyMatrix',
  contains = 'WoodburyMatrix',
  slots = list(
    A = 'Matrix',
    B = 'Matrix',
    U = 'Matrix',
    V = 'Matrix',
    O = 'Matrix'
  )
)

setClassUnion('symMatrix', c('diagonalMatrix', 'symmetricMatrix'))

#' @rdname WoodburyMatrix-class
#' @description Sub-class representing a symmetric matrix. Also subclasses
#' \link[Matrix:symmetricMatrix-class]{symmetricMatrix}.
#' @export
setClass(
  'SWoodburyMatrix',
  contains = c('WoodburyMatrix', 'symmetricMatrix'),
  slots = list(
    A = 'symMatrix',
    B = 'symMatrix',
    X = 'Matrix',
    O = 'symMatrix'
  )
)
