getCodes <- function(x, idx = 1:length(codes)) {
  codes <- x$codes
  if (length(idx) == 1) {
    codes[[idx]]
  } else {
    codes[idx]
  }
}
