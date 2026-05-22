#' Compute body length and speed summaries for one trial
#'
#' Uses markers "a" and "p" as head and tail to compute body length and
#' center-of-body speed summaries.
#'
#' @param x One element from \code{process_trk_file()}, i.e. a list with
#'   components `meta` and `coords`.
#'
#' @return A one-row data.frame with body length and speed metrics.
#' @export
compute_body_stats <- function(x) {
  meta <- x$meta

  out <- data.frame(
    file             = meta$file,
    instar           = meta$instar,
    nr               = meta$nr,
    id               = meta$id,
    mass             = meta$mass,
    length_unit      = meta$length_unit,
    delta_t          = meta$delta_t,
    body_length_mean = NA_real_,
    body_length_sd   = NA_real_,
    speed_mean       = NA_real_,
    speed_median     = NA_real_,
    speed_max        = NA_real_,
    speed_sd         = NA_real_
  )

  try({
    trj_a <- x$coords$a$trj[c("x", "y", "frame", "time")]
    trj_p <- x$coords$p$trj[c("x", "y", "frame", "time")]

    trj <- merge(trj_a, trj_p, by = "frame", suffixes = c("1", "2"), all = FALSE)

    l <- sqrt((trj$x2 - trj$x1)^2 + (trj$y2 - trj$y1)^2)
    out$body_length_mean <- mean(l, na.rm = TRUE)
    out$body_length_sd   <- stats::sd(l, na.rm = TRUE)

    midpoint <- trajr::TrajFromCoords(
      data.frame(
        x    = (trj$x1 + trj$x2) / 2,
        y    = (trj$y1 + trj$y2) / 2,
        time = trj$time1,
        frame = trj$frame
      ),
      timeCol = "time"
    )

    derivs_mid <- trajr::TrajDerivatives(
      trajr::TrajSmoothSG(midpoint, p = 3, n = 5)
    )

    v <- derivs_mid$speed

    out$speed_mean   <- mean(v, na.rm = TRUE)
    out$speed_median <- stats::median(v, na.rm = TRUE)
    out$speed_max    <- max(v, na.rm = TRUE)
    out$speed_sd     <- stats::sd(v, na.rm = TRUE)
  }, silent = TRUE)

  out
}

#' Summarise body kinematics for multiple trials
#'
#' Applies \code{compute_body_stats()} to a list of trial objects.
#'
#' @param trials List of objects as returned by \code{process_trk_file()}.
#'
#' @return A data.table with one row per trial.
#' @import data.table
#' @export
summarise_body <- function(trials) {
  stats <- lapply(trials, compute_body_stats)
  stats_dt <- data.table::rbindlist(stats, fill = TRUE)
  stats_dt$body_length_per_sec <- stats_dt$speed_mean / stats_dt$body_length_mean
  stats_dt
}

#' Rotate coordinates into body frame
#'
#' Applies a rotation matrix to translate/rotate coordinates such that the body
#' long axis is aligned, using per-frame angles.
#'
#' @param x Numeric vector of x coordinates (relative to COM).
#' @param y Numeric vector of y coordinates (relative to COM).
#' @param theta Numeric vector of angles (radians), one per row.
#' @param point Character scalar: marker name.
#' @param frame Integer vector of frame indices.
#' @param time Numeric vector of times.
#'
#' @return A data.frame with rotated coordinates `x`, `y` and metadata.
#' @export
rotate_body_frame <- function(x, y, theta, point, frame, time) {
  xr <- x * cos(-theta) - y * sin(-theta)
  yr <- x * sin(-theta) + y * cos(-theta)
  data.frame(
    x     = xr,
    y     = yr,
    frame = frame,
    point = point,
    theta = theta,
    time  = time
  )
}

#' Find full sign intervals in a time series
#'
#' Detects contiguous intervals where a signal stays strictly positive or
#' negative between zero-crossings.
#'
#' @param x Numeric vector (e.g. velocity).
#' @param crossing Direction of crossing defining the interval:
#'   `"up"` (from negative to positive) or `"down"` (positive to negative).
#'
#' @return A matrix with columns `start` and `end` (indices into `x`),
#'   possibly with zero rows if no intervals are found.
#' @export
find_full_intervals <- function(x, crossing = c("up", "down")) {
  crossing <- match.arg(crossing)
  up_crossings   <- which(diff(sign(x)) == 2)
  down_crossings <- which(diff(sign(x)) == -2)

  up_crossings   <- sort(up_crossings)
  down_crossings <- sort(down_crossings)

  intervals <- list()

  if (crossing == "up") {
    for (up in up_crossings) {
      downs_after <- down_crossings[down_crossings > up]
      if (length(downs_after) > 0) {
        intervals[[length(intervals) + 1]] <- c(start = up + 1, end = downs_after[1])
      }
    }
    if (length(intervals) > 0) {
      intervals <- do.call(rbind, intervals)
      intervals <- intervals[
        apply(intervals, 1, function(iv) all(x[iv[1]:iv[2]] > 0)),
        ,
        drop = FALSE
      ]
    } else {
      intervals <- matrix(ncol = 2, nrow = 0,
                          dimnames = list(NULL, c("start", "end")))
    }
  } else {
    for (down in down_crossings) {
      ups_after <- up_crossings[up_crossings > down]
      if (length(ups_after) > 0) {
        intervals[[length(intervals) + 1]] <- c(start = down + 1, end = ups_after[1])
      }
    }
    if (length(intervals) > 0) {
      intervals <- do.call(rbind, intervals)
      intervals <- intervals[
        apply(intervals, 1, function(iv) all(x[iv[1]:iv[2]] < 0)),
        ,
        drop = FALSE
      ]
    } else {
      intervals <- matrix(ncol = 2, nrow = 0,
                          dimnames = list(NULL, c("start", "end")))
    }
  }
  intervals
}

#' Compute leg protraction and retraction times
#'
#' For each marker in one trial, rotates trajectories into the body frame,
#' computes vertical velocity, and extracts protraction and retraction
#' intervals based on sign of velocity.
#'
#' @param trial One element as returned by \code{process_trk_file()}.
#'
#' @return The input trial object with added components
#'   `coords[[marker]]$protraction_time` and `$retraction_time` (vectors of durations).
#' @import data.table
#' @export
add_leg_timing <- function(trial) {
  meta <- trial$meta

  # head/tail for COM and body orientation
  trj_a <- trial$coords$a$trj[c("x", "y", "frame", "time")]
  trj_p <- trial$coords$p$trj[c("x", "y", "frame", "time")]
  trj   <- merge(trj_a, trj_p, by = "frame", suffixes = c("1", "2"), all = FALSE)

  R0 <- data.frame(
    x0    = (trj$x1 + trj$x2) / 2,
    y0    = (trj$y1 + trj$y2) / 2,
    frame = trj$frame,
    time  = trj$time1
  )

  delta_x <- trj$x2 - trj$x1
  delta_y <- trj$y2 - trj$y1
  theta   <- atan2(delta_y, delta_x)
  R0$theta <- theta + pi / 2

  point_names <- names(trial$coords)

  for (pt in point_names) {
    trj_pt <- trial$coords[[pt]]$trj[c("x", "y", "frame", "time")]
    # merge only by frame
    tmp    <- merge(trj_pt, R0, by = "frame")
    tmp$x  <- tmp$x.x - tmp$x0
    tmp$y  <- tmp$y.x - tmp$y0
    tmp$time <- tmp$time.x

    rot <- rotate_body_frame(
      x     = tmp$x,
      y     = tmp$y,
      theta = tmp$theta,
      point = pt,
      frame = tmp$frame,
      time  = tmp$time
    )

    v_y <- c(
      NA_real_,
      diff(gsignal::sgolayfilt(rot$y, p = 3, n = 7)) / (meta$delta_t / 1000)
    )

    pro <- as.data.frame(find_full_intervals(v_y, crossing = "up"))
    pro$duration <- (pro$end - pro$start) * meta$delta_t / 1000

    ret <- as.data.frame(find_full_intervals(v_y, crossing = "down"))
    ret$duration <- (ret$end - ret$start) * meta$delta_t / 1000

    trial$coords[[pt]]$protraction_time <- pro$duration
    trial$coords[[pt]]$retraction_time  <- ret$duration
  }

  trial
}

#' Compute leg protraction and retraction times for one trial
#'
#' For each marker in one trial, this function:
#' \itemize{
#'   \item Computes body COM and orientation from head ("a") and tail ("p").
#'   \item Rotates marker trajectories into the body frame.
#'   \item Computes vertical velocity via Savitzky–Golay filtering.
#'   \item Extracts protraction and retraction intervals from velocity
#'         zero-crossings and stores their durations.
#' }
#'
#' Protraction and retraction times are stored as vectors of durations (in seconds)
#' in \code{trial$coords[[marker]]$protraction_time} and
#' \code{trial$coords[[marker]]$retraction_time}. These vectors may have
#' length zero if no intervals are detected.
#'
#' @param trial One element as returned by \code{process_trk_file()}, i.e.
#'   a list with components \code{meta} and \code{coords}.
#'
#' @return The input trial object with added components
#'   \code{coords[[marker]]$protraction_time} and
#'   \code{coords[[marker]]$retraction_time}.
#' @export
add_leg_timing <- function(trial) {
  meta <- trial$meta

  # --- 1. Head/tail trajectories for COM and body orientation ----
  trj_a <- trial$coords$a$trj[, c("x", "y", "frame", "time")]
  trj_p <- trial$coords$p$trj[, c("x", "y", "frame", "time")]

  # Merge by frame only; time is derived from 'a'
  trj <- merge(trj_a, trj_p, by = "frame", suffixes = c("1", "2"), all = FALSE)

  R0 <- data.frame(
    x0    = (trj$x1 + trj$x2) / 2,
    y0    = (trj$y1 + trj$y2) / 2,
    frame = trj$frame,
    time  = trj$time1
  )

  delta_x <- trj$x2 - trj$x1
  delta_y <- trj$y2 - trj$y1
  theta   <- atan2(delta_y, delta_x)
  R0$theta <- theta + pi / 2

  # --- 2. For each marker: rotate into body frame and extract intervals ----
  point_names <- names(trial$coords)

  for (pt in point_names) {
    trj_pt <- trial$coords[[pt]]$trj[, c("x", "y", "frame", "time")]

    # Align marker with COM/orientation by frame
    tmp <- merge(trj_pt, R0, by = "frame")
    # After merge: x.x, y.x (marker); x0, y0 (COM); time.x (marker time)
    tmp$x_rel <- tmp$x.x - tmp$x0
    tmp$y_rel <- tmp$y.x - tmp$y0
    tmp$time  <- tmp$time.x

    rot <- rotate_body_frame(
      x     = tmp$x_rel,
      y     = tmp$y_rel,
      theta = tmp$theta,
      point = pt,
      frame = tmp$frame,
      time  = tmp$time
    )

    # vertical velocity in body frame (mm/s)
    dt_s <- meta$delta_t / 1000
    v_y <- c(
      NA_real_,
      diff(gsignal::sgolayfilt(rot$y, p = 3, n = 7)) / dt_s
    )

    # --- protraction (upward crossings, v_y > 0) ---
    pro_mat <- find_full_intervals(v_y, crossing = "up")
    if (nrow(pro_mat) > 0) {
      pro <- as.data.frame(pro_mat)
      pro$duration <- (pro$end - pro$start) * dt_s
      pro_duration <- pro$duration
    } else {
      pro_duration <- numeric(0)
    }

    # --- retraction (downward crossings, v_y < 0) ---
    ret_mat <- find_full_intervals(v_y, crossing = "down")
    if (nrow(ret_mat) > 0) {
      ret <- as.data.frame(ret_mat)
      ret$duration <- (ret$end - ret$start) * dt_s
      ret_duration <- ret$duration
    } else {
      ret_duration <- numeric(0)
    }

    trial$coords[[pt]]$protraction_time <- pro_duration
    trial$coords[[pt]]$retraction_time  <- ret_duration
  }

  trial
}


#' Rotate coordinates into body frame
#'
#' @param x Numeric vector of x coordinates (relative to COM).
#' @param y Numeric vector of y coordinates (relative to COM).
#' @param theta Numeric vector of angles (radians), same length as x/y.
#' @param point Character scalar: marker name.
#' @param frame Integer vector of frame indices.
#' @param time Numeric vector of times.
#'
#' @return A data.frame with rotated coordinates and metadata.
#' @export
rotate_body_frame <- function(x, y, theta, point, frame, time) {
  n <- length(x)
  if (!all(lengths(list(y, theta, frame, time)) == n)) {
    stop("rotate_body_frame: input vectors must have equal length")
  }
  if (n == 0L) {
    return(data.frame(
      x = numeric(0),
      y = numeric(0),
      frame = integer(0),
      point = character(0),
      theta = numeric(0),
      time = numeric(0),
      stringsAsFactors = FALSE
    ))
  }

  xr <- x * cos(-theta) - y * sin(-theta)
  yr <- x * sin(-theta) + y * cos(-theta)

  data.frame(
    x     = xr,
    y     = yr,
    frame = frame,
    point = rep(point, n),
    theta = theta,
    time  = time,
    stringsAsFactors = FALSE
  )
}


#' Compute leg protraction and retraction times for one trial
#'
#' For each marker in one trial, this function:
#' - Computes body COM and orientation from head ("a") and tail ("p").
#' - Rotates marker trajectories into the body frame.
#' - Computes vertical velocity via Savitzky–Golay filtering.
#' - Extracts protraction and retraction intervals from velocity
#'   zero-crossings and stores their durations (seconds).
#'
#' @param trial One element as returned by process_trk_file().
#' @return The input trial with protraction_time and retraction_time
#'         vectors added to each marker in coords.
#' @export
add_leg_timing <- function(trial) {
  meta <- trial$meta

  # --- 1. Head/tail trajectories for COM and body orientation ----
  trj_a <- trial$coords$a$trj[, c("x", "y", "frame", "time")]
  trj_p <- trial$coords$p$trj[, c("x", "y", "frame", "time")]

  trj <- merge(trj_a, trj_p, by = "frame", suffixes = c("1", "2"), all = FALSE)

  R0 <- data.frame(
    x0    = (trj$x1 + trj$x2) / 2,
    y0    = (trj$y1 + trj$y2) / 2,
    frame = trj$frame,
    time  = trj$time1,
    theta = atan2(trj$y2 - trj$y1, trj$x2 - trj$x1) + pi / 2
  )

  point_names <- names(trial$coords)

  dt_s <- meta$delta_t / 1000

  for (pt in point_names) {
    trj_pt <- trial$coords[[pt]]$trj[, c("x", "y", "frame", "time")]

    # Align marker with COM/orientation by frame
    tmp <- merge(trj_pt, R0, by = "frame")

    # If no overlapping frames, skip this marker
    if (nrow(tmp) == 0) {
      trial$coords[[pt]]$protraction_time <- numeric(0)
      trial$coords[[pt]]$retraction_time  <- numeric(0)
      next
    }

    # After merge: x.x, y.x (marker), x0, y0 (COM), time.x, theta
    x_rel <- tmp$x.x - tmp$x0
    y_rel <- tmp$y.x - tmp$y0
    time  <- tmp$time.x
    theta <- tmp$theta

    rot <- rotate_body_frame(
      x     = x_rel,
      y     = y_rel,
      theta = theta,
      point = pt,
      frame = tmp$frame,
      time  = time
    )

    # vertical velocity in body frame (mm/s)
    v_y <- c(
      NA_real_,
      diff(gsignal::sgolayfilt(rot$y, p = 3, n = 7)) / dt_s
    )

    # --- protraction intervals (v_y > 0) ---
    pro_mat <- find_full_intervals(v_y, crossing = "up")
    if (nrow(pro_mat) > 0) {
      pro <- as.data.frame(pro_mat)
      pro$duration <- (pro$end - pro$start) * dt_s
      pro_duration <- pro$duration
    } else {
      pro_duration <- numeric(0)
    }

    # --- retraction intervals (v_y < 0) ---
    ret_mat <- find_full_intervals(v_y, crossing = "down")
    if (nrow(ret_mat) > 0) {
      ret <- as.data.frame(ret_mat)
      ret$duration <- (ret$end - ret$start) * dt_s
      ret_duration <- ret$duration
    } else {
      ret_duration <- numeric(0)
    }

    trial$coords[[pt]]$protraction_time <- pro_duration
    trial$coords[[pt]]$retraction_time  <- ret_duration
  }

  trial
}
add_leg_timing <- function(trial) {
  meta <- trial$meta

  # 1. Head & tail to get COM and orientation
  trj_a <- trial$coords$a$trj[, c("x", "y", "frame", "time")]
  trj_p <- trial$coords$p$trj[, c("x", "y", "frame", "time")]

  trj <- merge(trj_a, trj_p, by = "frame", suffixes = c("1", "2"), all = FALSE)

  if (nrow(trj) == 0L) {
    # no overlap between a and p; no timing possible
    for (pt in names(trial$coords)) {
      trial$coords[[pt]]$protraction_time <- numeric(0)
      trial$coords[[pt]]$retraction_time  <- numeric(0)
    }
    return(trial)
  }

  R0 <- data.frame(
    frame = trj$frame,
    time  = trj$time1,
    x0    = (trj$x1 + trj$x2) / 2,
    y0    = (trj$y1 + trj$y2) / 2,
    theta = atan2(trj$y2 - trj$y1, trj$x2 - trj$x1) + pi / 2
  )

  dt_s <- meta$delta_t / 1000
  point_names <- names(trial$coords)

  for (pt in point_names) {
    trj_pt <- trial$coords[[pt]]$trj[, c("x", "y", "frame", "time")]

    tmp <- merge(trj_pt, R0, by = "frame", all = FALSE)

    # No overlapping frames for this marker
    if (nrow(tmp) == 0L) {
      trial$coords[[pt]]$protraction_time <- numeric(0)
      trial$coords[[pt]]$retraction_time  <- numeric(0)
      next
    }

    # All these come from tmp, so have equal length nrow(tmp)
    x_rel <- tmp$x.x - tmp$x0
    y_rel <- tmp$y.x - tmp$y0
    time  <- tmp$time.x
    theta <- tmp$theta
    frame <- tmp$frame

    rot <- rotate_body_frame(
      x     = x_rel,
      y     = y_rel,
      theta = theta,
      point = pt,
      frame = frame,
      time  = time
    )

    # If rot ended up empty (shouldn’t, but be safe)
    if (nrow(rot) == 0L) {
      trial$coords[[pt]]$protraction_time <- numeric(0)
      trial$coords[[pt]]$retraction_time  <- numeric(0)
      next
    }

    v_y <- c(
      NA_real_,
      diff(gsignal::sgolayfilt(rot$y, p = 3, n = 7)) / dt_s
    )

    # Protraction (up)
    pro_mat <- find_full_intervals(v_y, crossing = "up")
    if (nrow(pro_mat) > 0) {
      pro <- as.data.frame(pro_mat)
      pro$duration <- (pro$end - pro$start) * dt_s
      pro_duration <- pro$duration
    } else {
      pro_duration <- numeric(0)
    }

    # Retraction (down)
    ret_mat <- find_full_intervals(v_y, crossing = "down")
    if (nrow(ret_mat) > 0) {
      ret <- as.data.frame(ret_mat)
      ret$duration <- (ret$end - ret$start) * dt_s
      ret_duration <- ret$duration
    } else {
      ret_duration <- numeric(0)
    }

    trial$coords[[pt]]$protraction_time <- pro_duration
    trial$coords[[pt]]$retraction_time  <- ret_duration
  }

  trial
}

rotate_body_frame <- function(x, y, theta, point, frame, time) {
  # Coerce to vectors
  x     <- as.numeric(x)
  y     <- as.numeric(y)
  theta <- as.numeric(theta)
  frame <- as.integer(frame)
  time  <- as.numeric(time)

  n <- min(length(x), length(y), length(theta), length(frame), length(time))

  if (n == 0L) {
    return(data.frame(
      x = numeric(0),
      y = numeric(0),
      frame = integer(0),
      point = character(0),
      theta = numeric(0),
      time = numeric(0),
      stringsAsFactors = FALSE
    ))
  }

  # If any vector is longer, truncate it; this avoids recycling,
  # but keeps everything aligned.
  x     <- x[seq_len(n)]
  y     <- y[seq_len(n)]
  theta <- theta[seq_len(n)]
  frame <- frame[seq_len(n)]
  time  <- time[seq_len(n)]

  xr <- x * cos(-theta) - y * sin(-theta)
  yr <- x * sin(-theta) + y * cos(-theta)

  data.frame(
    x     = xr,
    y     = yr,
    frame = frame,
    point = rep(point, n),
    theta = theta,
    time  = time,
    stringsAsFactors = FALSE
  )
}

add_leg_timing <- function(trial) {
  meta <- trial$meta

  trj_a <- trial$coords$a$trj[, c("x", "y", "frame", "time")]
  trj_p <- trial$coords$p$trj[, c("x", "y", "frame", "time")]

  trj <- merge(trj_a, trj_p, by = "frame", suffixes = c("1", "2"), all = FALSE)
  if (nrow(trj) == 0L) {
    for (pt in names(trial$coords)) {
      trial$coords[[pt]]$protraction_time <- numeric(0)
      trial$coords[[pt]]$retraction_time  <- numeric(0)
    }
    return(trial)
  }

  R0 <- data.frame(
    frame = trj$frame,
    time  = trj$time1,
    x0    = (trj$x1 + trj$x2) / 2,
    y0    = (trj$y1 + trj$y2) / 2,
    theta = atan2(trj$y2 - trj$y1, trj$x2 - trj$x1) + pi / 2
  )

  dt_s <- meta$delta_t / 1000
  point_names <- names(trial$coords)

  for (pt in point_names) {
    trj_pt <- trial$coords[[pt]]$trj[, c("x", "y", "frame", "time")]

    tmp <- merge(trj_pt, R0, by = "frame", all = FALSE)
    if (nrow(tmp) == 0L) {
      trial$coords[[pt]]$protraction_time <- numeric(0)
      trial$coords[[pt]]$retraction_time  <- numeric(0)
      next
    }

    x_rel <- tmp$x.x - tmp$x0
    y_rel <- tmp$y.x - tmp$y0
    theta <- tmp$theta
    frame <- tmp$frame
    time  <- tmp$time.x

    rot <- rotate_body_frame(
      x     = x_rel,
      y     = y_rel,
      theta = theta,
      point = pt,
      frame = frame,
      time  = time
    )

    if (nrow(rot) == 0L) {
      trial$coords[[pt]]$protraction_time <- numeric(0)
      trial$coords[[pt]]$retraction_time  <- numeric(0)
      next
    }

    v_y <- c(
      NA_real_,
      diff(gsignal::sgolayfilt(rot$y, p = 3, n = 7)) / dt_s
    )

    pro_mat <- find_full_intervals(v_y, crossing = "up")
    if (nrow(pro_mat) > 0) {
      pro <- as.data.frame(pro_mat)
      pro$duration <- (pro$end - pro$start) * dt_s
      pro_duration <- pro$duration
    } else {
      pro_duration <- numeric(0)
    }

    ret_mat <- find_full_intervals(v_y, crossing = "down")
    if (nrow(ret_mat) > 0) {
      ret <- as.data.frame(ret_mat)
      ret$duration <- (ret$end - ret$start) * dt_s
      ret_duration <- ret$duration
    } else {
      ret_duration <- numeric(0)
    }

    trial$coords[[pt]]$protraction_time <- pro_duration
    trial$coords[[pt]]$retraction_time  <- ret_duration
  }

  trial
}

#' Protraction / retraction intervals from speed
#'
#' Uses \code{trajr::TrajSpeedIntervals()} to detect intervals where
#' leg-tip speed is above or below a threshold, and interprets these
#' as swing (protraction) and stance (retraction) phases.
#'
#' @param trj A \code{trajr} trajectory object for one leg marker.
#' @param swing_threshold Speed threshold above which a leg is considered
#'   to be in swing (protraction). Units follow \code{trj}.
#' @param method Velocity differencing method passed to
#'   \code{trajr::TrajSpeedIntervals()}.
#'
#' @return A list with components:
#'   \describe{
#'     \item{protraction}{data.frame with intervals where speed > threshold.}
#'     \item{retraction}{data.frame with intervals where speed <= threshold.}
#'   }
#'   Each data.frame has columns \code{startFrame}, \code{startTime},
#'   \code{stopFrame}, \code{stopTime}, \code{duration}.
#' @export
leg_speed_phases <- function(trj,
                             swing_threshold,
                             method = c("central", "backward", "forward")) {
  method <- match.arg(method)

  # Swing = fasterThan threshold
  swing <- trajr::TrajSpeedIntervals(
    trj,
    fasterThan     = swing_threshold,
    slowerThan     = NULL,
    interpolateTimes = TRUE,
    diff           = method
  )

  # Stance = slowerThan or equal threshold
  stance <- trajr::TrajSpeedIntervals(
    trj,
    fasterThan     = NULL,
    slowerThan     = swing_threshold,
    interpolateTimes = TRUE,
    diff           = method
  )

  list(
    protraction = swing,
    retraction  = stance
  )
}


#' Extract swing/stance durations for all legs and trials
#'
#' @param trials List of objects as returned by process_trk_file().
#' @param swing_threshold Speed threshold for swing (same units as coordinates/time).
#'
#' @return data.frame with columns:
#'   file, instar, id, point, phase (protraction/retraction), duration.
#' @export
summarise_leg_speed_phases <- function(trials, swing_threshold) {
  results <- list()

  for (i in seq_along(trials)) {
    tr <- trials[[i]]
    meta <- tr$meta
    file <- meta$file

    for (pt in names(tr$coords)) {
      leg_trj <- tr$coords[[pt]]$trj

      # skip head/tail etc. if you only want legs:
      # if (!grepl("^l[0-9]+[tr]$", pt)) next

      phases <- leg_speed_phases(leg_trj, swing_threshold = swing_threshold)

      if (nrow(phases$protraction) > 0) {
        results[[length(results) + 1]] <- data.frame(
          file    = file,
          instar  = meta$instar,
          id      = meta$id,
          point   = pt,
          phase   = "protraction",
          duration = phases$protraction$duration
        )
      }

      if (nrow(phases$retraction) > 0) {
        results[[length(results) + 1]] <- data.frame(
          file    = file,
          instar  = meta$instar,
          id      = meta$id,
          point   = pt,
          phase   = "retraction",
          duration = phases$retraction$duration
        )
      }
    }
  }

  if (length(results) == 0) {
    return(data.frame(
      file = character(0),
      instar = character(0),
      id = character(0),
      point = character(0),
      phase = character(0),
      duration = numeric(0),
      stringsAsFactors = FALSE
    ))
  }

  do.call(rbind, results)
}

