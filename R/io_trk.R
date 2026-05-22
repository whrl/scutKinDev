#' List Tracker .trk files
#'
#' Recursively finds all Tracker .trk files under a directory.
#'
#' @param path Root directory to search.
#' @param pattern File name pattern, default "\*.trk$".
#' @return Character vector of full file paths.
#' @export
discover_trk_files <- function(path = ".", pattern = "*.trk$") {
  list.files(
    path,
    pattern   = pattern,
    full.names = TRUE,
    recursive  = TRUE,
    ignore.case = TRUE
  )
}

#' Extract Tracker properties from a .trk XML file
#'
#' Extracts selected `<property>` values from a Tracker .trk XML file.
#'
#' @param xml_file Path to a .trk (XML) file.
#' @param prop_names Character vector of property names to extract.
#'
#' @return A one-row data.frame with columns named by `prop_names`.
#' @importFrom XML xmlParse getNodeSet xmlValue xmlGetAttr
#' @export
extract_props <- function(xml_file, prop_names) {
  doc   <- XML::xmlParse(xml_file)
  values <- setNames(vector("list", length(prop_names)), prop_names)
  types  <- setNames(vector("character", length(prop_names)), prop_names)

  for (prop in prop_names) {
    node_set <- XML::getNodeSet(doc, paste0('//property[@name="', prop, '"]'))
    node <- if (length(node_set) > 0) node_set[[1]] else NULL
    values[[prop]] <- if (!is.null(node)) XML::xmlValue(node) else NA_character_
    types[[prop]]  <- if (!is.null(node)) {
      XML::xmlGetAttr(node, "type", NA_character_)
    } else {
      NA_character_
    }
  }

  conversion_map <- list(
    int    = as.integer,
    double = as.double,
    string = as.character
  )

  for (prop in prop_names) {
    typ <- types[[prop]]
    val <- values[[prop]]
    if (!is.na(typ) && typ %in% names(conversion_map)) {
      values[[prop]] <- conversion_map[[typ]](val)
    } else {
      values[[prop]] <- val
    }
  }

  result <- as.data.frame(values, stringsAsFactors = FALSE)
  attr(result, "property_types") <- types
  result
}

#' Extract point-mass marker names from a .trk file
#'
#' Reads a Tracker .trk file and returns the names of all PointMass objects.
#'
#' @param xml_file Path to a .trk (XML) file.
#'
#' @return A character vector of unique marker names.
#' @importFrom XML xmlParse getNodeSet xmlValue
#' @export
extract_pointmass_names <- function(xml_file) {
  doc      <- XML::xmlParse(xml_file)
  pm_nodes <- XML::getNodeSet(
    doc,
    "//object[@class='org.opensourcephysics.cabrillo.tracker.PointMass']"
  )
  names <- vapply(pm_nodes, function(node) {
    name_node <- XML::getNodeSet(node, "property[@name='name']")
    if (length(name_node) > 0) XML::xmlValue(name_node[[1]]) else NA_character_
  }, FUN.VALUE = character(1))
  unique(stats::na.omit(names))
}

#' Extract Tracker point-mass coordinates
#'
#' Extracts frame-wise x/y coordinates for a single point mass from a .trk file.
#'
#' @param filename Path to a .trk (XML) file.
#' @param pointname Name of the point mass (as used in Tracker).
#'
#' @return A data.frame with columns `frame`, `x`, `y`.
#' @importFrom XML xmlParse getNodeSet xmlValue xmlGetAttr xmlParent
#' @export
extract_pointmass_coords <- function(filename, pointname) {
  doc <- XML::xmlParse(filename)

  node_pattern <- paste0(
    "/child::object//child::property[@name='name' and text()='",
    pointname,
    "']"
  )
  nodeset <- XML::getNodeSet(doc, node_pattern)
  if (length(nodeset) == 0) {
    stop("PointMass '", pointname, "' not found in file: ", filename)
  }
  parent <- XML::xmlParent(nodeset[[1]])

  x_nodes <- XML::getNodeSet(parent, "property//object//property[@name='x']")
  y_nodes <- XML::getNodeSet(parent, "property//object//property[@name='y']")
  f_nodes <- XML::getNodeSet(parent, "property/property[@type='object']")

  x <- as.numeric(vapply(x_nodes, XML::xmlValue, FUN.VALUE = character(1)))
  y <- as.numeric(vapply(y_nodes, XML::xmlValue, FUN.VALUE = character(1)))
  frames <- as.numeric(gsub("\\]|\\[", "", vapply(f_nodes, XML::xmlGetAttr,
                                                  FUN.VALUE = character(1),
                                                  "name")))

  df <- data.frame(frame = frames, x = x, y = y)
  df[order(df$frame), , drop = FALSE]
}

#' Parse metadata from Tracker .trk files
#'
#' Extracts a standard set of global properties from multiple .trk files.
#'
#' @param files Character vector of .trk file paths.
#'
#' @return A data.table with one row per file and standard meta columns.
#' @import data.table
#' @export
parse_trk_meta <- function(files) {
  props <- c(
    "semantic_version",
    "width",
    "height",
    "center_x",
    "center_y",
    "duration",
    "frame_count",
    "frame_rate",
    "video_framecount",
    "startframe",
    "stepsize",
    "stepcount",
    "starttime",
    "readout",
    "rate",
    "delta_t",
    "frame",
    "fixedorigin",
    "fixedangle",
    "fixedscale",
    "locked",
    "xorigin",
    "yorigin",
    "angle",
    "xscale",
    "yscale",
    "length_unit"
  )

  meta <- data.table::rbindlist(
    lapply(files, function(f) {
      m <- extract_props(f, props)
      m$filename <- basename(f)
      m$file     <- f
      m
    }),
    fill = TRUE
  )

  meta$delta_t_s     <- meta$delta_t / 1000
  meta$frame_rate_hz <- 1 / meta$delta_t_s
  meta
}

#' Process a single Tracker .trk file
#'
#' Reads one .trk file, extracts meta information and trajectories for all
#' point masses, converts to physical units, and computes basic speed summaries.
#'
#' @param file Path to a .trk file.
#'
#' @return A list with components:
#'   \describe{
#'     \item{meta}{data.frame with `file`, `instar`, `mass`, `nr`, `id`,
#'                 `xscale`, `yscale`, `delta_t`, `length_unit`.}
#'     \item{coords}{Named list; one element per marker. Each element is a list
#'                   with components `trj` (trajr object),
#'                   `derivs` (data.frame from \code{TrajDerivatives}),
#'                   `summary` (one-row data.frame with speed summaries).}
#'   }
#'
#' @importFrom trajr TrajFromCoords TrajSmoothSG TrajDerivatives
#' @export
process_trk_file <- function(file) {
  instar <- regmatches(file, gregexpr("instar_[0-9]+", file))[[1]][1]
  mass   <- as.numeric(regmatches(file, regexec("([0-9]+\\.?[0-9]*)mg", file))[[1]][2])
  nr     <- as.numeric(regmatches(file, regexec("mg_([0-9]+)", file))[[1]][2])
  id     <- regmatches(file, gregexpr("(vie|ibzi)[0-9]+", file))[[1]][1]

  markers <- extract_pointmass_names(file)

  meta <- extract_props(
    file,
    c("xscale", "yscale", "delta_t", "length_unit")
  )
  meta$file        <- file
  meta$instar      <- instar
  meta$mass        <- mass
  meta$nr          <- nr
  meta$id          <- id

  coords_list <- setNames(
    lapply(markers, function(pt) {
      coords <- extract_pointmass_coords(file, pt)
      coords$time <- coords$frame * meta$delta_t / 1000
      coords$x    <- coords$x / meta$xscale
      coords$y    <- coords$y / meta$yscale

      trj <- trajr::TrajFromCoords(
        coords,
        xCol        = "x",
        yCol        = "y",
        timeCol      = "time",
        spatialUnits = meta$length_unit
      )

      derivs <- trajr::TrajDerivatives(
        trajr::TrajSmoothSG(trj, p = 3, n = 5)
      )

      summary <- data.frame(
        speed_median = stats::median(derivs$speed, na.rm = TRUE),
        speed_mean   = mean(derivs$speed, na.rm = TRUE),
        speed_sd     = stats::sd(derivs$speed, na.rm = TRUE),
        speed_max    = max(derivs$speed, na.rm = TRUE)
      )

      list(trj = trj, derivs = derivs, summary = summary)
    }),
    markers
  )

  list(meta = meta, coords = coords_list)
}
