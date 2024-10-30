test_that("run_cosinor works", {
  period <- 24
  n_reps <- 1
  time <- rep(seq(0, 24, by = 2), each = n_reps)
  acrophases <- c(1, 10, 21)
  y <- lapply(acrophases, function(ac) make_data(time = time, acrophase = ac, error_sd = .1)$value)
  x <- matrix(unlist(y), byrow = TRUE, ncol = length(time))
  expect_true(is(run_cosinor(x = x, time = time), "data.frame"))

  xE <- new("EList", list(E = x))
  expect_true(is(run_cosinor(x = xE, time = time), "data.frame"))

  xE$weights <- rnorm(nrow(xE), 1, .01)
  expect_true(is(run_cosinor(x = xE, time = time), "data.frame"))

  xNA <- x
  xNA[1, 1] <- NA
  expect_error(run_cosinor(list(1), time = time))
  expect_error(run_cosinor(x = xNA, time = time))
  expect_error(run_cosinor(x = xNA, time = time[1:2]))
  expect_error(run_cosinor(x = xNA, time = as.character(time)))
  expect_error(run_cosinor(x = x, time = time, period = "24"))
  expect_error(run_cosinor(x = xNA, time = time, eBayes_args = "1"))
  expect_error(run_cosinor(x = xNA, time = time, lmFit_args = "1"))
})

test_that("make_design works", {
  expect_true(is(make_design(time = c(1, 2, 3), period = 24), "matrix"))
  expect_error(make_design(time = c(1, 2, "3")))
  expect_error(make_design(time = c(1, 2, 3), period = "24"))
  expect_error(make_design(time = c(1, 2, 3), period = -1))
})

test_that("circular_density works", {
  data <- abs(c(rnorm(1000, 12, 5), rnorm(500, 23))) %% 24
  expect_true(is(circular_density(data), "data.frame"))
  expect_error(circular_density(data, to = -1))
  expect_error(circular_density(LETTERS[1:10], from = -1))
  expect_error(circular_density(data, to = -1))
  expect_error(circular_density(data, to = NA))
  expect_error(circular_density(data, density_args = 1))
})

test_that("rad2period works", {
  expect_type(rad2period(pi), "double")
  expect_error(rad2period(2 * pi + 1))
  expect_error(rad2period(1, period = -1))
})

test_that("circular_distance works", {
  expect_type(circular_distance(1, 2), "double")
  expect_error(circular_distance(1, NA))
  expect_error(circular_distance(1, 1, NA))
  expect_error(circular_distance(NA, 1, 2))
  expect_error(circular_distance(1, 2, -3))
})

test_that("circular_distance_reference works", {
  expect_type(circular_distance_reference(1, c(1, 2, 3), period = 24), "double")
  expect_error(circular_distance_reference(1, 1, NA))
  expect_error(circular_distance_reference(-1, 1, 24))
  expect_error(circular_distance_reference(1, -1, 24))
})

test_that("circular_intersect works", {
  expect_type(circular_intersect(1, c(1, 2)), "logical")
  expect_error(circular_intersect(-1, c(1, 2)))
  expect_error(circular_intersect(1, c(1, -2)))
  expect_error(circular_intersect(1, c(-1, 2)))
  expect_error(circular_intersect(1, c(1, 1)))
})

test_that("run_model_selection works", {
  y <- rbind(
    # amplitude change
    make_data(time = seq(1, 21, 4), amplitude = 2),
    make_data(time = seq(1, 21, 4), amplitude = 1),
    # acrophase change
    make_data(time = seq(1, 21, 4), acrophase = 1),
    make_data(time = seq(1, 21, 4), acrophase = 13),
    # same rhythmicity
    make_data(time = seq(1, 21, 4), acrophase = 1),
    make_data(time = seq(1, 21, 4), acrophase = 1),
    # loss of rhythmicity
    make_data(time = seq(1, 21, 4), amplitude = 1),
    make_data(time = seq(1, 21, 4), amplitude = 0),
    # gain of rhythmicity
    make_data(time = seq(1, 21, 4), amplitude = 0),
    make_data(time = seq(1, 21, 4), amplitude = 1)
  )
  x <- matrix(y$value, byrow = TRUE, ncol = 12)
  rownames(x) <- paste0("gene_", 1:nrow(x))
  time <- y$time[1:12]
  condition <- rep(c("A", "B"), each = 6)
  r <- run_model_selection(x = x, condition = condition, time = time)
  expect_true(is(r, "data.frame"))
  expect_error(run_model_selection(x = x, condition = condition, time = time[1:2]))
  expect_error(run_model_selection(x = x, condition = condition[1:2], time = time))
  expect_error(run_model_selection(x = x, condition = condition, time = time, period = -1))
  expect_error(run_model_selection(x = x, condition = condition, time = time, sample_weights = FALSE))
})
