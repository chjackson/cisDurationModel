#' A model for p_{ia} with double interval censored data
#'
#' @param prev_neg times of the last negative test prior to each episode
#' @param first_pos times of the first positive test of each episode
#' @param last_pos times of the last positive test of each episode
#' @param first_neg times of the first negative test following each episode
#' @export
pa_double_censor = function(prev_neg, first_pos, last_pos, first_neg) {
  stopifnot(all.equal(length(prev_neg), length(first_pos), length(last_pos), length(first_pos)))
  stopifnot(all(prev_neg < first_pos))
  stopifnot(all(first_pos <= last_pos))
  stopifnot(all(last_pos < first_neg))
  n_a = length(prev_neg)
  list(
    data_args = list(
      E = prev_neg,
      R = first_pos,
      L = last_pos,
      U = first_neg,
      n_a = n_a,
      max_S = max(first_neg - prev_neg) + 1
    ),
    data_code = list(
      "int<lower=0> n_a; // Number of episodes",
      "int<lower=0> max_S; // Number of days to compute survival for",
      "array[n_a] int E; // Day of previous negative",
      "array[n_a] int R; // Day of first positive",
      "array[n_a] int L; // Day after last positive",
      "array[n_a] int U; // Day of first negative"
    ),
    transformed_data_declare = list(
      "array[n_a] int<lower=1, upper=max_S-1> Rdash; // Day of first positive",
      "array[n_a] int<lower=1, upper=max_S-1> Ldash; // Day after last positive",
      "array[n_a] int<lower=2, upper=max_S> Udash; // Day of first negative"
    ),
    transformed_data_code = list(
      "for (i in 1:n_a) {",
      "\tRdash[i] = R[i] - E[i];",
      "\tLdash[i] = L[i] - E[i];",
      "\tUdash[i] = U[i] - E[i];",
      "}"
    ),
    transformed_parameters_declare = list(
      "array[n_a] real ll_numerator;"
    ),
    transformed_parameters_code = list(
      "for (i in 1:n_a) {",
      "\tll_numerator[i] = 0;",
      "\tfor (t in 1:Rdash[i]) {",
      "\t\tint max_positive = Udash[i] - t + 1;",
      "\t\tint min_positive = Ldash[i] - t + 1;",
      "\t\tll_numerator[i] += S[min_positive] - sensitivity * S[max_positive];",
      "\t}",
      "}"
    ),
    model = list("target += sum(log(ll_numerator));")
  )
}
