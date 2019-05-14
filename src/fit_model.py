"""Fit a model for the given data using the requested predictors and evaluate the model by time series cross-validation.
"""
import argparse
import json
import numpy as np
import pandas as pd
from scipy.optimize import minimize
import sys

from forecast.fitness_model import get_train_validate_timepoints
from forecast.metrics import add_pseudocounts_to_frequencies, negative_information_gain
from forecast.metrics import mean_absolute_error, sum_of_squared_errors, root_mean_square_error
from weighted_distances import get_distances_by_sample_names


def sum_of_differences(observed, estimated, y_diff, **kwargs):
    """
    Calculates the sum of squared errors for observed and estimated values.

    Parameters
    ----------
    observed : numpy.ndarray
        observed values

    estimated : numpy.ndarray
        estimated values

    y_diff : numpy.ndarray
        differences between observed and estimated values

    Returns
    -------
    float :
        sum of differences between estimated and observed future values
    """
    return np.sum(y_diff)


class ExponentialGrowthModel(object):
    def __init__(self, predictors, delta_time, l1_lambda, cost_function):
        """Construct an empty exponential growth model instance.

        Parameters
        ----------
        predictors : list
            a list of predictors to estimate coefficients for

        delta_time : float
            number of years into the future to project frequencies

        l1_lambda : float
            hyperparameter to scale L1 regularization penalty for non-zero coefficients

        cost_function : callable
            function returning the error to be minimized between observed and estimated values

        Returns
        -------
        ExponentialGrowthModel
        """
        self.predictors = predictors
        self.delta_time = delta_time
        self.l1_lambda = l1_lambda
        self.cost_function = cost_function

    def calculate_mean_stds(self, X, predictors):
        """Calculate mean standard deviations of predictors by timepoints prior to
        fitting.

        Parameters
        ----------
        X : pandas.DataFrame
            standardized tip attributes by timepoint

        predictors : ndarray
            predictor values per sample (n x p matrix for p predictors and n samples)

        Returns
        -------
        ndarray :
            mean standard deviation per predictor across all timepoints

        """
        return X.loc[:, ["timepoint"] + predictors].groupby("timepoint").std().mean().values

    def standardize_predictors(self, predictors, mean_stds):
        """Standardize the values for the given predictors by centering on the mean of
        each predictor and scaling by the mean standard deviation provided.

        Parameters
        ----------
        predictors : ndarray
            matrix of values per sample (rows) and predictor (columns)

        mean_stds : ndarray
            mean standard deviations of predictors across all training
            timepoints

        Returns
        -------
        ndarray :
            standardized predictor values

        """
        centered_predictors = (predictors - np.mean(predictors, axis=0))
        nonzero_stds = np.where(mean_stds)[0]

        if len(nonzero_stds) == 0:
            print(
                "Warning: all mean standard deviations are zero, so predictors will not be scaled",
                file=sys.stderr
            )
            return centered_predictors

        standardized_predictors = centered_predictors
        standardized_predictors[:, nonzero_stds] = standardized_predictors[:, nonzero_stds] / mean_stds[nonzero_stds]

        return standardized_predictors

    def get_fitnesses(self, coefficients, predictors):
        """Apply the coefficients to the predictors and sum them to get strain
        fitnesses.

        Parameters
        ----------
        coefficients : ndarray or list
            coefficients for given predictors

        predictors : ndarray
            predictor values per sample (n x p matrix for p predictors and n samples)

        Returns
        -------
        ndarray :
            fitnesses per sample
        """
        return np.sum(predictors * coefficients, axis=-1)

    def project_frequencies(self, initial_frequencies, fitnesses, delta_time):
        """Project the given initial frequencies into the future by the given delta time
        based on the given fitnesses.

        Returns the projected frequencies normalized to sum to 1.

        Parameters
        ----------
        initial_frequencies : ndarray
            floating point frequencies for all samples in a timepoint

        fitnesses : ndarray
            floating point fitnesses for all samples in same order as given frequencies

        delta_time : float
            number of years to project into the future

        Returns
        -------
        ndarray :
            projected and normalized frequencies
        """
        # Exponentiate the fitnesses and multiply them by strain frequencies.
        projected_frequencies = initial_frequencies * np.exp(fitnesses * self.delta_time)

        # Sum the projected frequencies.
        total_projected_frequencies = projected_frequencies.sum()

        # Normalize the projected frequencies.
        projected_frequencies = projected_frequencies / total_projected_frequencies

        return projected_frequencies

    def _fit(self, coefficients, X, y):
        """Calculate the error between observed and estimated values for the given
        parameters and data.

        Parameters
        ----------
        coefficients : ndarray
            coefficients for each of the model's predictors

        X : pandas.DataFrame
            standardized tip attributes by timepoint

        y : pandas.DataFrame
            final clade frequencies at delta time in the future from each
            timepoint in the given tip attributes table

        Returns
        -------
        float :
            error between estimated values using the given coefficients and
            input data and the observed values
        """
        # Estimate final frequencies.
        y_hat = self.predict(X, coefficients)

        # Merge estimated and observed frequencies. The left join enables
        # tracking of clades that die in the future and are therefore not
        # observed in the future frequencies data frame.
        frequencies = y_hat.merge(
            y,
            how="left",
            on=["timepoint", "clade_membership"],
            suffixes=["_estimated", "_observed"]
        )
        frequencies["frequency_observed"] = frequencies["frequency_observed"].fillna(0.0)

        # Calculate initial frequencies for use by cost function.
        initial_frequencies = X.groupby([
            "timepoint",
            "clade_membership"
        ])["frequency"].sum().reset_index()

        # Annotate future frequencies with initial frequencies.
        frequencies = frequencies.merge(
            initial_frequencies,
            how="inner",
            on=["timepoint", "clade_membership"]
        )

        # Calculate the error between the observed and estimated frequencies.
        error = self.cost_function(
            frequencies["frequency_observed"],
            frequencies["frequency_estimated"],
            initial=frequencies["frequency"]
        )
        l1_penalty = self.l1_lambda * np.abs(coefficients).sum()

        return error + l1_penalty

    def fit(self, X, y):
        """Fit a model to the given input data, producing beta coefficients for each of
        the model's predictors.

        Coefficients are stored in the `coef_` attribute, after the pattern of
        scikit-learn models.

        Parameters
        ----------
        X : pandas.DataFrame
            standardized tip attributes by timepoint

        y : pandas.DataFrame
            final clade frequencies at delta time in the future from each
            timepoint in the given tip attributes table

        Returns
        -------
        float :
            model training error

        """
        # Calculate mean standard deviations of predictors by timepoints prior
        # to fitting.
        self.mean_stds_ = self.calculate_mean_stds(X, self.predictors)

        # Find coefficients that minimize the model's cost function.
        initial_coefficients = np.random.normal(size=len(self.predictors))
        results = minimize(
            self._fit,
            initial_coefficients,
            args=(X, y),
            method="Nelder-Mead",
            options={"disp": True}
        )
        self.coef_ = results.x

        training_error = self.score(X, y)

        return training_error

    def predict(self, X, coefficients=None, mean_stds=None):
        """Calculate the estimate final frequencies of all clades in the given tip
        attributes data frame using previously calculated beta coefficients.

        Parameters
        ----------
        X : pandas.DataFrame
            standardized tip attributes by timepoint

        coefficients : ndarray
            optional coefficients to use for each of the model's predictors
            instead of the model's currently defined coefficients

        mean_stds : ndarray
            optional mean standard deviations of predictors across all training
            timepoints

        Returns
        -------
        pandas.DataFrame
            estimated final clade frequencies at delta time in the future for
            each clade from each timepoint in the given tip attributes table

        """
        # Use model coefficients, if none are provided.
        if coefficients is None:
            coefficients = self.coef_

        if mean_stds is None:
            mean_stds = self.mean_stds_

        estimated_frequencies = []
        for timepoint, timepoint_df in X.groupby("timepoint"):
            # Select predictors from the timepoint.
            predictors = timepoint_df.loc[:, self.predictors].values

            # Standardize predictors by timepoint centering by means at
            # timepoint and mean standard deviation provided.
            standardized_predictors = self.standardize_predictors(predictors, mean_stds)

            # Select frequencies from timepoint.
            initial_frequencies = timepoint_df["frequency"].values

            # Calculate fitnesses.
            fitnesses = self.get_fitnesses(coefficients, standardized_predictors)

            # Project frequencies.
            projected_frequencies = self.project_frequencies(
                initial_frequencies,
                fitnesses,
                self.delta_time
            )

            # Sum the estimated frequencies by clade.
            projected_timepoint_df = timepoint_df[["timepoint", "clade_membership"]].copy()
            projected_timepoint_df["frequency"] = projected_frequencies
            projected_clade_frequencies = projected_timepoint_df.groupby([
                "timepoint",
                "clade_membership"
            ])["frequency"].sum().reset_index()

            estimated_frequencies.append(projected_clade_frequencies)

        # Collect all estimated frequencies by timepoint.
        estimated_frequencies = pd.concat(estimated_frequencies)
        return estimated_frequencies

    def score(self, X, y):
        """Calculate model error between the estimated final clade frequencies for the
        given tip attributes, `X`, and the observed final clade frequencies in
        `y`.

        Parameters
        ----------
        X : pandas.DataFrame
            standardized tip attributes by timepoint

        y : pandas.DataFrame
            final clade frequencies at delta time in the future from each
            timepoint in the given tip attributes table

        Returns
        -------
        float :
            model error
        """
        return self._fit(self.coef_, X, y)


class DistanceExponentialGrowthModel(ExponentialGrowthModel):
    def __init__(self, predictors, delta_time, l1_lambda, cost_function, distances):
        super().__init__(predictors, delta_time, l1_lambda, cost_function)
        self.distances = distances

    def _fit(self, coefficients, X, y):
        """Calculate the error between observed and estimated values for the given
        parameters and data.

        Parameters
        ----------
        coefficients : ndarray
            coefficients for each of the model's predictors

        X : pandas.DataFrame
            standardized tip attributes by timepoint

        y : pandas.DataFrame
            final weighted distances at delta time in the future from each
            timepoint in the given tip attributes table

        Returns
        -------
        float :
            error between estimated values using the given coefficients and
            input data and the observed values
        """
        # Estimate target values.
        y_hat = self.predict(X, coefficients)

        # Merge estimated and observed target values.
        targets = y_hat.merge(
            y,
            how="inner",
            on=["timepoint", "strain"],
            suffixes=["_estimated", "_observed"]
        )

        # Calculate the error between the observed and estimated targets.
        error = self.cost_function(
            targets["y_observed"],
            targets["y_estimated"],
            y_diff=targets["y_diff"]
        )
        l1_penalty = self.l1_lambda * np.abs(coefficients).sum()

        return error + l1_penalty

    def predict(self, X, coefficients=None):
        """Calculate the estimated final weighted distance between tips at each
        timepoint and at that timepoint plus delta months in the future.

        Parameters
        ----------
        X : pandas.DataFrame
            standardized tip attributes by timepoint

        coefficients : ndarray
            optional coefficients to use for each of the model's predictors
            instead of the model's currently defined coefficients

        Returns
        -------
        pandas.DataFrame
            estimated weighted distances at delta time in the future for
            each tip from each timepoint in the given tip attributes table

        """
        # Use model coefficients, if none are provided.
        if coefficients is None:
            coefficients = self.coef_
            model_is_fit = True
        else:
            model_is_fit = False

        estimated_targets = []
        for timepoint, timepoint_df in X.groupby("timepoint"):
            # Select predictors from the timepoint.
            predictors = timepoint_df.loc[:, self.predictors].values

            # Select frequencies from timepoint.
            initial_frequencies = timepoint_df["frequency"].values

            # Calculate fitnesses.
            fitnesses = self.get_fitnesses(coefficients, predictors)

            # Project frequencies.
            projected_frequencies = self.project_frequencies(
                initial_frequencies,
                fitnesses,
                self.delta_time
            )

            # Calculate observed distance between current tips and the future
            # using projected frequencies and weighted distances to the future.
            projected_timepoint_df = timepoint_df[["timepoint", "strain", "frequency", "weighted_distance_to_present", "weighted_distance_to_future"]].copy()
            projected_timepoint_df["projected_frequency"] = projected_frequencies
            projected_timepoint_df["y_diff"] = projected_timepoint_df["projected_frequency"] * projected_timepoint_df["weighted_distance_to_future"]

            if model_is_fit or self.cost_function != sum_of_differences:
                # Calculate estimate distance between current tips and future tips
                # based on projections of current tips.
                estimated_weighted_distance_to_future = []
                for current_tip, current_tip_frequency in projected_timepoint_df.loc[:, ["strain", "frequency"]].values:
                    weighted_distance_to_future = 0.0
                    for other_tip, other_tip_projected_frequency in projected_timepoint_df.loc[:, ["strain", "projected_frequency"]].values:
                        weighted_distance_to_future += other_tip_projected_frequency * self.distances[current_tip][other_tip]

                    estimated_weighted_distance_to_future.append(
                        current_tip_frequency * weighted_distance_to_future
                    )

                projected_timepoint_df["y"] = np.array(estimated_weighted_distance_to_future)
            else:
                projected_timepoint_df["y"] = np.nan

            estimated_targets.append(projected_timepoint_df)

        # Collect all estimated targets by timepoint.
        estimated_targets = pd.concat(estimated_targets, ignore_index=True)
        return estimated_targets


def cross_validate(model, data, targets, train_validate_timepoints, coefficients=None, group_by="clade_membership"):
    """Calculate cross-validation scores for the given data and targets across the
    given train/validate timepoints.

    Parameters
    ----------
    model : ExponentialGrowthModel
        an instance of a model with defined hyperparameters including a list of
        predictors to use for fitting

    data : pandas.DataFrame
        standardized input attributes to use for model fitting

    targets : pandas.DataFrame
        observed outputs to fit the model to

    train_validate_timepoints : list
        a list of dictionaries of lists indexed by "train" and "validate" keys
        and containing timepoints to use for model training and validation,
        respectively

    coefficients : ndarray
        an optional array of fixed coefficients for the given model's predictors
        to use when calculating cross-validation error for specific models
        (e.g., naive forecasts)

    Returns
    -------
    list
        a list of dictionaries containing cross-validation results with scores,
        training and validation results, and beta coefficients per timepoint

    """
    results = []

    for timepoints in train_validate_timepoints:
        # Get training and validation timepoints.
        training_timepoints = pd.to_datetime(timepoints["train"])
        validation_timepoint = pd.to_datetime(timepoints["validate"])

        # Get training data by timepoints.
        training_X = data[data["timepoint"].isin(training_timepoints)].copy()
        training_y = targets[targets["timepoint"].isin(training_timepoints)].copy()

        # Fit a model to the training data.
        if coefficients is None:
            training_error = model.fit(training_X, training_y)
        else:
            model.coef_ = coefficients
            model.mean_stds_ = model.calculate_mean_stds(training_X, model.predictors)
            training_error = model.score(training_X, training_y)

        # Get validation data by timepoints.
        validation_X = data[data["timepoint"] == validation_timepoint].copy()
        validation_y = targets[targets["timepoint"] == validation_timepoint].copy()

        # Calculate the model score for the validation data.
        validation_error = model.score(validation_X, validation_y)

        # Get the estimated frequencies for training and validation sets to export.
        training_y_hat = model.predict(training_X)
        validation_y_hat = model.predict(validation_X)

        # Convert timestamps to a serializable format.
        training_X["timepoint"] = training_X["timepoint"].dt.strftime("%Y-%m-%d")
        training_y["timepoint"] = training_y["timepoint"].dt.strftime("%Y-%m-%d")
        training_y_hat["timepoint"] = training_y_hat["timepoint"].dt.strftime("%Y-%m-%d")
        validation_X["timepoint"] = validation_X["timepoint"].dt.strftime("%Y-%m-%d")
        validation_y["timepoint"] = validation_y["timepoint"].dt.strftime("%Y-%m-%d")
        validation_y_hat["timepoint"] = validation_y_hat["timepoint"].dt.strftime("%Y-%m-%d")

        # Store training results, beta coefficients, and validation results.
        results.append({
            "predictors": model.predictors,
            "training_data": {
                "X": training_X.to_dict(orient="records"),
                "y": training_y.to_dict(orient="records"),
                "y_hat": training_y_hat.to_dict(orient="records")
            },
            "training_n": training_X[group_by].unique().shape[0],
            "training_error": training_error,
            "coefficients": model.coef_.tolist(),
            "mean_stds": model.mean_stds_.tolist(),
            "validation_data": {
                "X": validation_X.to_dict(orient="records"),
                "y": validation_y.to_dict(orient="records"),
                "y_hat": validation_y_hat.to_dict(orient="records")
            },
            "validation_n": validation_X[group_by].unique().shape[0],
            "validation_error": validation_error,
            "last_training_timepoint": training_timepoints[-1].strftime("%Y-%m-%d"),
            "validation_timepoint": validation_timepoint.strftime("%Y-%m-%d")
        })

    # Return results for all validation timepoints.
    return results


def summarize_cross_validation_scores(scores):
    """Summarize model errors including in-sample errors by AIC, out-of-sample
    errors by cross-validation, and beta parameters across timepoints.

    Parameters
    ----------
    scores : list
        a list of cross-validation results including training errors,
        cross-validation errors, and beta coefficients

    Returns
    -------
    dict :
        a dictionary of all cross-validation results plus summary statistics for
        training, cross-validation, and beta coefficients
    """
    summary = {
        "scores": scores,
        "predictors": scores[0]["predictors"]
    }

    validation_errors = [score["validation_error"] for score in scores]
    summary["cv_error_mean"] = np.mean(validation_errors)
    summary["cv_error_std"] = np.std(validation_errors)

    coefficients = np.array([
        np.array(score["coefficients"])
        for score in scores
    ])
    summary["coefficients_mean"] = coefficients.mean(axis=0).tolist()
    summary["coefficients_std"] = coefficients.std(axis=0).tolist()

    mean_stds = np.array([
        np.array(score["mean_stds"])
        for score in scores
    ])
    summary["mean_stds_mean"] = mean_stds.mean(axis=0).tolist()
    summary["mean_stds_std"] = mean_stds.std(axis=0).tolist()

    return summary


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--tip-attributes", required=True, help="tab-delimited file describing tip attributes at all timepoints with standardized predictors")
    parser.add_argument("--output", required=True, help="JSON representing the model fit with training and cross-validation results, beta coefficients for predictors, and summary statistics")
    parser.add_argument("--predictors", required=True, nargs="+", help="tip attribute columns to use as predictors of final clade frequencies")
    parser.add_argument("--delta-months", required=True, type=int, help="number of months to project clade frequencies into the future")
    parser.add_argument("--target", required=True, choices=["clades", "distances"], help="target for models to fit")
    parser.add_argument("--final-clade-frequencies", help="tab-delimited file of clades per timepoint and their corresponding tips and tip frequencies at the given delta time in the future")
    parser.add_argument("--distances", help="tab-delimited file of distances between pairs of samples")
    parser.add_argument("--training-window", type=int, default=4, help="number of years required for model training")
    parser.add_argument("--l1-lambda", type=float, default=0.0, help="L1 regularization lambda")
    parser.add_argument("--cost-function", default="sse", choices=["sse", "rmse", "mae", "information_gain", "diffsum"], help="name of the function that returns the error between observed and estimated values")
    parser.add_argument("--pseudocount", type=float, help="pseudocount numerator to adjust all frequencies by, enabling some information theoretic metrics like information gain")

    args = parser.parse_args()

    # Load standardized tip attributes subsetting to tip name, clade, frequency,
    # and requested predictors.
    tips = pd.read_csv(
        args.tip_attributes,
        sep="\t",
        parse_dates=["timepoint"]
    )

    if args.target == "clades":
        # Load final clade tip frequencies.
        final_clade_tip_frequencies = pd.read_csv(
            args.final_clade_frequencies,
            sep="\t",
            parse_dates=["initial_timepoint", "final_timepoint"]
        )

        # If a pseudocount numerator has been provided, update the given tip
        # frequencies both from current and future timepoints.
        if args.pseudocount is not None and args.pseudocount > 0.0:
            tips = add_pseudocounts_to_frequencies(tips, args.pseudocount)
            print("Sum of tip frequencies by timepoint: ",
                  tips.groupby("timepoint")["frequency"].sum())
            final_clade_tip_frequencies = add_pseudocounts_to_frequencies(
                final_clade_tip_frequencies,
                args.pseudocount,
                timepoint_column="initial_timepoint"
            )
            print("Sum of tip frequencies by timepoint: ",
                  final_clade_tip_frequencies.groupby("initial_timepoint")["frequency"].sum())

        # Aggregate final clade frequencies.
        final_clade_frequencies = final_clade_tip_frequencies.groupby([
            "initial_timepoint",
            "clade_membership"
        ])["frequency"].sum().reset_index()

        # Rename initial timepoint column for comparison with tip attribute data.
        targets = final_clade_frequencies.rename(
            columns={"initial_timepoint": "timepoint"}
        )
        model_class = ExponentialGrowthModel
        model_kwargs = {}
        group_by_attribute = "clade_membership"
    elif args.target == "distances":
        # Scale each tip's weighted distance to future populations by the tip's
        # current frequency.
        tips["y"] = tips["frequency"] * tips["weighted_distance_to_future"]
        targets = tips.loc[:, ["strain", "timepoint", "y"]].copy()
        model_class = DistanceExponentialGrowthModel
        distances = pd.read_csv(args.distances, sep="\t")
        distances_by_sample_names = get_distances_by_sample_names(distances)
        model_kwargs = {"distances": distances_by_sample_names}
        group_by_attribute = "strain"

    # Identify all available timepoints from tip attributes.
    timepoints = tips["timepoint"].dt.strftime("%Y-%m-%d").unique()

    # Identify train/validate splits from timepoints.
    train_validate_timepoints = get_train_validate_timepoints(
        timepoints,
        args.delta_months,
        args.training_window
    )

    # Select the cost function.
    if args.cost_function == "sse":
        cost_function = sum_of_squared_errors
    elif args.cost_function == "rmse":
        cost_function = root_mean_square_error
    elif args.cost_function == "mae":
        cost_function = mean_absolute_error
    elif args.cost_function == "information_gain":
        cost_function = negative_information_gain
    elif args.cost_function == "diffsum":
        cost_function = sum_of_differences

    # For each train/validate split, fit a model to the training data, and
    # evaluate the model with the validation data, storing the training results,
    # beta parameters, and validation results.
    delta_time = args.delta_months / 12.0
    model = model_class(
        predictors=args.predictors,
        delta_time=delta_time,
        l1_lambda=args.l1_lambda,
        cost_function=cost_function,
        **model_kwargs
    )

    # If this is a naive model, set the coefficients to zero so cross-validation
    # can run under naive model conditions.
    if "naive" in args.predictors:
        coefficients = np.zeros(len(args.predictors))
    else:
        coefficients = None

    scores = cross_validate(
        model,
        tips,
        targets,
        train_validate_timepoints,
        coefficients,
        group_by=group_by_attribute
    )

    # Summarize model errors including in-sample errors by AIC, out-of-sample
    # errors by cross-validation, and beta parameters across timepoints.
    model_results = summarize_cross_validation_scores(scores)

    # Annotate parameters used to produce models.
    model_results["cost_function"] = args.cost_function
    model_results["l1_lambda"] = args.l1_lambda
    model_results["delta_months"] = args.delta_months
    model_results["training_window"] = args.training_window
    model_results["pseudocount"] = args.pseudocount

    # Save model fitting hyperparameters, raw results, and summary of results to
    # JSON.
    with open(args.output, "w") as fh:
        json.dump(model_results, fh, indent=1)
