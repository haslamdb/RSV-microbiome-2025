#!/usr/bin/env python
# scripts/kraken/feature_selection.py
import os
import sys
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.ensemble import RandomForestRegressor
from sklearn.preprocessing import StandardScaler, LabelEncoder
from sklearn.model_selection import train_test_split
from sklearn.metrics import mean_squared_error, r2_score, mean_absolute_error
import scipy.spatial.distance as ssd
import logging
import traceback
import argparse
from pathlib import Path

# Setup simple logger for standalone use
logger = logging.getLogger('feature_selection')
logger.setLevel(logging.INFO)
if not logger.handlers:
    handler = logging.StreamHandler()
    handler.setFormatter(logging.Formatter('%(asctime)s - %(levelname)s - %(message)s'))
    logger.addHandler(handler)

def log_print(message, level="info"):
    """Log a message at specified level"""
    if level.lower() == "debug":
        logger.debug(message)
    elif level.lower() == "info":
        logger.info(message)
    elif level.lower() == "warning":
        logger.warning(message)
    elif level.lower() == "error":
        logger.error(message)
    elif level.lower() == "critical":
        logger.critical(message)
    else:
        logger.info(message)

def check_file_exists_with_logger(file_path, desc, logger_func=log_print):
    """Check if a file exists and log an error if it doesn't"""
    if not os.path.exists(file_path):
        logger_func(f"{desc} not found: {file_path}", level="error")
        return False
    return True

def transform_abundance_data(abundance_df, transform="clr", logger_func=log_print):
    """
    Apply transformation to abundance data
    
    Args:
        abundance_df: DataFrame with abundance data
        transform: Transformation to apply (clr, hellinger, log, none)
        logger_func: Function for logging
        
    Returns:
        Transformed abundance DataFrame
    """
    if transform.lower() == "clr":
        logger_func("Applying CLR transformation", level="info")
        # Add small pseudocount to zeros
        df_pseudo = abundance_df.replace(0, np.nextafter(0, 1))
        # Log transform
        df_log = np.log(df_pseudo)
        # Subtract column-wise mean (CLR transformation)
        abundance_transformed = df_log.subtract(df_log.mean(axis=0), axis=1)
    elif transform.lower() == "hellinger":
        logger_func("Applying Hellinger transformation", level="info")
        # Calculate sample totals
        sample_totals = abundance_df.sum(axis=0)
        # Divide each value by its sample total
        df_rel = abundance_df.div(sample_totals, axis=1)
        # Apply square root
        abundance_transformed = np.sqrt(df_rel)
    elif transform.lower() == "log":
        logger_func("Applying log transformation", level="info")
        abundance_transformed = np.log1p(abundance_df)
    else:
        logger_func("No transformation applied", level="info")
        abundance_transformed = abundance_df.copy()
    
    return abundance_transformed

def encode_categorical_variables(metadata_df, categorical_features, logger_func=log_print):
    """
    Encode categorical variables for Random Forest
    
    Args:
        metadata_df: DataFrame with metadata
        categorical_features: List of categorical variables to encode
        logger_func: Function for logging
        
    Returns:
        DataFrame with encoded variables
    """
    encoded_df = metadata_df.copy()
    
    for col in categorical_features:
        if col in encoded_df.columns:
            # Fill NAs
            encoded_df[col] = encoded_df[col].fillna('missing')
            
            try:
                # Use label encoding for Random Forest
                encoded_df[col] = LabelEncoder().fit_transform(encoded_df[col])
                logger_func(f"Encoded categorical variable: {col}", level="info")
            except Exception as e:
                logger_func(f"Error encoding {col}, dropping column: {str(e)}", level="warning")
                encoded_df = encoded_df.drop(columns=[col])
    
    return encoded_df

def scale_numerical_variables(metadata_df, numerical_features, logger_func=log_print):
    """
    Scale numerical variables for Random Forest
    
    Args:
        metadata_df: DataFrame with metadata
        numerical_features: List of numerical variables to scale
        logger_func: Function for logging
        
    Returns:
        DataFrame with scaled variables
    """
    scaled_df = metadata_df.copy()
    
    # Create scaler
    scaler = StandardScaler()
    
    # Get only numeric columns that exist in the DataFrame
    numeric_cols = [col for col in numerical_features if col in scaled_df.columns 
                   and np.issubdtype(scaled_df[col].dtype, np.number)]
    
    if numeric_cols:
        # Scale numeric columns
        scaled_df[numeric_cols] = scaler.fit_transform(scaled_df[numeric_cols])
        logger_func(f"Scaled numerical variables: {numeric_cols}", level="info")
    
    return scaled_df

def create_pairwise_differences(X_scaled, logger_func=log_print):
    """
    Create matrix of pairwise differences between samples
    
    Args:
        X_scaled: Scaled/encoded metadata DataFrame
        logger_func: Function for logging
        
    Returns:
        Tuple of (pairwise differences array, indices)
    """
    logger_func("Creating pairwise differences matrix", level="info")
    
    # Convert to numpy array
    X_array = X_scaled.values
    
    # Create pairwise differences
    X_diff = np.abs(X_array[:, None, :] - X_array[None, :, :])
    
    # Get upper triangle indices
    n_samples = X_array.shape[0]
    indices = np.triu_indices(n_samples, k=1)
    
    # Extract upper triangle (pairwise differences)
    X_diff_upper = X_diff[indices]
    
    logger_func(f"Created {X_diff_upper.shape[0]} pairwise differences with {X_diff_upper.shape[1]} features", level="info")
    
    return X_diff_upper, indices

def calculate_microbiome_distances(abundance_df, distance_metric="bray", logger_func=log_print):
    """
    Calculate pairwise distances between microbiome samples
    
    Args:
        abundance_df: DataFrame with abundance data
        distance_metric: Distance metric to use
        logger_func: Function for logging
        
    Returns:
        Array of pairwise distances
    """
    logger_func(f"Calculating microbiome distances using {distance_metric} metric", level="info")
    
    try:
        if distance_metric.lower() == "bray":
            distances = ssd.pdist(abundance_df.values, metric="braycurtis")
        elif distance_metric.lower() == "jaccard":
            distances = ssd.pdist(abundance_df.values, metric="jaccard")
        elif distance_metric.lower() == "euclidean":
            distances = ssd.pdist(abundance_df.values, metric="euclidean")
        else:
            logger_func(f"Unknown distance metric: {distance_metric}. Using Bray-Curtis.", level="warning")
            distances = ssd.pdist(abundance_df.values, metric="braycurtis")
        
        logger_func(f"Calculated {len(distances)} pairwise distances", level="info")
        return distances
    
    except Exception as e:
        logger_func(f"Error calculating distances: {str(e)}", level="error")
        logger_func(traceback.format_exc(), level="error")
        return None

def train_random_forest(X_diff, y_distances, feature_names, test_size=0.2, 
                       n_estimators=100, max_features="sqrt", random_state=42, logger_func=log_print):
    """
    Train Random Forest model to predict microbiome distances from metadata differences
    
    Args:
        X_diff: Pairwise differences in metadata
        y_distances: Pairwise microbiome distances
        feature_names: Names of features
        test_size: Proportion of data for testing
        n_estimators: Number of trees in the forest
        max_features: Max features to consider at each split
        random_state: Random seed for reproducibility
        logger_func: Function for logging
        
    Returns:
        Dictionary with model, performance metrics, and feature importance
    """
    # Split data
    X_train, X_test, y_train, y_test = train_test_split(
        X_diff, y_distances, test_size=test_size, random_state=random_state
    )
    
    logger_func(f"Training data: {X_train.shape[0]} pairs, Test data: {X_test.shape[0]} pairs", level="info")
    
    # Train Random Forest
    logger_func(f"Training Random Forest with {n_estimators} trees", level="info")
    rf = RandomForestRegressor(
        n_estimators=n_estimators,
        max_features=max_features,
        random_state=random_state,
        n_jobs=-1
    )
    
    rf.fit(X_train, y_train)
    
    # Predict on test set
    y_pred = rf.predict(X_test)
    
    # Calculate performance metrics
    mse = mean_squared_error(y_test, y_pred)
    rmse = np.sqrt(mse)
    r2 = r2_score(y_test, y_pred)
    mae = mean_absolute_error(y_test, y_pred)
    
    logger_func(f"Model performance: R² = {r2:.4f}, RMSE = {rmse:.4f}, MAE = {mae:.4f}", level="info")
    
    # Calculate feature importance
    importance = rf.feature_importances_
    
    # Create feature importance DataFrame
    importance_df = pd.DataFrame({
        'Feature': feature_names,
        'Importance': importance
    })
    
    # Sort by importance
    importance_df = importance_df.sort_values('Importance', ascending=False)
    
    return {
        'model': rf,
        'metrics': {
            'mse': mse,
            'rmse': rmse,
            'r2': r2,
            'mae': mae
        },
        'importance': importance_df,
        'predictions': {
            'y_test': y_test,
            'y_pred': y_pred
        }
    }

def plot_feature_importance(importance_df, output_dir, top_n=15, logger_func=log_print):
    """
    Create and save feature importance plots
    
    Args:
        importance_df: DataFrame with feature importance
        output_dir: Directory to save plots
        top_n: Number of top features to include
        logger_func: Function for logging
    """
    # Limit to top_n features
    top_features = importance_df.head(top_n)
    
    # Reverse for better visualization
    top_features = top_features.iloc[::-1]
    
    try:
        # Create figure
        plt.figure(figsize=(12, 8))
        
        # Create horizontal line + dot plot
        plt.hlines(
            y=range(len(top_features)), 
            xmin=0, 
            xmax=top_features["Importance"].values,  
            color="skyblue", 
            alpha=0.7, 
            linewidth=2
        )
        
        plt.plot(
            top_features["Importance"].values,  
            range(len(top_features)), 
            "o", 
            markersize=10, 
            color="blue", 
            alpha=0.8
        )
        
        # Add feature names and values
        plt.yticks(range(len(top_features)), top_features["Feature"].values)
        plt.xlabel("Feature Importance Score")
        plt.title("Variables Associated with Microbiome Composition Differences")
        
        # Add values next to dots
        for i, importance in enumerate(top_features["Importance"].values):  
            plt.text(importance + 0.001, i, f"{importance:.4f}", va='center')
        
        plt.tight_layout()
        
        # Save plot
        plot_file = os.path.join(output_dir, "feature_importance_plot.pdf")
        plt.savefig(plot_file, bbox_inches="tight")
        
        # Also save as PNG
        png_file = os.path.join(output_dir, "feature_importance_plot.png")
        plt.savefig(png_file, dpi=300, bbox_inches="tight")
        
        logger_func(f"Feature importance plot saved to {plot_file} and {png_file}", level="info")
        plt.close()
        
        # Create a regular bar plot as an alternative visualization
        plt.figure(figsize=(10, 6))
        sns.barplot(x="Importance", y="Feature", data=top_features[::-1], palette="viridis")
        plt.xlabel("Feature Importance Score")
        plt.ylabel("Features")
        plt.title("Top Features Driving Microbiome Differences")
        
        # Save bar plot
        bar_plot_file = os.path.join(output_dir, "feature_importance_barplot.pdf")
        plt.savefig(bar_plot_file, bbox_inches="tight")
        plt.close()
        
    except Exception as e:
        logger_func(f"Error creating feature importance plot: {str(e)}", level="error")
        logger_func(traceback.format_exc(), level="error")

def save_model_results(result_dict, output_dir, logger_func=log_print):
    """
    Save model results to files
    
    Args:
        result_dict: Dictionary with model results
        output_dir: Directory to save results
        logger_func: Function for logging
    """
    try:
        # Save feature importance
        importance_file = os.path.join(output_dir, "feature_importance.csv")
        result_dict['importance'].to_csv(importance_file, index=False)
        logger_func(f"Feature importance saved to {importance_file}", level="info")
        
        # Save model summary
        metrics = result_dict['metrics']
        with open(os.path.join(output_dir, "rf_model_summary.txt"), 'w') as f:
            f.write("Random Forest Model Summary\n")
            f.write("==========================\n\n")
            f.write(f"R² score: {metrics['r2']:.4f}\n")
            f.write(f"RMSE: {metrics['rmse']:.4f}\n")
            f.write(f"MAE: {metrics['mae']:.4f}\n\n")
            f.write("Top 10 features:\n")
            for i, row in result_dict['importance'].head(10).iterrows():
                f.write(f"  {row['Feature']}: {row['Importance']:.4f}\n")
        
        # Save predictions
        pred_df = pd.DataFrame({
            'y_true': result_dict['predictions']['y_test'],
            'y_pred': result_dict['predictions']['y_pred']
        })
        pred_file = os.path.join(output_dir, "rf_predictions.csv")
        pred_df.to_csv(pred_file, index=False)
        
    except Exception as e:
        logger_func(f"Error saving model results: {str(e)}", level="error")
        logger_func(traceback.format_exc(), level="error")

def run_feature_selection(abundance_file, metadata_file, output_dir, predictors=None,
                        n_estimators=100, max_features="sqrt", distance_metric="bray",
                        transform="clr", test_size=0.2, random_state=42, log_file=None):
    """
    Run feature selection analysis to identify important variables
    
    Args:
        abundance_file: Path to abundance file
        metadata_file: Path to metadata file
        output_dir: Directory to save output files
        predictors: List of predictor variables to use
        n_estimators: Number of trees in the Random Forest
        max_features: Max features to consider at each split
        distance_metric: Distance metric to use
        transform: Transformation to apply to abundance data
        test_size: Proportion of data for testing
        random_state: Random seed for reproducibility
        log_file: Path to log file
        
    Returns:
        DataFrame with feature importance
    """
    # Setup logging
    if log_file:
        import logging.handlers
        file_handler = logging.handlers.RotatingFileHandler(log_file, maxBytes=10_485_760, backupCount=5)
        file_handler.setFormatter(logging.Formatter('%(asctime)s - %(levelname)s - %(message)s'))
        logger.addHandler(file_handler)
    
    log_print("Starting feature selection analysis", level="info")
    
    # Validate input files
    if not check_file_exists_with_logger(abundance_file, "Abundance file", log_print):
        return None
    
    if not check_file_exists_with_logger(metadata_file, "Metadata file", log_print):
        return None
    
    # Create output directory
    os.makedirs(output_dir, exist_ok=True)
    
    try:
        # Read abundance data
        abundance_df = pd.read_csv(abundance_file, sep='\t', index_col=0)
        log_print(f"Loaded abundance data: {abundance_df.shape[0]} taxa, {abundance_df.shape[1]} samples", level="info")
        
        # Read metadata
        metadata_df = pd.read_csv(metadata_file)
        if 'SampleID' in metadata_df.columns:
            metadata_df.set_index('SampleID', inplace=True)
        log_print(f"Loaded metadata: {metadata_df.shape[0]} samples, {metadata_df.shape[1]} variables", level="info")
        
        # Ensure indices are strings for consistent joining
        abundance_df.index = abundance_df.index.astype(str)
        metadata_df.index = metadata_df.index.astype(str)
        
        # Identify common samples
        common_samples = set(abundance_df.columns) & set(metadata_df.index)
        if len(common_samples) < 10:
            log_print(f"Not enough common samples between abundance and metadata: {len(common_samples)}", level="error")
            return None
        
        log_print(f"Found {len(common_samples)} common samples between abundance and metadata", level="info")
        
        # Subset to common samples
        abundance_subset = abundance_df[list(common_samples)]
        metadata_subset = metadata_df.loc[list(common_samples)]
        
        # Transform abundance data
        abundance_transformed = transform_abundance_data(abundance_subset, transform, log_print)
        
        # Get predictors
        if predictors:
            if isinstance(predictors, str):
                predictors = [p.strip() for p in predictors.split(',')]
            
            # Filter to existing columns
            predictor_cols = [p for p in predictors if p in metadata_subset.columns]
            if not predictor_cols:
                log_print(f"None of the specified predictors exist in metadata: {predictors}", level="error")
                return None
                
            log_print(f"Using specified predictors: {predictor_cols}", level="info")
        else:
            # Use all columns except obviously problematic ones
            exclude_cols = ['sample_id', 'sampleid', 'sample_name', 'samplename']
            predictor_cols = [col for col in metadata_subset.columns if col.lower() not in exclude_cols]
            log_print(f"Using all available metadata columns as predictors: {predictor_cols}", level="info")
        
        # Identify categorical and numerical columns
        categorical_cols = metadata_subset[predictor_cols].select_dtypes(include=['object', 'category']).columns.tolist()
        numerical_cols = [col for col in predictor_cols if col not in categorical_cols]
        
        log_print(f"Categorical predictors: {categorical_cols}", level="info")
        log_print(f"Numerical predictors: {numerical_cols}", level="info")
        
        # Prepare metadata
        encoded_df = encode_categorical_variables(metadata_subset, categorical_cols, log_print)
        
        # Scale numerical variables
        if numerical_cols:
            processed_df = scale_numerical_variables(encoded_df, numerical_cols, log_print)
        else:
            processed_df = encoded_df
        
        # Get final list of features
        feature_names = [col for col in predictor_cols if col in processed_df.columns]
        if not feature_names:
            log_print("No valid features available after processing", level="error")
            return None
            
        log_print(f"Final feature set: {feature_names}", level="info")
        
        # Create pairwise differences
        X_diff, indices = create_pairwise_differences(processed_df[feature_names], log_print)
        
        # Calculate microbiome distances
        y_distances = calculate_microbiome_distances(abundance_transformed.T, distance_metric, log_print)
        if y_distances is None:
            log_print("Failed to calculate microbiome distances", level="error")
            return None
        
        # Train Random Forest model
        result_dict = train_random_forest(
            X_diff, y_distances, feature_names,
            test_size=test_size,
            n_estimators=n_estimators,
            max_features=max_features,
            random_state=random_state,
            logger_func=log_print
        )
        
        # Plot feature importance
        plot_feature_importance(result_dict['importance'], output_dir, logger_func=log_print)
        
        # Save results
        save_model_results(result_dict, output_dir, log_print)
        
        log_print("Feature selection analysis completed successfully", level="info")
        return result_dict['importance']
        
    except Exception as e:
        log_print(f"Error in feature selection analysis: {str(e)}", level="error")
        log_print(traceback.format_exc(), level="error")
        return None

def parse_arguments():
    parser = argparse.ArgumentParser(description="Run feature selection analysis to identify clinical variables associated with microbiome differences")
    
    parser.add_argument(
        "--abundance-file", 
        required=True,
        help="Path to abundance file (filtered_S_abundance.tsv)"
    )
    
    parser.add_argument(
        "--metadata", 
        required=True,
        help="Path to metadata CSV file"
    )
    
    parser.add_argument(
        "--output-dir",
        default="results/feature_selection",
        help="Directory for output files (default: results/feature_selection)"
    )
    
    parser.add_argument(
        "--predictors", 
        default=None,
        help="Comma-separated list of predictor variables to test (default: all available variables)"
    )
    
    parser.add_argument(
        "--distance-metric", 
        default="bray", 
        choices=["bray", "jaccard", "euclidean"],
        help="Distance metric to use (default: bray)"
    )
    
    parser.add_argument(
        "--transform", 
        default="clr", 
        choices=["clr", "hellinger", "log", "none"],
        help="Transformation to apply to abundance data (default: clr)"
    )
    
    parser.add_argument(
        "--n-estimators", 
        type=int, 
        default=100, 
        help="Number of trees in the Random Forest (default: 100)"
    )
    
    parser.add_argument(
        "--test-size", 
        type=float, 
        default=0.2, 
        help="Proportion of data to use for testing (default: 0.2)"
    )
    
    parser.add_argument(
        "--random-state", 
        type=int, 
        default=42, 
        help="Random seed for reproducibility (default: 42)"
    )
    
    parser.add_argument(
        "--log-file",
        default=None,
        help="Path to log file (default: log to console only)"
    )
    
    parser.add_argument(
        "--log-level",
        default="INFO",
        choices=["DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"],
        help="Logging level (default: INFO)"
    )
    
    return parser.parse_args()

def main():
    # Parse command line arguments
    args = parse_arguments()
    
    # Configure logging level
    logger.setLevel(getattr(logging, args.log_level))
    
    log_print("Starting feature selection analysis", level="info")
    
    # Check if files exist
    if not os.path.exists(args.abundance_file):
        log_print(f"Error: Abundance file not found: {args.abundance_file}", level="error")
        sys.exit(1)
    
    if not os.path.exists(args.metadata):
        log_print(f"Error: Metadata file not found: {args.metadata}", level="error")
        sys.exit(1)
    
    # Create output directory
    os.makedirs(args.output_dir, exist_ok=True)
    
    # Parse predictors
    predictors = None
    if args.predictors:
        predictors = args.predictors.split(',')
        log_print(f"Using specified predictors: {', '.join(predictors)}", level="info")
    
    # Run analysis
    result = run_feature_selection(
        abundance_file=args.abundance_file,
        metadata_file=args.metadata,
        output_dir=args.output_dir,
        predictors=predictors,
        n_estimators=args.n_estimators,
        distance_metric=args.distance_metric,
        transform=args.transform,
        test_size=args.test_size,
        random_state=args.random_state,
        log_file=args.log_file
    )
    
    if result is not None:
        log_print("Feature selection completed successfully", level="info")
        log_print(f"Results saved to {args.output_dir}", level="info")
        log_print("\nTop 10 features:", level="info")
        for i, (_, row) in enumerate(result.head(10).iterrows()):
            log_print(f"{i+1}. {row['Feature']}: {row['Importance']:.4f}", level="info")
    else:
        log_print("Feature selection failed", level="error")
        sys.exit(1)

if __name__ == "__main__":
    main()