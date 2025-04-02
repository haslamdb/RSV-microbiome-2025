# kraken_tools/analysis/rf_shap.py
import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import shap
from sklearn.ensemble import RandomForestRegressor
from sklearn.preprocessing import LabelEncoder
from sklearn.model_selection import train_test_split
import statsmodels.api as sm
from statsmodels.formula.api import mixedlm
import logging
import traceback

# Import locally for the standalone version
import sys
sys.path.append(os.path.expanduser("~/Documents/Code/kraken_tools"))
from kraken_tools.utils.file_utils import check_file_exists_with_logger
from kraken_tools.analysis.permanova import transform_abundance_data

def select_target_taxa(abundance_df, target_taxa=None, top_n=10, logger=None):
    """
    Select target taxa for analysis
    
    Args:
        abundance_df: DataFrame with abundance data
        target_taxa: List of specific taxa to analyze
        top_n: Number of top taxa to analyze if target_taxa not specified
        logger: Logger instance
        
    Returns:
        List of target taxa
    """
    if logger is None:
        logger = logging.getLogger('kraken_analysis')
    
    if target_taxa:
        if isinstance(target_taxa, str):
            target_taxa = [t.strip() for t in target_taxa.split(',')]
        
        # Check which specified taxa exist in the data
        existing_taxa = [t for t in target_taxa if t in abundance_df.index]
        if not existing_taxa:
            logger.warning(f"None of the specified target taxa found in data: {target_taxa}")
            logger.info(f"Using top {top_n} most abundant taxa instead")
            # Fall back to top taxa
            mean_abundance = abundance_df.mean(axis=1)
            existing_taxa = mean_abundance.sort_values(ascending=False).head(top_n).index.tolist()
        else:
            if len(existing_taxa) < len(target_taxa):
                logger.warning(f"Some specified taxa not found: {set(target_taxa) - set(existing_taxa)}")
            logger.info(f"Using {len(existing_taxa)} specified taxa: {existing_taxa}")
    else:
        # Use top most abundant taxa
        mean_abundance = abundance_df.mean(axis=1)
        existing_taxa = mean_abundance.sort_values(ascending=False).head(top_n).index.tolist()
        logger.info(f"Using top {len(existing_taxa)} most abundant taxa")
    
    return existing_taxa

def prepare_metadata(metadata_df, predictor_vars, random_effects=None, logger=None):
    """
    Prepare metadata for modeling
    
    Args:
        metadata_df: DataFrame with metadata
        predictor_vars: List of predictor variables
        random_effects: List of random effect variables
        logger: Logger instance
        
    Returns:
        DataFrame with prepared metadata
    """
    if logger is None:
        logger = logging.getLogger('kraken_analysis')
    
    processed_df = metadata_df.copy()
    
    # Handle random effects
    if random_effects:
        if isinstance(random_effects, str):
            random_effects = [r.strip() for r in random_effects.split(',')]
        
        # Check which random effects exist
        existing_random = [r for r in random_effects if r in processed_df.columns]
        if existing_random:
            logger.info(f"Using random effects: {existing_random}")
        else:
            logger.warning(f"None of the specified random effects found in data: {random_effects}")
    else:
        existing_random = []
    
    # Process predictor variables
    if predictor_vars:
        if isinstance(predictor_vars, str):
            predictor_vars = [p.strip() for p in predictor_vars.split(',')]
        
        # Check which predictors exist
        existing_predictors = [p for p in predictor_vars if p in processed_df.columns]
        if not existing_predictors:
            logger.error(f"None of the specified predictors found in data: {predictor_vars}")
            return None
            
        if len(existing_predictors) < len(predictor_vars):
            logger.warning(f"Some predictors not found: {set(predictor_vars) - set(existing_predictors)}")
        
        logger.info(f"Using predictors: {existing_predictors}")
    else:
        # Use all columns except obvious non-predictors
        exclude_cols = ['sample_id', 'sampleid', 'sample_name', 'samplename']
        existing_predictors = [col for col in processed_df.columns if col.lower() not in exclude_cols]
        logger.info(f"Using all available columns as predictors: {existing_predictors}")
    
    # Identify categorical predictors
    categorical_cols = processed_df[existing_predictors].select_dtypes(include=['object', 'category']).columns.tolist()
    
    # Convert categorical variables to categories for formulas
    for col in categorical_cols:
        # Fill NAs
        processed_df[col] = processed_df[col].fillna('missing')
        processed_df[col] = processed_df[col].astype('category')
    
    # For random effects, ensure they're properly formatted
    for col in existing_random:
        if col in processed_df.columns:
            processed_df[col] = processed_df[col].astype(str)
    
    return processed_df, existing_predictors, existing_random

def train_rf_model(X, y, n_estimators=100, random_state=42, logger=None):
    """
    Train Random Forest model for a single taxon
    
    Args:
        X: Feature matrix
        y: Target vector
        n_estimators: Number of trees in the forest
        random_state: Random seed
        logger: Logger instance
        
    Returns:
        Trained model
    """
    if logger is None:
        logger = logging.getLogger('kraken_analysis')
    
    try:
        # Create and train the model
        model = RandomForestRegressor(
            n_estimators=n_estimators,
            random_state=random_state,
            n_jobs=-1
        )
        
        model.fit(X, y)
        
        return model
    
    except Exception as e:
        logger.error(f"Error training Random Forest model: {str(e)}")
        logger.error(traceback.format_exc())
        return None

def calculate_shap_values(model, X, feature_names, logger=None):
    """
    Calculate SHAP values for a trained model
    
    Args:
        model: Trained model
        X: Feature matrix
        feature_names: List of feature names
        logger: Logger instance
        
    Returns:
        SHAP values object
    """
    if logger is None:
        logger = logging.getLogger('kraken_analysis')
    
    try:
        # Create a TreeExplainer for Random Forest
        explainer = shap.TreeExplainer(model)
        
        # Calculate SHAP values
        shap_values = explainer(X)
        
        return shap_values
    
    except Exception as e:
        logger.error(f"Error calculating SHAP values: {str(e)}")
        logger.error(traceback.format_exc())
        return None

def create_shap_plots(shap_values, X, feature_names, taxon, output_dir, logger=None):
    """
    Create and save SHAP plots
    
    Args:
        shap_values: SHAP values object
        X: Feature matrix
        feature_names: List of feature names
        taxon: Name of the taxon
        output_dir: Directory to save plots
        logger: Logger instance
    """
    if logger is None:
        logger = logging.getLogger('kraken_analysis')
    
    try:
        # Create taxon-specific output directory
        taxon_dir = os.path.join(output_dir, f"shap_summary")
        os.makedirs(taxon_dir, exist_ok=True)
        
        # 1. Summary Plot
        plt.figure(figsize=(12, 8))
        shap.summary_plot(shap_values, X, feature_names=feature_names, show=False)
        plt.title(f"SHAP Summary: {taxon}")
        plt.tight_layout()
        plt.savefig(os.path.join(taxon_dir, f"shap_summary_{taxon.replace('.', '_')}.pdf"), bbox_inches="tight")
        plt.savefig(os.path.join(taxon_dir, f"shap_summary_{taxon.replace('.', '_')}.png"), dpi=300, bbox_inches="tight")
        plt.close()
        
        # 2. Feature Importance Plot (Bar Plot)
        plt.figure(figsize=(12, 8))
        shap.plots.bar(shap_values, max_display=15, show=False)
        plt.title(f"SHAP Feature Importance: {taxon}")
        plt.tight_layout()
        plt.savefig(os.path.join(taxon_dir, f"shap_importance_bar_{taxon.replace('.', '_')}.pdf"), bbox_inches="tight")
        plt.close()
        
        # 3. PyCaret-style line+dot plot
        # Calculate mean absolute SHAP values for each feature
        shap_importance = np.abs(shap_values.values).mean(axis=0)
        importance_df = pd.DataFrame({
            "Feature": feature_names,
            "SHAP Importance": shap_importance
        }).sort_values(by="SHAP Importance", ascending=True)
        
        # Create line+dot plot
        plt.figure(figsize=(12, 8))
        top_n = min(15, len(importance_df))
        top_features = importance_df.tail(top_n)
        
        # Create horizontal line + dot plot
        plt.hlines(
            y=range(top_n), 
            xmin=0, 
            xmax=top_features["SHAP Importance"].values,
            color="skyblue", 
            alpha=0.7, 
            linewidth=2
        )
        
        plt.plot(
            top_features["SHAP Importance"].values, 
            range(top_n), 
            "o", 
            markersize=10, 
            color="blue", 
            alpha=0.8
        )
        
        # Add feature names and values
        plt.yticks(range(top_n), top_features["Feature"].values)
        plt.xlabel("SHAP Importance Score")
        plt.title(f"Clinical Variables Associated with {taxon} Abundance")
        
        # Add values next to dots
        for i, importance in enumerate(top_features["SHAP Importance"].values):
            plt.text(importance + 0.001, i, f"{importance:.4f}", va='center')
            
        plt.tight_layout()
        plt.savefig(os.path.join(taxon_dir, f"shap_feature_importance_{taxon.replace('.', '_')}.pdf"), bbox_inches="tight")
        plt.savefig(os.path.join(taxon_dir, f"shap_feature_importance_{taxon.replace('.', '_')}.png"), dpi=300, bbox_inches="tight")
        plt.close()
        
        # 4. Dependence plots for top features (if there are at least 3 features)
        if len(feature_names) >= 3:
            # Get top 3 features by importance
            top3_features = importance_df.tail(3)['Feature'].tolist()
            
            for feature in top3_features:
                plt.figure(figsize=(10, 6))
                feature_idx = feature_names.index(feature)
                shap.dependence_plot(
                    feature_idx, 
                    shap_values.values, 
                    X, 
                    feature_names=feature_names,
                    show=False
                )
                plt.title(f"Dependence Plot: {feature} effect on {taxon}")
                plt.tight_layout()
                plt.savefig(os.path.join(taxon_dir, f"shap_dependence_{taxon.replace('.', '_')}_{feature.replace('.', '_')}.pdf"), bbox_inches="tight")
                plt.close()
        
        return importance_df
    
    except Exception as e:
        logger.error(f"Error creating SHAP plots: {str(e)}")
        logger.error(traceback.format_exc())
        return None

def fit_mixed_model(taxon, abundance, metadata_df, predictors, random_effects, output_dir, logger=None):
    """
    Fit a mixed effects model for a single taxon
    
    Args:
        taxon: Name of the taxon
        abundance: Abundance vector for the taxon
        metadata_df: DataFrame with metadata
        predictors: List of predictor variables
        random_effects: List of random effect variables
        output_dir: Directory to save results
        logger: Logger instance
        
    Returns:
        Dictionary with model results
    """
    if logger is None:
        logger = logging.getLogger('kraken_analysis')
    
    try:
        # Create data for modeling
        model_data = metadata_df.copy()
        model_data['MicrobeAbundance'] = abundance
        
        # Create mixed model directory
        model_dir = os.path.join(output_dir, "mixed_models")
        os.makedirs(model_dir, exist_ok=True)
        
        # Create formula
        # For categorical variables, use C(var)
        formula_parts = []
        for pred in predictors:
            if pred in model_data.columns and pred not in random_effects:
                if model_data[pred].dtype.name == 'category':
                    formula_parts.append(f"C({pred})")
                else:
                    formula_parts.append(pred)
        
        if not formula_parts:
            logger.warning(f"No valid predictors for {taxon}, skipping mixed model")
            return None
            
        fixed_effects = " + ".join(formula_parts)
        
        # Add random effects if any
        random_formula = None
        if random_effects:
            # For now, we only support simple random intercepts
            random_parts = []
            for re in random_effects:
                if re in model_data.columns:
                    # Check if there are enough groups
                    n_groups = model_data[re].nunique()
                    if n_groups >= 5:  # Minimum for meaningful random effects
                        random_parts.append(f"(1|{re})")
                    else:
                        logger.warning(f"Not enough groups for random effect {re} (found {n_groups}), skipping")
            
            if random_parts:
                random_formula = " + ".join(random_parts)
        
        # Full formula
        if random_formula:
            formula = f"MicrobeAbundance ~ {fixed_effects} + {random_formula}"
        else:
            formula = f"MicrobeAbundance ~ {fixed_effects}"
            
        logger.info(f"Fitting mixed model for {taxon} with formula: {formula}")
        
        # Fit the model
        result = None
        model_type = None
        
        try:
            # First try mixed model if we have random effects
            if random_formula:
                for re in random_effects:
                    if re in model_data.columns:
                        # Try with mixed model
                        try:
                            model = mixedlm(f"MicrobeAbundance ~ {fixed_effects}", model_data, groups=model_data[re])
                            result = model.fit(maxiter=1000)
                            model_type = "MixedLM"
                            break
                        except:
                            logger.warning(f"Mixed model failed with random effect {re}, trying another")
                            continue
            
            # Fall back to OLS if no random effects or mixed model failed
            if result is None:
                logger.info(f"Fitting OLS model for {taxon}")
                model = sm.formula.ols(f"MicrobeAbundance ~ {fixed_effects}", data=model_data)
                result = model.fit()
                model_type = "OLS"
        
        except Exception as e:
            logger.error(f"Error fitting model for {taxon}: {str(e)}")
            return None
        
        # Save model results
        with open(os.path.join(model_dir, f"model_results_{taxon.replace('.', '_')}.txt"), "w") as f:
            f.write(str(result.summary()))
            f.write(f"\nModel type: {model_type}")
        
        # Extract coefficients
        coef_df = pd.DataFrame({
            'Variable': result.params.index.tolist(),
            'Coefficient': result.params.values,
            'P_value': result.pvalues.values if hasattr(result, 'pvalues') else np.nan,
            'Model_Type': model_type,
            'Taxon': taxon
        })
        
        # Save coefficients
        coef_df.to_csv(os.path.join(model_dir, f"coefficients_{taxon.replace('.', '_')}.csv"), index=False)
        
        # Create coefficient plot
        plt.figure(figsize=(10, 6))
        # Filter out intercept and random effects
        plot_coefs = coef_df[~coef_df['Variable'].str.contains('Intercept|groups')].copy()
        if not plot_coefs.empty:
            # Sort by absolute coefficient value
            plot_coefs['AbsCoef'] = np.abs(plot_coefs['Coefficient'])
            plot_coefs = plot_coefs.sort_values('AbsCoef')
            
            # Plot coefficients
            bars = plt.barh(
                plot_coefs['Variable'], 
                plot_coefs['Coefficient'],
                color=[
                    'green' if (c > 0 and p < 0.05) else 
                    'red' if (c < 0 and p < 0.05) else 
                    'lightgreen' if c > 0 else 'lightcoral' 
                    for c, p in zip(plot_coefs['Coefficient'], plot_coefs['P_value'])
                ],
                alpha=0.7
            )
            
            # Add p-value annotations
            for i, (_, row) in enumerate(plot_coefs.iterrows()):
                if row['P_value'] < 0.001:
                    sig = "***"
                elif row['P_value'] < 0.01:
                    sig = "**"
                elif row['P_value'] < 0.05:
                    sig = "*"
                else:
                    sig = ""
                
                if sig:
                    offset = 0.1 if row['Coefficient'] >= 0 else -0.3
                    plt.text(row['Coefficient'] + offset, i, sig, va='center')
            
            plt.axvline(0, color='gray', linestyle='--')
            plt.xlabel('Coefficient')
            plt.title(f'Mixed Model Coefficients: {taxon}')
            plt.tight_layout()
            plt.savefig(os.path.join(model_dir, f"coefficient_plot_{taxon.replace('.', '_')}.pdf"), bbox_inches="tight")
            plt.savefig(os.path.join(model_dir, f"coefficient_plot_{taxon.replace('.', '_')}.png"), dpi=300, bbox_inches="tight")
            plt.close()
        
        return {
            'coefficients': coef_df,
            'model_type': model_type,
            'formula': formula
        }
    
    except Exception as e:
        logger.error(f"Error in mixed model analysis for {taxon}: {str(e)}")
        logger.error(traceback.format_exc())
        return None

def run_rf_shap_analysis(abundance_file, metadata_file, output_dir, target_taxa=None,
                       predictors=None, random_effects=None, transform="clr",
                       n_estimators=100, test_size=0.2, top_n=10, 
                       mixed_model="lmer", log_file=None):
    """
    Run Random Forest with SHAP analysis for target taxa
    
    Args:
        abundance_file: Path to abundance file
        metadata_file: Path to metadata file
        output_dir: Directory to save output files
        target_taxa: List of specific taxa to analyze
        predictors: List of predictor variables
        random_effects: Variables to treat as random effects
        transform: Transformation to apply to abundance data
        n_estimators: Number of trees in the Random Forest
        test_size: Proportion of data to use for testing
        top_n: Number of top taxa to analyze if target_taxa not specified
        mixed_model: Type of mixed model ("lmer" or "glmm")
        log_file: Path to log file
        
    Returns:
        Dictionary with results for each taxon
    """
    # Setup logging
    import logging as logging_module  # Fixed import
    
    if log_file:
        logger = logging_module.getLogger('kraken_analysis')
        handler = logging_module.handlers.RotatingFileHandler(log_file, maxBytes=10_485_760, backupCount=5)
        handler.setFormatter(logging_module.Formatter('%(asctime)s - %(levelname)s - %(message)s'))
        logger.addHandler(handler)
    else:
        logger = logging_module.getLogger('kraken_analysis')
    
    logger.info("Starting Random Forest with SHAP analysis")
    
    # Validate input files
    if not check_file_exists_with_logger(abundance_file, "Abundance file", logger):
        return None
    
    if not check_file_exists_with_logger(metadata_file, "Metadata file", logger):
        return None
    
    # Create output directory
    os.makedirs(output_dir, exist_ok=True)
    
    try:
        # Read abundance data
        abundance_df = pd.read_csv(abundance_file, sep='\t', index_col=0)
        logger.info(f"Loaded abundance data: {abundance_df.shape[0]} taxa, {abundance_df.shape[1]} samples")
        
        # Read metadata
        metadata_df = pd.read_csv(metadata_file, index_col=0)
        logger.info(f"Loaded metadata: {metadata_df.shape[0]} samples, {metadata_df.shape[1]} variables")
        
        # Ensure indices are strings for consistent joining
        abundance_df.index = abundance_df.index.astype(str)
        metadata_df.index = metadata_df.index.astype(str)
        
        # Identify common samples
        common_samples = list(set(abundance_df.columns) & set(metadata_df.index))
        if len(common_samples) < 10:
            logger.error(f"Not enough common samples between abundance and metadata: {len(common_samples)}")
            return None
        
        logger.info(f"Found {len(common_samples)} common samples between abundance and metadata")
        
        # Subset to common samples
        abundance_subset = abundance_df[common_samples]
        metadata_subset = metadata_df.loc[common_samples]
        
        # Transform abundance data
        abundance_transformed = transform_abundance_data(abundance_subset, transform, logger)
        
        # Select target taxa
        selected_taxa = select_target_taxa(abundance_transformed, target_taxa, top_n, logger)
        if not selected_taxa:
            logger.error("No taxa selected for analysis")
            return None
        
        # Prepare metadata
        prepared_data = prepare_metadata(metadata_subset, predictors, random_effects, logger)
        if prepared_data is None:
            logger.error("Failed to prepare metadata")
            return None
        
        processed_df, existing_predictors, existing_random = prepared_data
        
        # Store results for all taxa
        all_results = {}
        all_shap_importance = []
        all_model_results = []
        
        # Process each taxon
        for taxon in selected_taxa:
            logger.info(f"\nAnalyzing taxon: {taxon}")
            
            # Extract abundance for this taxon
            taxon_abundance = abundance_transformed.loc[taxon]
            
            # Create metadata for model training (encode categorical variables)
            X_data = processed_df.copy()
            
            # Identify categorical predictors
            categorical_cols = [col for col in existing_predictors 
                               if col in X_data.columns and X_data[col].dtype.name == 'category']
            
            # Label encode for Random Forest
            X_encoded = X_data.copy()
            for col in categorical_cols:
                try:
                    X_encoded[col] = LabelEncoder().fit_transform(X_encoded[col])
                except:
                    logger.warning(f"Error encoding {col}, dropping column")
                    X_encoded = X_encoded.drop(columns=[col])
            
            # Get final feature list
            feature_names = [col for col in existing_predictors if col in X_encoded.columns]
            if not feature_names:
                logger.warning(f"No valid features for {taxon} after encoding, skipping")
                continue
                
            # Split data for model training
            X = X_encoded[feature_names].values
            y = taxon_abundance.values
            
            X_train, X_test, y_train, y_test = train_test_split(
                X, y, test_size=test_size, random_state=42
            )
            
            # Train Random Forest model
            model = train_rf_model(X_train, y_train, n_estimators, 42, logger)
            if model is None:
                logger.warning(f"Failed to train Random Forest for {taxon}, skipping")
                continue
                
            # Calculate SHAP values
            shap_values = calculate_shap_values(model, X_test, feature_names, logger)
            if shap_values is None:
                logger.warning(f"Failed to calculate SHAP values for {taxon}, skipping")
                continue
                
            # Create SHAP plots
            importance_df = create_shap_plots(shap_values, X_test, feature_names, taxon, output_dir, logger)
            if importance_df is not None:
                # Add taxon to importance data
                importance_df['Taxon'] = taxon
                all_shap_importance.append(importance_df)
            
            # Fit mixed model
            model_results = fit_mixed_model(
                taxon, 
                taxon_abundance, 
                processed_df, 
                existing_predictors, 
                existing_random, 
                output_dir, 
                logger
            )
            
            if model_results is not None:
                all_model_results.append(model_results['coefficients'])
            
            # Store results for this taxon
            all_results[taxon] = {
                'shap_importance': importance_df,
                'model_results': model_results,
                'top_predictors': importance_df['Feature'].values[-5:].tolist() if importance_df is not None else []
            }
            
            logger.info(f"Completed analysis for {taxon}")
        
        # Combine all SHAP results
        if all_shap_importance:
            combined_shap = pd.concat(all_shap_importance)
            combined_shap.to_csv(os.path.join(output_dir, "shap_feature_importance_all_microbes.csv"), index=False)
            logger.info(f"Saved combined SHAP importance to {os.path.join(output_dir, 'shap_feature_importance_all_microbes.csv')}")
        
        # Combine all model results
        if all_model_results:
            combined_models = pd.concat(all_model_results)
            combined_models.to_csv(os.path.join(output_dir, "model_summary_all_microbes.csv"), index=False)
            logger.info(f"Saved combined model results to {os.path.join(output_dir, 'model_summary_all_microbes.csv')}")
        
        logger.info("RF-SHAP analysis completed successfully")
        return all_results
        
    except Exception as e:
        logger.error(f"Error in RF-SHAP analysis: {str(e)}")
        logger.error(traceback.format_exc())
        return None