# RSV Microbiome Feature Selection Analysis

## Overview
This analysis identifies which clinical variables are most associated with differences in microbiome composition in the RSV dataset.

## Methodology
- We used Random Forest to predict pairwise microbiome distances from pairwise differences in clinical variables
- The importance scores indicate which clinical variables best explain microbiome composition differences
- Higher importance scores (closer to 1.0) indicate stronger relationships between the variable and microbiome composition

## Results

### Method 1: CLR Transformation with Euclidean Distance
| Feature  | Importance |
|----------|------------|
| Timing   | 0.5324     |
| Severity | 0.3466     |
| Symptoms | 0.1211     |

Model Performance:
- R² score: 0.0050
- RMSE: 486.1728
- MAE: 378.9784

### Method 2: Hellinger Transformation with Bray-Curtis Distance
| Feature  | Importance |
|----------|------------|
| Severity | 0.3902     |
| Timing   | 0.3750     |
| Symptoms | 0.2348     |

Model Performance:
- R² score: 0.0009
- RMSE: 0.1668
- MAE: 0.1333

## Key Findings

1. **Timing** was identified as the most important clinical variable when using CLR transformation (53.2% of importance), suggesting that the disease stage (prior, acute, post) is strongly associated with microbiome differences.

2. **Severity** was identified as the most important variable when using Hellinger transformation (39.0% of importance), followed closely by Timing (37.5%).

3. **Symptoms** consistently ranked as the least important variable across both methods, but still captured 12.1% - 23.5% of the total importance.

4. The low R² values (0.0009 - 0.0050) suggest that while these clinical variables do help explain microbiome composition differences, there are likely many other factors influencing microbiome variation that are not captured in our current set of clinical variables.

## Conclusion
The timing of sample collection relative to RSV infection and the severity of illness appear to be the most strongly associated with differences in microbiome composition. This aligns with our previous PERMANOVA findings, which also showed Timing and Severity as significant factors. The biological interpretation could be that different stages of infection and different severity levels create distinct environments in the respiratory tract that select for different microbial communities.