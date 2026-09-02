# SurvAssist

**SurvAssist** is an interactive R Shiny web application designed to enable
clinicians and biomedical researchers to perform survival analysis without
programming expertise.

## Web application

SurvAssist is freely accessible through a web browser:

https://survassist.shinyapps.io/SurvAssist/

## Features

SurvAssist provides an integrated workflow for:

- Data import from CSV and XLSX files
- Descriptive statistics
- Kaplan-Meier survival analysis
- Log-rank testing
- Survival probability estimation
- Median survival estimation
- Follow-up adequacy assessment
- Univariable Cox proportional hazards regression
- Multivariable Cox proportional hazards regression
- Kaplan-Meier plot customization
- Downloadable statistical results
- Downloadable publication-ready plots
- Automated Microsoft Word report generation

## Intended users

SurvAssist is primarily intended for:

- Clinicians
- Biomedical researchers
- Epidemiologists
- Medical researchers
- Students learning survival analysis

No programming knowledge is required to use the web application.

## Input data

The application accepts datasets in:

- CSV format
- Microsoft Excel XLSX format

Users should provide a survival time variable and an event/status variable.
Additional grouping variables and covariates can be selected depending on
the analysis.

## Example dataset

The `sample_data/data.xlsx` file contains the publicly available
NCCTG lung cancer dataset distributed with the R `survival` package.

The dataset can be used to test the SurvAssist workflow. For example:

- `time` or `time_in_months` can be selected as the survival time variable.
- `status` can be selected as the event/status variable.
- Additional clinical variables can be selected as covariates for Kaplan-Meier analysis or Cox
  proportional hazards regression.

The dataset is included only as an example for demonstrating and testing
the application.
## Statistical methods

SurvAssist implements commonly used survival analysis methods using the
R statistical environment, including Kaplan-Meier estimation, log-rank
testing, reverse Kaplan-Meier follow-up assessment, and Cox proportional
hazards regression.

## Software architecture

SurvAssist was developed using R and the Shiny web application framework.
Statistical analysis and visualization are implemented using established
R packages including `survival` and `survminer`.

## Installation

SurvAssist can be used directly through the web application without local
installation.

For local deployment, download or clone this repository and run `app.R`
using R/RStudio with the required R packages installed.

## Reproducibility

A sample dataset is provided in the repository to allow users to test the
application and reproduce the example analyses described in the associated
publication.

## Citation

If you use SurvAssist in your research, please cite the software using the
citation information provided in `CITATION.cff` and the associated software
DOI.

## License

SurvAssist is distributed under the MIT License. See the `LICENSE` file
for details.

## Support

For questions, bug reports, or feature requests, please use the GitHub
Issues section of this repository.
