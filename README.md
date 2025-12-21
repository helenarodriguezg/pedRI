# pedRI

**pedRI** offers tools for calculating pediatric reference intervals adjusted for age and inflammation, facilitating accurate interpretation. Serum copper is used as a proof-of-concept.

This package **does not require advanced programming skills**.  

---

## Table of Contents

1. [Installation](#installation)  
2. [Getting Started](#getting-started)  
3. [Step-by-Step Usage](#step-by-step-usage)  
    - [1. Load the package](#1-load-the-package)  
    - [2. Load your data](#2-load-your-data)  
    - [3. Filter samples without inflammation](#3-filter-samples-without-inflammation)  
    - [4. Fit age-dependent copper model](#4-fit-age-dependent-copper-model)  
    - [5. Fit reference intervals](#5-fit-reference-intervals)  
    - [6. Fit inflammation score](#6-fit-inflammation-score)  
    - [7. Build the complete copper model](#7-build-the-complete-copper-model)  
    - [8. Interpret individual patient results](#8-interpret-individual-patient-results)  
4. [Example Dataset](#example-dataset)  

---

## Installation

You can install `pedRI` directly from GitHub:

```R
# Install devtools if not already installed
install.packages("devtools")

# Install pedRI from GitHub
devtools::install_github("helenarodriguezg/pedRI")
````

Load the package:

```R
library(pedRI)
```

---

## Getting Started

For the serum copper proof-of-concept, you will need a dataset containing:

* **Copper concentrations (µg/L)**
* **Patient age (months)**
* **Inflammation markers:** CRP (mg/L), ESR (mm/h), Fibrinogen (g/L)

An example dataset is included at `example_data\data.xlsx`.

---

## Step-by-Step Usage

### 1. Load the package

```R
library(pedRI)
```

### 2. Load your data

```R
read.csv2("data.csv")
head(data)
```

Make sure your column names match:

* Copper: `"Cu"`
* Age: `"age.months"`
* CRP: `"CRP"`
* ESR: `"ESR"`
* Fibrinogen: `"fibrinogen"`

---

### 3. Filter samples without inflammation

```R
data$infl_flag <- inflammation_flag(data)
normal_data <- subset(data, infl_flag == FALSE)
```

This ensures the age model is fitted only on samples **without systemic inflammation**.

---

### 4. Fit age-dependent copper model

```R
age_mod <- fit_age_model(normal_data, cu_var = "Cu", age_var = "age.months")
```

This models copper levels as a function of age.

---

### 5. Fit reference intervals

```R
ri <- fit_reference_intervals(age_mod, data = normal_data, n_boot = 1000)
```
This provides coefficients to compute lower and upper reference limits.
---

### 6. Fit inflammation score

```R
infl_score <- fit_inflammation_score(age_mod, data, 
                                     fibrinogen_col = "fibrinogen",
                                     crp_col = "CRP",
                                     esr_col = "ESR")
infl_score$B3        # coefficient for inflammation adjustment
infl_score$weights   # weights for each marker
```

This calculates a combined inflammation score for adjusting copper levels.

---

### 7. Build the complete copper model

```R
copper_mod <- build_copper_model(age_mod, ri, infl_score)
print(copper_mod)
```

This object combines the age model, reference intervals, and inflammation adjustment.

---

### 8. Interpret individual patient results

```R
result <- interpret_copper(
  cu = 950,           # patient's copper (µg/L)
  age_months = 24,    # patient's age
  crp = 12,           # CRP value
  esr = 25,           # ESR value
  fibrinogen = 3.5,   # Fibrinogen value
  model = copper_mod
)

print(result)
```

**Output explanation:**

* `Cu (inflammation adjusted)`: copper adjusted for age and inflammation
* `Reference values`: age-adjusted lower and upper limits
* `Interpretation`: `LOW`, `NORMAL`, or `HIGH`
* `Inflammation?`: whether any marker indicates inflammation

---

## Example Dataset

An example dataset `data.xlsx` is provided with columns:

| Cu  | age.months | CRP | ESR | fibrinogen |
| --- | ---------- | --- | --- | ---------- |
| 950 | 24         | 12  | 25  | 3.5        |
| ... | ...        | ... | ... | ...        |

---

