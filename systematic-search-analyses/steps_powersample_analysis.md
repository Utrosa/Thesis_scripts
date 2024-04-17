# Power and meta-analysis steps

## Prepare the data
1. Create an Excel spreadsheet that systemizes information about all records selected in the systematic literature search. This information should include at least:
- record identifiers (e.g.: Lan_2023),
- sample size of autistic and control participants,
- mean age of the sample,
- type of effect size measure (e.g.: Cohen's d, correlation coefficient, ...), 
- effect size estimate, and
- significance level (`&alpha`).
2. Check that the spreadsheet has a default value of `1.0` for all empty cells in the `reliability` column.
3. Sort the spreadsheet by age category: infant, toddler, preschool children, school children, teen, young adult, adult. 
4. Save the file as: `Meta-Analysis-Static-FC.xlsx`.
5. Create a copy of the Excel sheet `data`. Save it as a Tab Delimited file: `Meta-Analysis-Static-FC.tsv`.

## Run the analysis
- Run Python notebook on `Meta-Analysis-Static-FC.xlsx`.
- Run R code on `Meta-Analysis-Static-FC.tsv`.