# {{data_source}} Data Quality Report

Generated on: {{timestamp}}

## Executive Summary

Overall data quality assessment: **{{overall_assessment}}**

{{summary_text}}

### Database Statistics

| Metric | Count |
| ------ | ----- |
{{#counts}}
| {{name}} | {{value}} |
{{/counts}}

## Detailed Assessments

{{#assessments}}
### {{name}}: {{status}}

{{details}}

{{#has_visualization}}
![{{name}} Visualization]({{visualization_path}})
{{/has_visualization}}
{{/assessments}}

## Property Completeness

Average property coverage: **{{average_coverage}}%**

| Property | Coverage | Molecules with Property | Total Molecules |
| -------- | -------- | ----------------------- | --------------- |
{{#properties}}
| {{name}} | {{coverage}}% | {{molecules_with_property}} | {{total_molecules}} |
{{/properties}}

## Reference Compounds

- Total reference compounds: {{total_reference_compounds}}
- Present in database: {{present_reference_compounds}}
- Missing from database: {{missing_reference_compounds}}
- Success rate: {{reference_success_rate}}%

{{#has_missing_compounds}}
### Missing Reference Compounds

```
{{missing_compounds_list}}
```
{{/has_missing_compounds}}

## Recommendations

{{#recommendations}}
- {{priority}}: {{text}}
{{#sub_recommendations}}
  - {{text}}
{{/sub_recommendations}}
{{/recommendations}}

## Next Steps

{{next_steps}}