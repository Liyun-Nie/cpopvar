# cpopvar Chloroplast Genome Variation Report

**Project**: {{PROJECT_NAME}}  
**Session ID**: {{SESSION_ID}}  
**Analysis date**: {{REPORT_DATE}}  
**Duration**: {{ANALYSIS_DURATION}}  
**R version**: {{R_VERSION}}  

---

## Summary

This report collects results from the modules enabled in the cpopvar configuration.

## Analysis flowchart

{{ANALYSIS_FLOWCHART}}

### Data statistics
- **Raw records**: {{RAW_DATA_COUNT}} variant records
- **Processed records**: {{PROCESSED_DATA_COUNT}} retained variant records
- **Species analysed**: {{SPECIES_COUNT}}
- **Genome regions**: {{REGION_TYPES}}

---

## Result files

{{MODULE_RESULTS}}

---

## Technical notes

### Processing steps
1. **Preprocessing**: region annotation, IR coordinate conversion, and deduplication
2. **Variant filtering**: frequency-threshold filtering
3. **Frequency normalisation**: variants per kilobase
4. **Visualisation**: M01/M02/M03 module analyses

### Configuration
- **Execution mode**: {{EXECUTION_MODE}}
- **Session ID**: `{{SESSION_ID}}`
- **Output location**: determined jointly by runtime `output_dir` and the session directory

### Data quality
- Input, configuration, and output contracts are those of the installed version
- Statistical limitations and scope should be interpreted for the enabled modules

---

**Generated**: {{GENERATION_TIME}}  
**Pipeline**: cpopvar  
**System**: {{SYSTEM_INFO}}
