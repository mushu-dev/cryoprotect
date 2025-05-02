# Reports Organization Plan

## Current Situation

Verification reports are scattered across the repository with timestamped filenames:

### RLS Verification Reports
```
./rls_verification_report_20250418_131730.json
./rls_verification_report_20250418_131904.json
./rls_verification_report_20250418_132127.json
./rls_verification_report_20250418_132704.json
./rls_verification_report_20250418_132744.json
```

### RLS Verification Results
```
./rls_verification_results_20250418_115950.json
./rls_verification_results_20250418_132009.json
./rls_verification_results_20250418_132051.json
```

### API Verification Files
```
./reports/api_verification_standalone_20250417_163050.json
./reports/api_verification_standalone_20250417_170311.json
./reports/api_verification_standalone_20250417_171622.json
./reports/api_verification_standalone_20250417_193217.json
./reports/api_verification_standalone_20250417_211621.json
```

## Recommended Structure

Create a structured `reports/` directory with subdirectories by type:

```
reports/
├── api/
│   └── verification/
│       └── latest.json  (copy of most recent report)
├── rls/
│   ├── reports/
│   │   └── latest.json  (copy of most recent report)
│   └── results/
│       └── latest.json  (copy of most recent report)
└── archives/
    ├── api/
    ├── rls/
    └── README.md  (explain archiving policy)
```

## Immediate Actions

1. Create the directory structure
2. Move the most recent reports to the appropriate locations
3. Move older reports to archives
4. Add reports directory to `.gitignore` to prevent future commit of generated reports

## Commands

```bash
# Create directory structure
mkdir -p reports/api/verification
mkdir -p reports/rls/reports
mkdir -p reports/rls/results
mkdir -p reports/archives/api
mkdir -p reports/archives/rls

# Move most recent reports (assuming latest has highest timestamp)
cp ./reports/api_verification_standalone_20250417_211621.json ./reports/api/verification/latest.json
cp ./rls_verification_report_20250418_132744.json ./reports/rls/reports/latest.json
cp ./rls_verification_results_20250418_132051.json ./reports/rls/results/latest.json

# Archive older reports
mv ./reports/api_verification_standalone_*.json ./reports/archives/api/
mv ./rls_verification_report_*.json ./reports/archives/rls/
mv ./rls_verification_results_*.json ./reports/archives/rls/

# Add to .gitignore
echo "# Ignore report files except latest" >> .gitignore
echo "reports/archives/" >> .gitignore
echo "reports/**/*.json" >> .gitignore
echo "!reports/**/latest.json" >> .gitignore

# Create README for archives
echo "# Archives\n\nThis directory contains archived reports that should generally not be committed to version control." > reports/archives/README.md

# Commit changes
git add reports/
git commit -m "Organize report files into structured directory"
```